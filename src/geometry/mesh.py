import numpy as np
import warp as wp
import meshio
import os
from collections import defaultdict

class Mesh:
    """
    Manages the computational mesh (topology and geometry).

    This class handles the creation or loading of the mesh, including vertices, connectivity,
    and geometric metrics required for the DG solver. It supports both structured Cartesian meshes
    (generated internally) and unstructured quadrilateral meshes loaded via meshio (e.g., from Gmsh).

    Attributes:
        nx (int): Number of elements in x-direction (Cartesian only).
        ny (int): Number of elements in y-direction (Cartesian only).
        num_elements (int): Total number of elements in the mesh.
        device (str): Warp compute device.
        x_min (float): Domain minimum x-coordinate.
        x_max (float): Domain maximum x-coordinate.
        y_min (float): Domain minimum y-coordinate.
        y_max (float): Domain maximum y-coordinate.
        dx (float): Characteristic element size (used for CFL).
        vertices (wp.array): Vertex coordinates (NumElements, 4, 2).
        connectivity (wp.array): Neighbor connectivity (NumElements, 4).
        boundary_tags (wp.array): Physical boundary tags per face (NumElements, 4).
        rx, ry, sx, sy (wp.array): Inverse Jacobian matrix components.
        J (wp.array): Jacobian determinant.
        Js (wp.array): Surface Jacobian (scaling factors) for faces.
        physical_groups (dict): Mapping from boundary name to integer tag.
    """
    def __init__(self, nx=None, ny=None, x_min=0.0, x_max=1.0, y_min=0.0, y_max=1.0, filename=None, device="cuda"):
        """
        Initializes the Mesh.

        Args:
            nx (int, optional): Number of elements in X (for Cartesian mesh).
            ny (int, optional): Number of elements in Y (for Cartesian mesh).
            x_min (float, optional): Min X coordinate. Defaults to 0.0.
            x_max (float, optional): Max X coordinate. Defaults to 1.0.
            y_min (float, optional): Min Y coordinate. Defaults to 0.0.
            y_max (float, optional): Max Y coordinate. Defaults to 1.0.
            filename (str, optional): Path to a mesh file to load. Overrides nx/ny.
            device (str, optional): Compute device ("cpu" or "cuda"). Defaults to "cuda".

        Raises:
            ValueError: If neither filename nor (nx, ny) are provided.
        """
        self.device = device
        
        if filename and os.path.exists(filename):
            print(f"Loading mesh from {filename}...")
            self._load_from_file(filename)
        elif nx is not None and ny is not None:
            self.nx = nx
            self.ny = ny
            self.num_elements = nx * ny
            
            self.x_min = x_min
            self.x_max = x_max
            self.y_min = y_min
            self.y_max = y_max

            self.dx = (x_max - x_min) / nx
            self.dy = (y_max - y_min) / ny

            # --- Mesh Data (Host) ---
            self.vertices_host = np.zeros((self.num_elements, 4, 2), dtype=np.float32)
            self.connectivity_host = -np.ones((self.num_elements, 4, 2), dtype=np.int32)
            
            self._create_cartesian_mesh()
            self._create_cartesian_connectivity()
        else:
            raise ValueError("Either filename or (nx, ny) must be provided.")

        # --- Geometric Factors Arrays (initialized to 0, computed by Solver) ---
        # We prepare these arrays here so they are available as member variables.
        # Unlike the Cartesian case, we don't preset them with constant values here
        # because for unstructured grids they vary element-by-element.
        
        # Inverse Jacobian Matrix components: 
        # [dr/dx  dr/dy]
        # [ds/dx  ds/dy]
        self.rx_host = np.zeros(self.num_elements, dtype=np.float32)
        self.ry_host = np.zeros(self.num_elements, dtype=np.float32)
        self.sx_host = np.zeros(self.num_elements, dtype=np.float32)
        self.sy_host = np.zeros(self.num_elements, dtype=np.float32)
        
        # Determinant of Jacobian (dx/dr * dy/ds - ...)
        self.J_host  = np.zeros(self.num_elements, dtype=np.float32)
        
        # Surface Jacobians (scaling factors for face integration)
        # Store as (NumElements, 4) since they can vary per face in unstructured
        self.Js_host = np.zeros((self.num_elements, 4), dtype=np.float32)

        # --- Transfer to Warp ---
        self.vertices = wp.array(self.vertices_host, dtype=wp.vec2, device=self.device)
        
        # Connectivity: Solver expects (NumElements, 4) -> NeighborID
        self.connectivity = wp.array(self.connectivity_host[:, :, 0], dtype=wp.int32, device=self.device)
        self.connectivity_face_indices = wp.array(self.connectivity_host[:, :, 1], dtype=wp.int32, device=self.device)
        
        # Boundary Tags: (NumElements, 4). 0 = Internal, >0 = Physical Tag
        if not hasattr(self, 'boundary_tags_host'):
             self.boundary_tags_host = np.zeros((self.num_elements, 4), dtype=np.int32)
             
        self.boundary_tags = wp.array(self.boundary_tags_host, dtype=wp.int32, device=self.device)
        
        self.rx = wp.array(self.rx_host, dtype=wp.float32, device=self.device)
        self.ry = wp.array(self.ry_host, dtype=wp.float32, device=self.device)
        self.sx = wp.array(self.sx_host, dtype=wp.float32, device=self.device)
        self.sy = wp.array(self.sy_host, dtype=wp.float32, device=self.device)
        self.J  = wp.array(self.J_host, dtype=wp.float32, device=self.device)
        
        # New Surface Jacobian array
        self.Js = wp.array(self.Js_host, dtype=wp.float32, device=self.device)

    def _create_cartesian_mesh(self):
        """Generates vertices for a structured Cartesian grid."""
        for j in range(self.ny):
            for i in range(self.nx):
                element_id = j * self.nx + i
                x0 = self.x_min + i * self.dx
                y0 = self.y_min + j * self.dy
                
                # Counter-Clockwise ordering
                self.vertices_host[element_id, 0, :] = [x0, y0]
                self.vertices_host[element_id, 1, :] = [x0 + self.dx, y0]
                self.vertices_host[element_id, 2, :] = [x0 + self.dx, y0 + self.dy]
                self.vertices_host[element_id, 3, :] = [x0, y0 + self.dy]

    def _create_cartesian_connectivity(self):
        """
        Generates connectivity and boundary tags for a structured Cartesian grid.
        
        Faces are ordered: 0: Bottom, 1: Right, 2: Top, 3: Left.
        Assigns standard physical tags: 1=Bottom, 2=Right, 3=Top, 4=Left.
        """
        self.boundary_tags_host = np.zeros((self.num_elements, 4), dtype=np.int32)
        
        for j in range(self.ny):
            for i in range(self.nx):
                element_id = j * self.nx + i

                # Left (Face 3)
                if i > 0:
                    ni, nj = i - 1, j
                    self.connectivity_host[element_id, 3] = [nj * self.nx + ni, 1]
                else:
                    self.boundary_tags_host[element_id, 3] = 4 # Left Tag

                # Right (Face 1)
                if i < self.nx - 1:
                    ni, nj = i + 1, j
                    self.connectivity_host[element_id, 1] = [nj * self.nx + ni, 3]
                else:
                    self.boundary_tags_host[element_id, 1] = 2 # Right Tag

                # Bottom (Face 0)
                if j > 0:
                    ni, nj = i, j - 1
                    self.connectivity_host[element_id, 0] = [nj * self.nx + ni, 2]
                else:
                    self.boundary_tags_host[element_id, 0] = 1 # Bottom Tag

                # Top (Face 2)
                if j < self.ny - 1:
                    ni, nj = i, j + 1
                    self.connectivity_host[element_id, 2] = [nj * self.nx + ni, 0]
                else:
                    self.boundary_tags_host[element_id, 2] = 3 # Top Tag

    def _load_from_file(self, filename):
        """
        Loads an unstructured quadrilateral mesh from a file using meshio.

        Args:
            filename (str): Path to the mesh file (e.g., .msh).
        
        Raises:
            ValueError: If no quadrilateral cells are found.
        """
        mesh = meshio.read(filename)
        
        # Store Physical Groups Mapping: Name -> Tag
        self.physical_groups = {}
        if hasattr(mesh, 'field_data'):
            for name, data in mesh.field_data.items():
                # data is [tag, dim]
                self.physical_groups[name] = data[0]
        
        print(f"  Physical Groups: {self.physical_groups}")
        
        # 1. Extract Quads
        quads = None
        for cell_block in mesh.cells:
            if cell_block.type == "quad":
                quads = cell_block.data
                break
        
        if quads is None:
            raise ValueError("No 'quad' cells found in the mesh file. This solver supports only Quadrilaterals.")
        
        self.num_elements = len(quads)
        print(f"  Found {self.num_elements} quadrilateral elements.")
        
        # 2. Vertices
        # mesh.points usually 3D (x, y, z). We take x, y.
        points_2d = mesh.points[:, :2]
        
        self.vertices_host = np.zeros((self.num_elements, 4, 2), dtype=np.float32)
        # Copy vertex coordinates for each element
        # quads is (NumElements, 4) indices
        for i, cell in enumerate(quads):
            self.vertices_host[i] = points_2d[cell]

        # 3. Identify Boundary Edges from 'line' cells
        # Map: tuple(sorted_nodes) -> physical_tag
        boundary_edges = {}
        
        # Extract lines
        lines = []
        line_tags = []
        
        # Check where tags are stored (gmsh:physical is common)
        # meshio structures varies: 
        # mesh.cell_data["gmsh:physical"] is a list of arrays corresponding to mesh.cells blocks
        
        for i, cell_block in enumerate(mesh.cells):
            if cell_block.type == "line":
                data = cell_block.data
                lines.append(data)
                
                # Try to get tags
                if "gmsh:physical" in mesh.cell_data:
                    line_tags.append(mesh.cell_data["gmsh:physical"][i])
                else:
                    # Fallback or different format
                    pass
        
        if lines:
            lines_concat = np.concatenate(lines) if len(lines) > 1 else lines[0]
            tags_concat = np.concatenate(line_tags) if len(line_tags) > 1 else line_tags[0]
            
            for edge, tag in zip(lines_concat, tags_concat):
                key = tuple(sorted(edge))
                boundary_edges[key] = tag
                
        print(f"  Found {len(boundary_edges)} boundary edges with tags.")

        # 4. Connectivity (Dual Graph)
        self.connectivity_host = -np.ones((self.num_elements, 4, 2), dtype=np.int32)
        self.boundary_tags_host = np.zeros((self.num_elements, 4), dtype=np.int32)
        
        # Helper to identify faces: defined by sorted tuple of 2 node indices
        # Face definitions for a quad (0,1,2,3):
        # 0: 0-1 (Bottom)
        # 1: 1-2 (Right)
        # 2: 2-3 (Top)
        # 3: 3-0 (Left)
        face_defs = [(0, 1), (1, 2), (2, 3), (3, 0)]
        
        # Map: tuple(sorted_nodes) -> list of (element_id, face_idx)
        face_map = defaultdict(list)
        
        for elem_idx, cell in enumerate(quads):
            for local_face_idx, (n1_idx, n2_idx) in enumerate(face_defs):
                # Global node indices
                g1 = cell[n1_idx]
                g2 = cell[n2_idx]
                face_key = tuple(sorted((g1, g2)))
                face_map[face_key].append((elem_idx, local_face_idx))
        
        # Build connectivity
        for face_key, connections in face_map.items():
            if len(connections) == 2:
                # Internal face
                (e1, f1), (e2, f2) = connections
                self.connectivity_host[e1, f1] = [e2, f2]
                self.connectivity_host[e2, f2] = [e1, f1]
            elif len(connections) == 1:
                # Boundary face
                (e1, f1) = connections[0]
                tag = boundary_edges.get(face_key, 0) # 0 if not found (should not happen if mesh is closed)
                self.boundary_tags_host[e1, f1] = tag
            else:
                print(f"Warning: Face {face_key} shared by {len(connections)} elements. Mesh might be non-manifold.")
        
        # Calculate approximate grid size for dt estimation
        # We take the sqrt of the average element area
        # (This is a rough heuristic, max_wave_speed kernel handles local speed, but we need a length scale)
        total_area = 0.0
        for i in range(self.num_elements):
            v = self.vertices_host[i]
            # Area using Shoelace formula / cross product for convex quads
            # d1 = v2 - v0, d2 = v3 - v1. Area = 0.5 * |d1 x d2|
            d1 = v[2] - v[0]
            d2 = v[3] - v[1]
            area = 0.5 * np.abs(d1[0]*d2[1] - d1[1]*d2[0])
            total_area += area
            
        avg_area = total_area / self.num_elements
        self.dx = np.sqrt(avg_area) # Characteristic length
        print(f"  Approximate element size (dx): {self.dx:.4e}")

    def apply_periodic_condition(self, tag1, tag2, axis, tol=1e-5):
        """
        Applies periodic boundary conditions by linking faces with tag1 to faces with tag2.
        
        Args:
            tag1 (int): Physical tag of the first boundary (e.g., Left).
            tag2 (int): Physical tag of the second boundary (e.g., Right).
            axis (str): The axis of periodicity ('x' or 'y'). Used to ignore the coordinate
                        along the period for matching.
        """
        print(f"Applying Periodic BC: Tag {tag1} <-> Tag {tag2} (Axis {axis})")
        
        # 1. Collect faces for each tag
        # List of (element_idx, face_idx, centroid_coord_to_match)
        faces1 = []
        faces2 = []
        
        # Helper to get face centroid (approximate for matching)
        def get_face_centroid(e, f):
            # Face vertices
            # Connectivity doesn't give vertices directly, need to infer from element
            # Standard order: 0:0-1, 1:1-2, 2:2-3, 3:3-0
            v = self.vertices_host[e]
            if f == 0:   p1, p2 = v[0], v[1]
            elif f == 1: p1, p2 = v[1], v[2]
            elif f == 2: p1, p2 = v[2], v[3]
            elif f == 3: p1, p2 = v[3], v[0]
            
            center = (p1 + p2) * 0.5
            return center

        for e in range(self.num_elements):
            for f in range(4):
                tag = self.boundary_tags_host[e, f]
                if tag == tag1 or tag == tag2:
                    center = get_face_centroid(e, f)
                    
                    # If periodic in X, we match based on Y (and vice versa)
                    match_coord = center[1] if axis == 'x' else center[0]
                    
                    item = (e, f, match_coord)
                    if tag == tag1: faces1.append(item)
                    else: faces2.append(item)
        
        print(f"  Found {len(faces1)} faces for Tag {tag1} and {len(faces2)} for Tag {tag2}.")
        
        if len(faces1) != len(faces2):
            print("  Warning: Number of faces do not match! Periodicity might be incomplete.")
        
        # 2. Match and Link
        # Sort by match coordinate to make matching O(N log N) or simple O(N) linear scan
        faces1.sort(key=lambda x: x[2])
        faces2.sort(key=lambda x: x[2])
        
        linked_count = 0
        
        # Simple greedy matching (assuming sorted)
        # For robust matching, we might need a KD-tree or similar, but 1D sort is fine here.
        idx2 = 0
        for e1, f1, coord1 in faces1:
            # Find closest in faces2
            best_idx = -1
            min_dist = float('inf')
            
            # Search around current index (optimization)
            # Since they are sorted, we can advance idx2
            while idx2 < len(faces2) and faces2[idx2][2] < coord1 - tol:
                idx2 += 1
            
            # Check a window
            start_search = max(0, idx2 - 5)
            end_search = min(len(faces2), idx2 + 10)
            
            for i in range(start_search, end_search):
                dist = abs(faces2[i][2] - coord1)
                if dist < min_dist:
                    min_dist = dist
                    best_idx = i
            
            if min_dist < tol:
                e2, f2, _ = faces2[best_idx]
                
                # Link them!
                # Update Connectivity (Dual Graph)
                # Store (NeighborID, NeighborFaceID)
                # Note: connectivity_host stores [neighbor_elem, neighbor_face_index]
                self.connectivity_host[e1, f1] = [e2, f2]
                self.connectivity_host[e2, f2] = [e1, f1]
                
                # Clear Boundary Tags (Mark as Internal)
                self.boundary_tags_host[e1, f1] = 0
                self.boundary_tags_host[e2, f2] = 0
                
                linked_count += 1
            else:
                print(f"  Warning: No match found for Face {f1} of Element {e1} at coord {coord1}")

        print(f"  Linked {linked_count} periodic pairs.")
        
        # 3. Update Device Arrays
        self.connectivity = wp.array(self.connectivity_host[:, :, 0], dtype=wp.int32, device=self.device)
        self.connectivity_face_indices = wp.array(self.connectivity_host[:, :, 1], dtype=wp.int32, device=self.device)
        self.boundary_tags = wp.array(self.boundary_tags_host, dtype=wp.int32, device=self.device)
