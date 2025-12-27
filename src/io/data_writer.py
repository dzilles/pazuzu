import h5py
import numpy as np
import os
from src.physics.laws.common import conservative_to_primitive

class HDF5Writer:
    """
    Handles writing simulation data to HDF5 format with an accompanying XDMF file for visualization.

    This writer saves the mesh topology/geometry once and appends time-step data as it becomes available.
    Supports high-order sub-cell tessellation for detailed visualization of DG solutions.
    """
    def __init__(self, filename, mesh, basis=None, node_coords=None):
        """
        Initializes the HDF5 Writer and writes the mesh geometry.

        Args:
            filename (str): The path where the .h5 file will be created.
            mesh (Mesh): The simulation mesh.
            basis (Basis, optional): The DG basis. If provided, enables sub-cell visualization.
            node_coords (tuple, optional): (x, y) arrays of shape (NumElements, Np) with physical coordinates.
                                           Required if basis is provided.
        """
        self.filename = filename
        self.mesh = mesh
        self.basis = basis
        self.xmf_filename = filename.replace(".h5", ".xmf")
        self.steps = []
        
        self.use_high_order = (basis is not None and node_coords is not None)
        
        # Ensure output directory exists
        out_dir = os.path.dirname(filename)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        # Explicitly remove old file to ensure no stale datasets remain
        if os.path.exists(self.filename):
            try: os.remove(self.filename)
            except: pass

        # Initialize HDF5 file and write static mesh data
        with h5py.File(self.filename, 'w') as f:
            # Create Mesh Group
            mesh_grp = f.create_group("mesh")
            
            if self.use_high_order:
                # --- High-Order Visualization (Sub-cell tessellation) ---
                N = self.basis.N
                N1 = N + 1
                Np = self.basis.Np
                n_elems = self.mesh.num_elements
                
                # 1. Vertices: All GLL nodes
                # Input coords: (NumElements, Np)
                x_all, y_all = node_coords
                verts_x = x_all.flatten()
                verts_y = y_all.flatten()
                n_points = verts_x.shape[0]
                
                points = np.zeros((n_points, 3), dtype=np.float32)
                points[:, 0] = verts_x
                points[:, 1] = verts_y
                mesh_grp.create_dataset("points", data=points)
                
                # 2. Connectivity: N*N sub-quads per element
                # Each sub-quad (i,j) connects (i,j), (i+1,j), (i+1,j+1), (i,j+1)
                # Local GLL indices
                sub_quads = []
                for j in range(N):
                    for i in range(N):
                        n0 = j * N1 + i
                        n1 = j * N1 + (i + 1)
                        n2 = (j + 1) * N1 + (i + 1)
                        n3 = (j + 1) * N1 + i
                        sub_quads.append([n0, n1, n2, n3])
                sub_quads = np.array(sub_quads, dtype=np.int32) # Shape (N*N, 4)
                
                n_sub_cells = N * N
                
                # Replicate for all elements
                # Global offsets: elem_idx * Np
                offsets = np.arange(n_elems) * Np
                # Broadcast addition: (NumElements, 1, 1) + (1, n_sub_cells, 4) -> (NumElems, n_sub_cells, 4)
                all_conn = offsets[:, None, None] + sub_quads[None, :, :]
                
                conn = all_conn.reshape(-1, 4)
                mesh_grp.create_dataset("connectivity", data=conn)
                
                self.n_vis_elements = conn.shape[0]
                self.n_vis_points = n_points
                
            else:
                # --- Low-Order Visualization (Cell Average) ---
                # Vertices are stored as (NumElements, 4, 2) in the mesh object (Discontinuous).
                verts_2d = mesh.vertices_host.reshape(-1, 2)
                n_points = verts_2d.shape[0]
                points = np.zeros((n_points, 3), dtype=np.float32)
                points[:, 0] = verts_2d[:, 0]
                points[:, 1] = verts_2d[:, 1]
                
                mesh_grp.create_dataset("points", data=points)
                
                # Connectivity: linear
                conn = np.arange(n_points, dtype=np.int32).reshape(-1, 4)
                mesh_grp.create_dataset("connectivity", data=conn)
                
                self.n_vis_elements = conn.shape[0]
                self.n_vis_points = n_points
            
            # Create Data Group for time steps
            f.create_group("data")

    def write_step(self, step, time, Q):
        """
        Writes a single simulation time step to the HDF5 file.
        """
        with h5py.File(self.filename, 'a') as f:
            grp = f["data"].create_group(f"step_{step}")
            grp.attrs["time"] = time
            grp.attrs["step"] = step
            
            if self.use_high_order:
                # --- Write Nodal Data (Point Data) ---
                # Q: (NumElements, Np, 4) -> (NumElements * Np, 4)
                # Convert to Primitive
                q_flat = Q.reshape(-1, 4).T # (4, TotalPoints)
                prim = conservative_to_primitive(q_flat)
                
                grp.create_dataset("rho", data=prim[0])
                grp.create_dataset("u", data=prim[1])
                grp.create_dataset("v", data=prim[2])
                grp.create_dataset("p", data=prim[3])
                
            else:
                # --- Write Cell Average Data (Cell Data) ---
                q_reshaped = Q.transpose(2, 0, 1).reshape(4, -1)
                prim_reshaped = conservative_to_primitive(q_reshaped)
                prim = prim_reshaped.reshape(4, Q.shape[0], Q.shape[1])
                avgs = np.mean(prim, axis=2)
                
                grp.create_dataset("rho", data=avgs[0])
                grp.create_dataset("u", data=avgs[1])
                grp.create_dataset("v", data=avgs[2])
                grp.create_dataset("p", data=avgs[3])
            
        self.steps.append((step, time))
        self._write_xmf()

    def _write_xmf(self):
        """Generates/Updates the XDMF file."""
        with open(self.xmf_filename, 'w') as f:
            f.write('<?xml version="1.0" ?>\n')
            f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
            f.write('<Xdmf Version="3.0">\n')
            f.write(' <Domain>\n')
            f.write('  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n')
            
            h5_rel = os.path.basename(self.filename)
            
            # Determine Attribute Type and Center
            if self.use_high_order:
                attr_type = "Node" # Point Data
                data_dim = self.n_vis_points
            else:
                attr_type = "Cell" # Cell Data
                data_dim = self.n_vis_elements
            
            for step, time in self.steps:
                f.write(f'   <Grid Name="Step_{step}" GridType="Uniform">\n')
                f.write(f'    <Time Value="{time}"/>\n')
                
                # --- Topology ---
                f.write(f'    <Topology TopologyType="Quadrilateral" NumberOfElements="{self.n_vis_elements}">\n')
                f.write(f'     <DataItem Dimensions="{self.n_vis_elements} 4" NumberType="Int" Format="HDF">\n')
                f.write(f'      {h5_rel}:/mesh/connectivity\n')
                f.write(f'     </DataItem>\n')
                f.write(f'    </Topology>\n')
                
                # --- Geometry ---
                f.write(f'    <Geometry GeometryType="XYZ">\n')
                f.write(f'     <DataItem Dimensions="{self.n_vis_points} 3" NumberType="Float" Precision="4" Format="HDF">\n')
                f.write(f'      {h5_rel}:/mesh/points\n')
                f.write(f'     </DataItem>\n')
                f.write(f'    </Geometry>\n')
                
                # --- Attributes ---
                for var in ["rho", "u", "v", "p"]:
                    f.write(f'    <Attribute Name="{var}" AttributeType="Scalar" Center="{attr_type}">\n')
                    f.write(f'     <DataItem Dimensions="{data_dim}" NumberType="Float" Precision="4" Format="HDF">\n')
                    f.write(f'      {h5_rel}:/data/step_{step}/{var}\n')
                    f.write(f'     </DataItem>\n')
                    f.write(f'    </Attribute>\n')
                
                f.write(f'   </Grid>\n')
            
            f.write('  </Grid>\n')
            f.write(' </Domain>\n')
            f.write('</Xdmf>\n')
