import h5py  # type: ignore # Missing type stubs
import numpy as np
import os
import warp as wp
from src.physics.laws.common import conservative_to_primitive

class HDF5Writer:
    """
    Handles writing simulation data to HDF5 format with an accompanying XDMF file for visualization.
    
    Updated for AMR: Each step now writes its own mesh geometry to handle adaptive topology changes.
    """
    def __init__(self, filename, solver):
        self.filename = filename
        self.solver = solver
        self.xmf_filename = filename.replace(".h5", ".xmf")
        self.steps = []
        
        # Ensure output directory exists
        out_dir = os.path.dirname(filename)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        # Precision strings
        self.precision_str = "8" if solver.state.scalar_dtype == wp.float64 else "4"
        self.numpy_dtype = np.float64 if solver.state.scalar_dtype == wp.float64 else np.float32

        # Explicitly remove old file to ensure no stale datasets remain
        if os.path.exists(self.filename):
            try:
                os.remove(self.filename)
            except Exception:
                pass

        # Initialize HDF5 file structure
        with h5py.File(self.filename, 'w') as f:
            f.create_group("data")
            
        # --- Macro Grid (Grid Only) Support ---
        self.grid_filename = self.filename.replace(".h5", "_grid.h5")
        self.grid_xmf_filename = self.grid_filename.replace(".h5", ".xmf")
        self.grid_steps = []
        
        if os.path.exists(self.grid_filename):
            try:
                os.remove(self.grid_filename)
            except Exception:
                pass
                
        with h5py.File(self.grid_filename, 'w') as f:
            f.create_group("data")
        # --------------------------------------

    def _get_mesh_data(self):
        """Generates point and connectivity data for the CURRENT active mesh."""
        solver = self.solver
        num_blocks = solver.quadtree.num_blocks
        active_indices = solver.state.active_block_indices.numpy()[:num_blocks]
        
        N = solver.basis.N
        N1 = N + 1
        Np = solver.basis.Np
        
        # 1. Geometry (Points)
        # Fetch only active blocks from pool
        x_pool = solver.state.x.numpy()
        y_pool = solver.state.y.numpy()
        
        x_active = x_pool[active_indices, :]
        y_active = y_pool[active_indices, :]
        
        verts_x = x_active.flatten()
        verts_y = y_active.flatten()
        n_points = verts_x.shape[0]
        
        points = np.zeros((n_points, 3), dtype=self.numpy_dtype)
        points[:, 0] = verts_x
        points[:, 1] = verts_y
        points[:, 2] = 0.0
        
        # 2. Topology (Connectivity)
        sub_quads = []
        for j in range(N):
            for i in range(N):
                n0 = j * N1 + i
                n1 = j * N1 + (i + 1)
                n2 = (j + 1) * N1 + (i + 1)
                n3 = (j + 1) * N1 + i
                sub_quads.append([n0, n1, n2, n3])
        sub_quads = np.array(sub_quads, dtype=np.int32)
        
        offsets = np.arange(num_blocks) * Np
        all_conn = offsets[:, None, None] + sub_quads[None, :, :]
        connectivity = all_conn.reshape(-1, 4)
        n_elements = connectivity.shape[0]
        
        return points, connectivity, n_points, n_elements

    def write_macro_structure(self, step, time):
        """
        Writes a lightweight 'Grid Only' file containing one element per block.
        Useful for visualizing AMR structure without the cost of high-order sub-cells.
        """
        # Reuse point generation logic (all nodes) to avoid code duplication
        # We accept that we write Np nodes but only reference 4 per block.
        points, _, n_points, _ = self._get_mesh_data()
        
        solver = self.solver
        num_blocks = solver.quadtree.num_blocks
        active_indices = solver.state.active_block_indices.numpy()[:num_blocks]
        
        # Construct Macro Connectivity (1 Quad per Block)
        # Corners relative to start of block node list
        N = solver.basis.N
        N1 = N + 1
        Np = solver.basis.Np
        
        # CCW winding: (0,0), (N,0), (N,N), (0,N)
        # Flat indices: 0, N, N*(N+1)+N, N*(N+1)
        c0 = 0
        c1 = N
        c2 = N * N1 + N
        c3 = N * N1
        corners = np.array([c0, c1, c2, c3], dtype=np.int32)
        
        offsets = np.arange(num_blocks, dtype=np.int32) * Np
        macro_conn = offsets[:, None] + corners[None, :]
        macro_conn = macro_conn.reshape(-1, 4)
        n_elements = num_blocks # One element per block
        
        # Fetch Metadata
        levels = solver.quadtree.block_levels.numpy()
        block_levels = levels[active_indices] # (num_blocks,)
        
        # Write to Grid HDF5
        with h5py.File(self.grid_filename, 'a') as f:
            step_grp = f["data"].create_group(f"step_{step}")
            step_grp.attrs["time"] = time
            step_grp.attrs["step"] = step
            
            # Write Mesh
            mesh_grp = step_grp.create_group("mesh")
            mesh_grp.create_dataset("points", data=points)
            mesh_grp.create_dataset("connectivity", data=macro_conn)
            
            # Write Meta Attributes
            step_grp.create_dataset("level", data=block_levels)
            step_grp.create_dataset("element_id", data=active_indices)
            
        self.grid_steps.append({
            "step": step,
            "time": time,
            "n_points": n_points,
            "n_elements": n_elements
        })
        self._write_macro_xmf()

    def _write_macro_xmf(self):
        """Generates the XDMF for the macro grid."""
        with open(self.grid_xmf_filename, 'w') as f:
            f.write('<?xml version="1.0" ?>\n')
            f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
            f.write('<Xdmf Version="3.0">\n')
            f.write(' <Domain>\n')
            f.write('  <Grid Name="GridSeries" GridType="Collection" CollectionType="Temporal">\n')
            
            h5_rel = os.path.basename(self.grid_filename)
            
            for info in self.grid_steps:
                step = info["step"]
                time = info["time"]
                n_p = info["n_points"]
                n_e = info["n_elements"]
                
                f.write(f'   <Grid Name="Step_{step}" GridType="Uniform">\n')
                f.write(f'    <Time Value="{time}"/>\n')
                
                # --- Topology (Macro Quads) ---
                f.write(f'    <Topology TopologyType="Quadrilateral" NumberOfElements="{n_e}">\n')
                f.write(f'     <DataItem Dimensions="{n_e} 4" NumberType="Int" Format="HDF">\n')
                f.write(f'      {h5_rel}:/data/step_{step}/mesh/connectivity\n')
                f.write('     </DataItem>\n')
                f.write('    </Topology>\n')
                
                # --- Geometry ---
                f.write('    <Geometry GeometryType="XYZ">\n')
                f.write(f'     <DataItem Dimensions="{n_p} 3" NumberType="Float" Precision="{self.precision_str}" Format="HDF">\n')
                f.write(f'      {h5_rel}:/data/step_{step}/mesh/points\n')
                f.write('     </DataItem>\n')
                f.write('    </Geometry>\n')
                
                # --- Attributes (Cell-Centered) ---
                for var in ["level", "element_id"]:
                    f.write(f'    <Attribute Name="{var}" AttributeType="Scalar" Center="Cell">\n')
                    f.write(f'     <DataItem Dimensions="{n_e}" NumberType="Int" Format="HDF">\n')
                    f.write(f'      {h5_rel}:/data/step_{step}/{var}\n')
                    f.write('     </DataItem>\n')
                    f.write('    </Attribute>\n')
                
                f.write('   </Grid>\n')
            
            f.write('  </Grid>\n')
            f.write(' </Domain>\n')
            f.write('</Xdmf>\n')

    def write_step(self, step, time):
        """Writes current mesh and state data."""
        # Generate mesh for THIS step
        points, connectivity, n_points, n_elements = self._get_mesh_data()
        
        # Fetch current state
        solver = self.solver
        num_blocks = solver.quadtree.num_blocks
        active_indices = solver.state.active_block_indices.numpy()[:num_blocks]
        
        q_pool = solver.state.q.numpy()
        q_active = q_pool[active_indices, :, :]
        q_flat = q_active.reshape(-1, 4)
        
        # Convert to Primitive
        prim = conservative_to_primitive(q_flat.T)
        
        # Fetch phi (SDF)
        phi_pool = solver.state.phi.numpy()
        phi_active = phi_pool[active_indices, :]
        phi_flat = phi_active.flatten()

        # --- [NEW] Generate Boundary ID Field (bc_id) ---
        # Initialize all nodes to -1 (Internal)
        Np = solver.basis.Np
        bc_id_active = np.full((num_blocks, Np), -1.0, dtype=np.float32)
        
        # Fetch bc_mask (Shape: max_blocks x 4)
        if hasattr(solver.state, 'bc_mask'):
            bc_mask_pool = solver.state.bc_mask.numpy()
            bc_mask_active = bc_mask_pool[active_indices, :] # (num_blocks, 4)
            
            # Helper: face_nodes maps face_idx (0..3) to list of node indices
            # Shape: (4, Nfp)
            face_nodes = solver.basis.face_nodes.numpy()
            
            # Vectorized tagging
            for f in range(4):
                # Identify blocks that have a BC on face f
                # Mask value > -1 implies a BC ID
                mask_f = bc_mask_active[:, f]
                has_bc = mask_f > -1
                
                if np.any(has_bc):
                    block_idxs = np.where(has_bc)[0]
                    bc_values = mask_f[block_idxs]
                    nodes_on_face = face_nodes[f]
                    
                    # Broadcast assignment:
                    # For every block in block_idxs, set the nodes 'nodes_on_face' to the corresponding bc_value
                    # Use advanced indexing: [rows, columns]
                    # rows: block_idxs (broadcasted against nodes)
                    # cols: nodes_on_face
                    bc_id_active[block_idxs[:, None], nodes_on_face] = bc_values[:, None]

        bc_id_flat = bc_id_active.flatten()
        # -----------------------------------------------

        with h5py.File(self.filename, 'a') as f:
            step_grp = f["data"].create_group(f"step_{step}")
            step_grp.attrs["time"] = time
            step_grp.attrs["step"] = step
            
            # Write Mesh
            mesh_grp = step_grp.create_group("mesh")
            mesh_grp.create_dataset("points", data=points)
            mesh_grp.create_dataset("connectivity", data=connectivity)
            
            # Write Fields
            step_grp.create_dataset("rho", data=prim[0])
            step_grp.create_dataset("u", data=prim[1])
            step_grp.create_dataset("v", data=prim[2])
            step_grp.create_dataset("p", data=prim[3])
            step_grp.create_dataset("phi", data=phi_flat)
            # [NEW] Write BC ID
            step_grp.create_dataset("bc_id", data=bc_id_flat)
            
        self.steps.append({
            "step": step,
            "time": time,
            "n_points": n_points,
            "n_elements": n_elements
        })
        self._write_xmf()

    def _write_xmf(self):
        """Generates/Updates the temporal XDMF file."""
        with open(self.xmf_filename, 'w') as f:
            f.write('<?xml version="1.0" ?>\n')
            f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
            f.write('<Xdmf Version="3.0">\n')
            f.write(' <Domain>\n')
            f.write('  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n')
            
            h5_rel = os.path.basename(self.filename)
            
            for info in self.steps:
                step = info["step"]
                time = info["time"]
                n_p = info["n_points"]
                n_e = info["n_elements"]
                
                f.write(f'   <Grid Name="Step_{step}" GridType="Uniform">\n')
                f.write(f'    <Time Value="{time}"/>\n')
                
                # --- Topology ---
                f.write(f'    <Topology TopologyType="Quadrilateral" NumberOfElements="{n_e}">\n')
                f.write(f'     <DataItem Dimensions="{n_e} 4" NumberType="Int" Format="HDF">\n')
                f.write(f'      {h5_rel}:/data/step_{step}/mesh/connectivity\n')
                f.write('     </DataItem>\n')
                f.write('    </Topology>\n')
                
                # --- Geometry ---
                f.write('    <Geometry GeometryType="XYZ">\n')
                f.write(f'     <DataItem Dimensions="{n_p} 3" NumberType="Float" Precision="{self.precision_str}" Format="HDF">\n')
                f.write(f'      {h5_rel}:/data/step_{step}/mesh/points\n')
                f.write('     </DataItem>\n')
                f.write('    </Geometry>\n')
                
                # --- Attributes ---
                for var in ["rho", "u", "v", "p", "phi", "bc_id"]:
                    f.write(f'    <Attribute Name="{var}" AttributeType="Scalar" Center="Node">\n')
                    f.write(f'     <DataItem Dimensions="{n_p}" NumberType="Float" Precision="{self.precision_str}" Format="HDF">\n')
                    f.write(f'      {h5_rel}:/data/step_{step}/{var}\n')
                    f.write('     </DataItem>\n')
                    f.write('    </Attribute>\n')
                
                f.write('   </Grid>\n')
            
            f.write('  </Grid>\n')
            f.write(' </Domain>\n')
            f.write('</Xdmf>\n')


class HDF5Reader:
    """
    Handles reading simulation state from HDF5 files for restarting.
    """
    @staticmethod
    def load_checkpoint(filename, step=None):
        """
        Loads the simulation state from an HDF5 checkpoint.
        
        Args:
            filename (str): Path to the .h5 file.
            step (int, optional): Specific step to load. If None, loads the last available step.
            
        Returns:
            dict: {
                "time": float,
                "step": int,
                "q_prim_flat": np.array (shape (4, TotalPoints))
            }
        """
        if not os.path.exists(filename):
            raise FileNotFoundError(f"Checkpoint file not found: {filename}")
            
        with h5py.File(filename, 'r') as f:
            if "data" not in f:
                raise ValueError(f"Invalid checkpoint file structure: 'data' group missing in {filename}")

            # Find step
            if step is None:
                steps = []
                for k in f['data'].keys():
                    if k.startswith('step_'):
                        try:
                            steps.append(int(k.split('_')[1]))
                        except ValueError:
                            pass
                
                if not steps:
                    raise ValueError("No time steps found in checkpoint file.")
                step = max(steps)
            
            grp_name = f"step_{step}"
            grp_path = f"data/{grp_name}"
            
            if grp_path not in f:
                 raise ValueError(f"Step {step} not found in {filename}")
            
            grp = f[grp_path]
            time = grp.attrs.get("time", 0.0)
            loaded_step = grp.attrs.get("step", step)
            
            # Read variables
            # Assuming high-order storage: (TotalPoints,) for each variable
            if "rho" not in grp or "u" not in grp or "v" not in grp or "p" not in grp:
                raise ValueError(f"Missing flow variables (rho, u, v, p) in step {step}")

            rho = grp["rho"][:]
            u = grp["u"][:]
            v = grp["v"][:]
            p = grp["p"][:]
            
            # Stack: (4, TotalPoints)
            # Ensure they are 1D flat arrays
            prim_flat = np.vstack([
                rho.flatten(), 
                u.flatten(), 
                v.flatten(), 
                p.flatten()
            ])
            
            return {
                "time": time,
                "step": loaded_step,
                "q_prim_flat": prim_flat
            }