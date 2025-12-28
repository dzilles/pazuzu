import h5py
import numpy as np
import os
import warp as wp
from src.physics.laws.common import conservative_to_primitive

class HDF5Writer:
    """
    Handles writing simulation data to HDF5 format with an accompanying XDMF file for visualization.
    
    Adapted for the Quadtree-based PazuzuSolver.
    Writes high-order GLL nodes as explicit points, allowing for detailed sub-cell visualization in Paraview.
    """
    def __init__(self, filename, solver):
        """
        Initializes the HDF5 Writer and writes the mesh geometry.

        Args:
            filename (str): The path where the .h5 file will be created.
            solver (PazuzuSolver): The solver instance containing state and geometry.
        """
        self.filename = filename
        self.solver = solver
        self.xmf_filename = filename.replace(".h5", ".xmf")
        self.steps = []
        
        # Ensure output directory exists
        out_dir = os.path.dirname(filename)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        # Explicitly remove old file to ensure no stale datasets remain
        if os.path.exists(self.filename):
            try: os.remove(self.filename)
            except: pass

        # Extract Geometry from Solver
        # We assume the mesh is static for now (no adaptive refinement *during* the run that changes the writer)
        # If AMR happens, we'd need to re-write geometry or use XMF's temporal grid support fully (heavy).
        # For Phase 2, we assume the initial mesh (after refinement) stays fixed or we only write the initial state layout.
        # NOTE: If AMR changes the mesh, this writer needs to be re-initialized or updated per step. 
        # For now, we capture the mesh at __init__.
        
        num_blocks = solver.quadtree.num_blocks
        N = solver.basis.N
        N1 = N + 1
        Np = solver.basis.Np
        
        # Fetch coordinates for active blocks
        # state.x: (MAX_BLOCKS, Np) -> slice to (num_blocks, Np)
        x_wp = solver.state.x
        y_wp = solver.state.y
        
        # We need to copy to host
        # Slicing in Warp python API: array[start:end]
        # But Warp arrays are 1D in some contexts or multidimensional. 
        # state.x is (MAX_BLOCKS, Np). We can't easily slice distinct rows if they are not contiguous in memory 
        # (they are in the pool), but we can just copy the first num_blocks * Np elements if we view it as flat,
        # OR copy the whole thing and slice in numpy.
        # Since MAX_BLOCKS isn't huge (10k * 16 * 4 bytes ~ 640KB), copying all is fine.
        
        x_all = x_wp.numpy()[:num_blocks, :] # (num_blocks, Np)
        y_all = y_wp.numpy()[:num_blocks, :]
        
        verts_x = x_all.flatten()
        verts_y = y_all.flatten()
        n_points = verts_x.shape[0]
        
        # Initialize HDF5 file and write static mesh data
        with h5py.File(self.filename, 'w') as f:
            # Create Mesh Group
            mesh_grp = f.create_group("mesh")
            
            # 1. Vertices (Geometry)
            points = np.zeros((n_points, 3), dtype=np.float32)
            points[:, 0] = verts_x
            points[:, 1] = verts_y
            points[:, 2] = 0.0
            mesh_grp.create_dataset("points", data=points)
            
            # 2. Connectivity (Topology)
            # We construct sub-quads for each high-order element.
            # Grid is N x N sub-cells per block.
            # Local node indices in a block:
            # (j, i) -> j * N1 + i
            
            sub_quads = []
            for j in range(N):
                for i in range(N):
                    n0 = j * N1 + i
                    n1 = j * N1 + (i + 1)
                    n2 = (j + 1) * N1 + (i + 1)
                    n3 = (j + 1) * N1 + i
                    sub_quads.append([n0, n1, n2, n3])
            sub_quads = np.array(sub_quads, dtype=np.int32) # Shape (N*N, 4)
            
            n_sub_cells_per_block = N * N
            
            # Replicate for all blocks
            # Global offsets: block_idx * Np
            offsets = np.arange(num_blocks) * Np
            # Broadcast addition: (num_blocks, 1, 1) + (1, n_sub_cells, 4)
            all_conn = offsets[:, None, None] + sub_quads[None, :, :]
            
            conn = all_conn.reshape(-1, 4)
            mesh_grp.create_dataset("connectivity", data=conn)
            
            self.n_vis_elements = conn.shape[0]
            self.n_vis_points = n_points
            
            # Create Data Group for time steps
            f.create_group("data")

    def write_step(self, step, time):
        """
        Writes a single simulation time step to the HDF5 file.
        """
        # Fetch current state
        # q: (MAX_BLOCKS, Np, 4)
        num_blocks = self.solver.quadtree.num_blocks
        q_wp = self.solver.state.q
        
        # Copy to host
        # Shape (MAX_BLOCKS, Np, 4)
        q_all = q_wp.numpy()
        q_active = q_all[:num_blocks, :, :] # (num_blocks, Np, 4)
        
        # Flatten to (TotalPoints, 4)
        q_flat = q_active.reshape(-1, 4) # (n_points, 4)
        
        # Convert to Primitive variables for visualization
        # Helper expects (4, N)
        q_transposed = q_flat.T # (4, n_points)
        prim = conservative_to_primitive(q_transposed) # (4, n_points)
        
        with h5py.File(self.filename, 'a') as f:
            grp = f["data"].create_group(f"step_{step}")
            grp.attrs["time"] = time
            grp.attrs["step"] = step
            
            grp.create_dataset("rho", data=prim[0])
            grp.create_dataset("u", data=prim[1])
            grp.create_dataset("v", data=prim[2])
            grp.create_dataset("p", data=prim[3])
            
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
            
            # Data is Nodal (Point)
            attr_type = "Node" 
            data_dim = self.n_vis_points
            
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