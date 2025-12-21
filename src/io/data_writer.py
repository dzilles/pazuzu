import h5py
import numpy as np
import os
from src.physics.equations import conservative_to_primitive

class HDF5Writer:
    """
    Handles writing simulation data to HDF5 format with an accompanying XDMF file for visualization.

    This writer saves the mesh topology/geometry once and appends time-step data as it becomes available.
    The output allows for time-series visualization in ParaView.
    Currently, it saves cell-averaged values for visualization on the linear mesh.

    Attributes:
        filename (str): Path to the HDF5 output file (.h5).
        mesh (Mesh): The mesh object containing geometry and connectivity.
        xmf_filename (str): Path to the XDMF metadata file (.xmf).
        steps (list): List of tuples (step_index, time) recorded so far.
    """
    def __init__(self, filename, mesh):
        """
        Initializes the HDF5 Writer and writes the mesh geometry.

        Args:
            filename (str): The path where the .h5 file will be created.
            mesh (Mesh): The simulation mesh.
        """
        self.filename = filename
        self.mesh = mesh
        self.xmf_filename = filename.replace(".h5", ".xmf")
        self.steps = []
        
        # Ensure output directory exists
        out_dir = os.path.dirname(filename)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        # Initialize HDF5 file and write static mesh data
        with h5py.File(self.filename, 'w') as f:
            # Create Mesh Group
            mesh_grp = f.create_group("mesh")
            
            # --- Write Geometry (Points) ---
            # Vertices are stored as (NumElements, 4, 2) in the mesh object (Discontinuous).
            # We flatten this to (NumElements * 4, 2) and add a z-coordinate (0.0) for 3D compatibility.
            verts_2d = mesh.vertices_host.reshape(-1, 2)
            n_points = verts_2d.shape[0]
            points = np.zeros((n_points, 3), dtype=np.float32)
            points[:, 0] = verts_2d[:, 0]
            points[:, 1] = verts_2d[:, 1]
            
            mesh_grp.create_dataset("points", data=points)
            
            # --- Write Topology (Connectivity) ---
            # Since vertices are explicit and duplicated for DG, the connectivity is simply linear:
            # Cell 0 uses points [0, 1, 2, 3], Cell 1 uses [4, 5, 6, 7], etc.
            conn = np.arange(n_points, dtype=np.int32).reshape(-1, 4)
            mesh_grp.create_dataset("connectivity", data=conn)
            
            # Create Data Group for time steps
            f.create_group("data")

    def write_step(self, step, time, Q):
        """
        Writes a single simulation time step to the HDF5 file and updates the XDMF.

        Calculates cell-averaged primitive variables (density, velocity, pressure) from the
        high-order DG conservative state and saves them.

        Args:
            step (int): The current time step index.
            time (float): The current simulation time.
            Q (np.array): The conservative state vector (NumElements, Np, 4).
        """
        # Calculate cell averages for visualization
        # Q shape: (NumElements, Np, 4)
        
        # 1. Convert Conservative (rho, rhou, rhov, E) -> Primitive (rho, u, v, p)
        # Transpose to (4, NumElements, Np) for easy reshaping -> (4, TotalNodes)
        q_reshaped = Q.transpose(2, 0, 1).reshape(4, -1)
        prim_reshaped = conservative_to_primitive(q_reshaped)
        
        # Reshape back to (4, NumElements, Np) -> (NumElements, Np, 4) if needed, 
        # but here we keep (4, NumElements, Np) to average over axis 2 (Np)
        prim = prim_reshaped.reshape(4, Q.shape[0], Q.shape[1])
        
        # 2. Average over all Np nodes in each element to get one value per cell
        # Result shape: (4, NumElements)
        avgs = np.mean(prim, axis=2)
        
        rho = avgs[0]
        u   = avgs[1]
        v   = avgs[2]
        p   = avgs[3]
        
        # Write to HDF5
        with h5py.File(self.filename, 'a') as f:
            grp = f["data"].create_group(f"step_{step}")
            grp.attrs["time"] = time
            grp.attrs["step"] = step
            
            grp.create_dataset("rho", data=rho)
            grp.create_dataset("u", data=u)
            grp.create_dataset("v", data=v)
            grp.create_dataset("p", data=p)
            
        self.steps.append((step, time))
        
        # Regenerate XDMF file to include the new step
        self._write_xmf()

    def _write_xmf(self):
        """
        Generates/Updates the XDMF file that points to the HDF5 data.
        
        The XDMF file defines the structure of the data for ParaView, linking the 
        static mesh topology with the time-dependent attribute data.
        """
        with open(self.xmf_filename, 'w') as f:
            f.write('<?xml version="1.0" ?>\n')
            f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
            f.write('<Xdmf Version="3.0">\n')
            f.write(' <Domain>\n')
            f.write('  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n')
            
            h5_rel = os.path.basename(self.filename)
            n_elems = self.mesh.num_elements
            n_points = n_elems * 4
            
            for step, time in self.steps:
                f.write(f'   <Grid Name="Step_{step}" GridType="Uniform">\n')
                f.write(f'    <Time Value="{time}"/>\n')
                
                # --- Topology ---
                f.write(f'    <Topology TopologyType="Quadrilateral" NumberOfElements="{n_elems}">\n')
                f.write(f'     <DataItem Dimensions="{n_elems} 4" NumberType="Int" Format="HDF">\n')
                f.write(f'      {h5_rel}:/mesh/connectivity\n')
                f.write(f'     </DataItem>\n')
                f.write(f'    </Topology>\n')
                
                # --- Geometry ---
                f.write(f'    <Geometry GeometryType="XYZ">\n')
                f.write(f'     <DataItem Dimensions="{n_points} 3" NumberType="Float" Precision="4" Format="HDF">\n')
                f.write(f'      {h5_rel}:/mesh/points\n')
                f.write(f'     </DataItem>\n')
                f.write(f'    </Geometry>\n')
                
                # --- Attributes (Cell Data) ---
                for var in ["rho", "u", "v", "p"]:
                    f.write(f'    <Attribute Name="{var}" AttributeType="Scalar" Center="Cell">\n')
                    f.write(f'     <DataItem Dimensions="{n_elems}" NumberType="Float" Precision="4" Format="HDF">\n')
                    f.write(f'      {h5_rel}:/data/step_{step}/{var}\n')
                    f.write(f'     </DataItem>\n')
                    f.write(f'    </Attribute>\n')
                
                f.write(f'   </Grid>\n')
            
            f.write('  </Grid>\n')
            f.write(' </Domain>\n')
            f.write('</Xdmf>\n')
