import h5py
import numpy as np
import os
from src.physics.equations import conservative_to_primitive

class HDF5Writer:
    def __init__(self, filename, mesh):
        self.filename = filename
        self.mesh = mesh
        self.xmf_filename = filename.replace(".h5", ".xmf")
        self.steps = []
        
        # Ensure output directory exists
        out_dir = os.path.dirname(filename)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        # Initialize HDF5 file and write mesh
        with h5py.File(self.filename, 'w') as f:
            # Create Mesh Group
            mesh_grp = f.create_group("mesh")
            
            # Points: (NumElements * 4, 3)
            # Flatten vertices (NumElems, 4, 2) -> (NumElems*4, 2) -> add z -> (NumElems*4, 3)
            verts_2d = mesh.vertices_host.reshape(-1, 2)
            n_points = verts_2d.shape[0]
            points = np.zeros((n_points, 3), dtype=np.float32)
            points[:, 0] = verts_2d[:, 0]
            points[:, 1] = verts_2d[:, 1]
            
            mesh_grp.create_dataset("points", data=points)
            
            # Connectivity: (NumElements, 4)
            # Just 0,1,2,3 ... because vertices are explicit per element
            conn = np.arange(n_points, dtype=np.int32).reshape(-1, 4)
            mesh_grp.create_dataset("connectivity", data=conn)
            
            # Create Data Group
            f.create_group("data")

    def write_step(self, step, time, Q):
        # Calculate cell averages
        # Q is (NumElements, Np, 4)
        # We want (NumElements, 1) per variable for CellData
        
        # 1. Conservative to Primitive
        # Shape (NumElements, Np, 4) -> (4, NumElements, Np) -> (4, TotalNodes)
        q_reshaped = Q.transpose(2, 0, 1).reshape(4, -1)
        prim = conservative_to_primitive(q_reshaped)
        # -> (4, NumElements, Np)
        prim = prim.reshape(4, Q.shape[0], Q.shape[1])
        
        # 2. Average over Np
        # Shape (4, NumElements)
        avgs = np.mean(prim, axis=2)
        
        rho = avgs[0]
        u   = avgs[1]
        v   = avgs[2]
        p   = avgs[3]
        
        with h5py.File(self.filename, 'a') as f:
            grp = f["data"].create_group(f"step_{step}")
            grp.attrs["time"] = time
            grp.attrs["step"] = step
            
            grp.create_dataset("rho", data=rho)
            grp.create_dataset("u", data=u)
            grp.create_dataset("v", data=v)
            grp.create_dataset("p", data=p)
            
        self.steps.append((step, time))
        self._write_xmf()

    def _write_xmf(self):
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
                
                # Topology
                f.write(f'    <Topology TopologyType="Quadrilateral" NumberOfElements="{n_elems}">\n')
                f.write(f'     <DataItem Dimensions="{n_elems} 4" NumberType="Int" Format="HDF">\n')
                f.write(f'      {h5_rel}:/mesh/connectivity\n')
                f.write(f'     </DataItem>\n')
                f.write(f'    </Topology>\n')
                
                # Geometry
                f.write(f'    <Geometry GeometryType="XYZ">\n')
                f.write(f'     <DataItem Dimensions="{n_points} 3" NumberType="Float" Precision="4" Format="HDF">\n')
                f.write(f'      {h5_rel}:/mesh/points\n')
                f.write(f'     </DataItem>\n')
                f.write(f'    </Geometry>\n')
                
                # Attributes (Cell Data)
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
