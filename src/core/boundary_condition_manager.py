import numpy as np
from src.kernels import boundary_conditions as bc

class BoundaryConditionManager:
    """
    Manages the setup and mapping of boundary conditions.
    Centralizes the logic for mapping configuration strings to integer IDs used by the solver kernels.
    """
    def __init__(self, mesh, config):
        self.mesh = mesh
        self.config = config
        
    def setup_boundary_conditions(self):
        """
        Parses the configuration to setup the boundary condition mask.
        Maps physical tags from the mesh to solver-specific BC IDs.
        
        Returns:
            np.array: Boundary condition mask (host side) with shape (num_elements, 4).
        """
        bc_mask_host = np.zeros((self.mesh.num_elements, 4), dtype=np.int32)
        
        if self.config.boundaries:
            tag_to_bc = {}
            for name, bc_conf in self.config.boundaries.items():
                tag = -1
                # Try to map name to tag using mesh info
                if hasattr(self.mesh, 'physical_groups') and name in self.mesh.physical_groups:
                    tag = self.mesh.physical_groups[name]
                else:
                    try: 
                        tag = int(name) 
                    except ValueError: 
                        pass
                
                if tag != -1:
                    bc_type_str = bc_conf.get('type')
                    bc_id = self._get_bc_id_from_type(bc_type_str)
                    tag_to_bc[tag] = bc_id
            
            # Apply to mask
            for e in range(self.mesh.num_elements):
                for f in range(4):
                    tag = self.mesh.boundary_tags_host[e, f]
                    if tag > 0:
                        # Default to Farfield (Freestream) if tag exists but mapping not found
                        bc_mask_host[e, f] = tag_to_bc.get(tag, bc.BC_FARFIELD)
                        
        return bc_mask_host

    def _get_bc_id_from_type(self, bc_type):
        """
        Maps a string boundary type (from config) to an integer ID (from kernels).
        """
        if bc_type == "slip_wall": 
            return bc.BC_WALL
        elif bc_type == "cylinder_wall": 
            return bc.BC_CYLINDER_WALL
        elif bc_type == "farfield": 
            return bc.BC_FARFIELD
        elif bc_type in ["outflow", "extrapolation"]: 
            return bc.BC_EXTRAPOLATION
        elif bc_type == "inlet": 
            return bc.BC_INLET
        elif bc_type == "outlet": 
            return bc.BC_OUTLET
        elif bc_type == "dmr_exact":
            return bc.BC_DOUBLE_MACH_EXACT
        
        # Default fallback
        return bc.BC_FARFIELD
