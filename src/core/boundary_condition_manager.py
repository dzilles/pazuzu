import numpy as np
from src.kernels import boundary_conditions as bc

class BoundaryConditionManager:
    """
    Manages the setup, validation, and mapping of boundary conditions.
    Centralizes the logic for mapping configuration strings to integer IDs and parameters
    used by the solver kernels.
    """
    REQUIRED_PARAMS = {
        bc.BC_INLET: ["rho", "u", "v", "p"],
        bc.BC_OUTLET: ["p_back"],
        bc.BC_ISOTHERMAL_WALL: ["T_wall"]
    }

    def __init__(self, mesh, config):
        self.mesh = mesh
        self.config = config
        
    def setup_boundary_conditions(self):
        """
        Parses and validates the configuration to setup the boundary condition mask and data.
        
        Returns:
            tuple: (bc_mask_host, bc_data_list)
                bc_mask_host (np.array): Index into bc_data_list per face (num_elements, 4).
                bc_data_list (list): List of dictionaries containing BC type and parameters.
        """
        # bc_mask will store the index into bc_data_list
        bc_mask_host = -np.ones((self.mesh.num_elements, 4), dtype=np.int32)
        bc_data_list = []
        
        # Map: tag -> index in bc_data_list
        tag_to_index = {}

        if self.config.boundaries:
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
                    
                    # Validate and extract parameters
                    params = bc_conf.get('params', {})
                    validated_params = self._validate_params(name, bc_id, params)
                    
                    # Create data entry
                    bc_data_list.append({
                        'type': bc_id,
                        'params': validated_params
                    })
                    tag_to_index[tag] = len(bc_data_list) - 1
            
            # Apply to mask
            for e in range(self.mesh.num_elements):
                for f in range(4):
                    tag = self.mesh.boundary_tags_host[e, f]
                    if tag > 0:
                        if tag in tag_to_index:
                            # If it's periodic but tag is still > 0, it means linking failed or wasn't performed
                            data_idx = tag_to_index[tag]
                            if bc_data_list[data_idx]['type'] == bc.BC_PERIODIC:
                                raise ValueError(f"Boundary Tag {tag} is marked as periodic but remains unlinked in the mesh. Ensure 'linked_to' is specified correctly and 'apply_periodic_condition' was called.")
                            
                            bc_mask_host[e, f] = data_idx
                        else:
                            # Fallback if tag is in mesh but not in config: use Farfield with defaults
                            # Create a unique entry for this tag
                            bc_data_list.append({
                                'type': bc.BC_FARFIELD,
                                'params': self._get_farfield_defaults()
                            })
                            tag_to_index[tag] = len(bc_data_list) - 1
                            bc_mask_host[e, f] = tag_to_index[tag]
                        
        return bc_mask_host, bc_data_list

    def _validate_params(self, name, bc_id, params):
        """
        Validates that all required parameters for a BC type are present.
        Raises ValueError if any are missing.
        """
        required = self.REQUIRED_PARAMS.get(bc_id, [])
        for p in required:
            if p not in params:
                raise ValueError(f"Boundary '{name}' of type '{self._get_type_name(bc_id)}' is missing required parameter: '{p}'")
        
        # Specific handling for Farfield: default to freestream if not provided
        if bc_id == bc.BC_FARFIELD:
            full_params = self._get_farfield_defaults()
            full_params.update(params)
            return full_params
            
        return params

    def _get_farfield_defaults(self):
        """Returns default values for farfield boundaries from physics config."""
        return {
            "rho": self.config.physics.rho_inf,
            "u": self.config.physics.u_inf,
            "v": self.config.physics.v_inf,
            "p": self.config.physics.p_inf
        }

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
        elif bc_type == "periodic":
            return bc.BC_PERIODIC
        elif bc_type in ["outflow", "extrapolation"]: 
            return bc.BC_EXTRAPOLATION
        elif bc_type == "inlet": 
            return bc.BC_INLET
        elif bc_type == "outlet": 
            return bc.BC_OUTLET
        elif bc_type == "dmr_exact":
            return bc.BC_DOUBLE_MACH_EXACT
        elif bc_type == "no_slip_wall":
            return bc.BC_NO_SLIP_WALL
        elif bc_type == "isothermal_wall":
            return bc.BC_ISOTHERMAL_WALL
        
        # Default fallback
        return bc.BC_FARFIELD

    def _get_type_name(self, bc_id):
        """Helper to get string name from ID for error messages."""
        for name, val in bc.__dict__.items():
            if name.startswith("BC_") and val == bc_id:
                return name[3:].lower()
        return "unknown"
