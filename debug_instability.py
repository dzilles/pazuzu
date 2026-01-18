import warp as wp
import numpy as np
import sys
import os
import math

# Add source to path
sys.path.append(os.getcwd())

from solver import PazuzuSolver

def analyze_ghost_node(solver, pool_idx, node_idx):
    """
    Performs a forensic analysis of a specific ghost node to see 
    where it is getting its value from (Feedback Loop Check).
    """
    # 1. Fetch Data
    x_arr = solver.state.x.numpy()
    y_arr = solver.state.y.numpy()
    phi_arr = solver.state.phi.numpy()
    q_arr = solver.state.q.numpy()
    nodes_1d = solver.basis.nodes_1d.numpy()
    
    # 2. Identify the Ghost Node
    gx = x_arr[pool_idx, node_idx]
    gy = y_arr[pool_idx, node_idx]
    gphi = phi_arr[pool_idx, node_idx]
    g_rho = q_arr[pool_idx, node_idx][0]
    
    print(f"\n[Forensic Analysis] Ghost Node (Block {pool_idx}, Node {node_idx})")
    print(f"  Position: ({gx:.6f}, {gy:.6f})")
    print(f"  SDF (phi): {gphi:.6f} (Should be < 0)")
    print(f"  Current Density: {g_rho:.6f}")

    # 3. Re-calculate Image Point Geometry (Python implementation of Kernel Logic)
    # Note: This assumes Cylinder params from config for simplicity. 
    # If using STL, this part is an approximation or requires mesh query.
    # We'll assume Analytical Cylinder for the debug trace logic:
    geom = solver.config.ibm.geometric_params
    cx = geom.get("center_x", 0.0)
    cy = geom.get("center_y", 0.0)
    radius = geom.get("radius", 0.5)
    
    dx = gx - cx
    dy = gy - cy
    dist_center = math.sqrt(dx*dx + dy*dy)
    nx = dx / dist_center
    ny = dy / dist_center
    
    # Wall Distance d = r - dist_center (since inside)
    # Image point is at d into the fluid
    # But wait, SDF is dist - radius. Inside, SDF is negative.
    # d_surf = |phi|
    d_surf = abs(gphi)
    
    # Image Point = Surface Point + Normal * d_surf
    # Surface Point = Center + Radius * Normal
    sx = cx + nx * radius
    sy = cy + ny * radius
    
    ix = sx + nx * d_surf
    iy = sy + ny * d_surf
    
    print(f"  Image Point (Calculated): ({ix:.6f}, {iy:.6f})")
    
    # 4. Find which block contains the Image Point
    # Brute-force search in active blocks (slow but accurate for debug)
    found_block = -1
    
    # Get active block indices
    active_indices = solver.state.active_block_indices.numpy()[:solver.quadtree.num_blocks]
    
    for b_idx in active_indices:
        # Check bounding box of the block nodes (approximate)
        bx = x_arr[b_idx]
        by = y_arr[b_idx]
        min_x, max_x = np.min(bx), np.max(bx)
        min_y, max_y = np.min(by), np.max(by)
        
        # Add a tiny tolerance
        tol = 1e-4
        if (ix >= min_x - tol and ix <= max_x + tol and
            iy >= min_y - tol and iy <= max_y + tol):
            found_block = b_idx
            break
            
    if found_block == -1:
        print("  [CRITICAL] Image Point falls OUTSIDE all active blocks!")
        return

    print(f"  Image Point falls into Block {found_block}")
    
    # 5. Check Neighbors in that block
    # Find the 4 closest nodes (Bilinear Stencil)
    # We do this by mapping (ix, iy) to reference space (-1, 1)
    bx = x_arr[found_block]
    by = y_arr[found_block]
    
    # Simple bounds (assuming Cartesian aligned for debug simplicity)
    # For general meshes, we'd need the Jacobian, but let's look at the PHI of the whole block
    b_phi = phi_arr[found_block]
    
    print(f"  [Stencil Check] Neighbors in Block {found_block}:")
    
    solid_neighbors = 0
    fluid_neighbors = 0
    total_nodes = len(b_phi)
    
    min_dist = 1e9
    closest_node_phi = 0.0
    closest_node_idx = -1
    
    for i in range(total_nodes):
        nx_pos = bx[i]
        ny_pos = by[i]
        dist = math.sqrt((nx_pos - ix)**2 + (ny_pos - iy)**2)
        
        if dist < min_dist:
            min_dist = dist
            closest_node_phi = b_phi[i]
            closest_node_idx = i
            
    print(f"    Closest Node Index: {closest_node_idx}")
    print(f"    Closest Node Distance: {min_dist:.6f}")
    print(f"    Closest Node PHI: {closest_node_phi:.6f}")
    
    if closest_node_phi < 1e-9:
        print("    [ALERT] The closest interpolation node is SOLID/GHOST (Phi <= 0).")
        print("            This is a high risk for Feedback Loop.")
    else:
        print("    [OK] The closest interpolation node is FLUID (Phi > 0).")

    # Check for NaNs in the block
    b_q = q_arr[found_block]
    if np.any(np.isnan(b_q)):
        print("    [CRITICAL] The Target Block contains NaNs!")

def run_debug_simulation(config_path):
    print(f"Initializing Debug Session for: {config_path}")
    solver = PazuzuSolver(config_path)
    
    # Run until t > 45.0 or crash
    t_max = 45.0
    print(f"Running simulation step-by-step until t > {t_max}...")
    
    try:
        step = 0
        while solver.state.t < t_max:
            step += 1
            solver.state.step += 1
            t = solver.state.t
            
            # 1. Compute DT
            dt = solver.compute_dt()
            if dt < 1e-9:
                print(f"[STOP] Timestep collapsed (dt={dt:.2e}) at step {step}.")
                break
                
            # 2. Run Step
            if solver.config.numerics.time_integrator == "rk4":
                solver.integrator.step_rk4(
                    solver.compute_rhs,
                    dt,
                    t,
                    solver.quadtree.num_blocks
                )
            else:
                solver.integrator.step_ssp_rk3(
                    solver.compute_rhs,
                    dt,
                    t,
                    solver.quadtree.num_blocks
                )
            solver.state.t += dt
            
            # 3. Monitor Stability
            # Check for NaNs
            solver.nan_flag.zero_()
            from src.kernels.common_kernels import check_nan_indirect
            wp.launch(
                kernel=check_nan_indirect,
                dim=solver.quadtree.num_blocks * solver.basis.Np,
                inputs=[solver.state.q, solver.state.active_block_indices, solver.nan_flag],
                device=solver.device
            )
            
            if solver.nan_flag.numpy()[0] > 0:
                print(f"[CRASH] NaN detected at Step {step} (t={t:.4f})!")
                break
            
            # 4. Check Wave Speeds (Warning Signs)
            if step % 500 == 0:
                max_inv_dt = solver.max_inv_dt.numpy()[0]
                max_wave_speed = max_inv_dt * 0.01 # Approx dx
                print(f"Step {step}: t={t:.4f}, dt={dt:.4e}, MaxWave={max_wave_speed:.2f}")

    except Exception as e:
        print(f"\n[EXCEPTION] Simulation crashed: {e}")

    # --- POST-CRASH FORENSICS ---
    print("\n" + "="*50)
    print("STARTING FORENSIC ANALYSIS")
    print("="*50)
    
    # Pull data
    phi = solver.state.phi.numpy()
    q = solver.state.q.numpy()
    active_indices = solver.state.active_block_indices.numpy()[:solver.quadtree.num_blocks]
    
    # 1. Find the "Hottest" Node (Highest Velocity or NaN)
    max_vel = -1.0
    worst_block = -1
    worst_node = -1
    
    found_nan = False
    
    for b_idx in active_indices:
        b_q = q[b_idx]
        rho = b_q[:, 0]
        u = b_q[:, 1] / rho
        v = b_q[:, 2] / rho
        vel = np.sqrt(u*u + v*v)
        
        # Check NaN
        if np.any(np.isnan(vel)):
            nan_idx = np.where(np.isnan(vel))[0][0]
            worst_block = b_idx
            worst_node = nan_idx
            print(f"Found NaN at Block {b_idx}, Node {nan_idx}")
            found_nan = True
            break
            
        # Check Velocity Spike
        local_max = np.max(vel)
        if local_max > max_vel:
            max_vel = local_max
            idx = np.argmax(vel)
            worst_block = b_idx
            worst_node = idx

    print(f"Worst Node detected at Block {worst_block}, Node {worst_node}")
    if not found_nan:
        print(f"Max Velocity: {max_vel:.2f}")
        
    # 2. Check if this is a Ghost Node
    worst_phi = phi[worst_block, worst_node]
    print(f"Node Type: {'GHOST/SOLID' if worst_phi < 0 else 'FLUID'} (Phi={worst_phi:.4f})")
    
    if worst_phi < 0:
        print("The instability originated in a GHOST NODE.")
        analyze_ghost_node(solver, worst_block, worst_node)
    else:
        print("The instability originated in the FLUID.")
        print("Checking neighbors...")
        # Check if any neighbor is a ghost node (Bounday Layer instability)
        # (Simplified: just analyzing the ghost node logic anyway)
        pass

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python debug_instability.py <config.yaml>")
    else:
        run_debug_simulation(sys.argv[1])
