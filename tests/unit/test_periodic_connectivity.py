import pytest
import warp as wp
import numpy as np
from src.core.basis import Basis
from src.core.simulation_state import SimulationState
from src.geometry.quadtree import Quadtree

def test_periodic_connectivity_2x2(device="cpu"):
    wp.init()
    basis = Basis(polynomial_degree=1, device=device)
    max_blocks = 16
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=max_blocks)
    
    # Enable Periodicity
    quadtree = Quadtree(device=device, max_blocks=max_blocks, periodic_x=True, periodic_y=True)
    
    # 2x2 Grid (Level 1)
    # Morton: BL=4, BR=5, TL=6, TR=7
    # (ix, iy): (0,0)=4, (1,0)=5, (0,1)=6, (1,1)=7
    
    quadtree.uniform_refine(1, state, basis)
    
    neighbors = state.neighbors.numpy()
    codes = quadtree.block_morton_codes.numpy()
    
    # Identify pool indices
    idx_4 = -1; idx_5 = -1; idx_6 = -1; idx_7 = -1
    for i in range(4):
        if codes[i] == 4: idx_4 = i
        if codes[i] == 5: idx_5 = i
        if codes[i] == 6: idx_6 = i
        if codes[i] == 7: idx_7 = i
        
    # Block 4 (0,0) Neighbors:
    # Left (-x) -> (1,0) = 5
    # Right (+x) -> (1,0) = 5
    # Bottom (-y) -> (0,1) = 6
    # Top (+y) -> (0,1) = 6
    assert neighbors[idx_4, 0] == idx_5
    assert neighbors[idx_4, 1] == idx_5
    assert neighbors[idx_4, 2] == idx_6
    assert neighbors[idx_4, 3] == idx_6

    # Block 7 (1,1) Neighbors:
    # Left (-x) -> (0,1) = 6
    # Right (+x) -> (0,1) = 6
    # Bottom (-y) -> (1,0) = 5
    # Top (+y) -> (1,0) = 5
    assert neighbors[idx_7, 0] == idx_6
    assert neighbors[idx_7, 1] == idx_6
    assert neighbors[idx_7, 2] == idx_5
    assert neighbors[idx_7, 3] == idx_5

def test_periodic_amr_connectivity(device="cpu"):
    wp.init()
    basis = Basis(polynomial_degree=1, device=device)
    max_blocks = 32
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=max_blocks)
    
    # Enable Periodicity
    quadtree = Quadtree(device=device, max_blocks=max_blocks, periodic_x=True, periodic_y=True)
    
    # 1. Start with 2x2 Grid (Level 1)
    # BL=4 (0,0), BR=5 (1,0), TL=6 (0,1), TR=7 (1,1)
    quadtree.uniform_refine(1, state, basis)
    
    # 2. Refine Block 5 (BR) to Level 2
    # It has kids (1<<4)|(interleave(2*1+[0,1], 2*0+[0,1]))
    # ix ranges [2,3], iy ranges [0,1]
    # Kids are: (2,0), (3,0), (2,1), (3,1)
    # L2 Codes: 
    # (2,0) -> 16 | (part1(2)|part1(0)<<1) = 16 | 4 = 20
    # (3,0) -> 16 | (part1(3)|part1(0)<<1) = 16 | 5 = 21
    # (2,1) -> 16 | (part1(2)|part1(1)<<1) = 16 | 6 = 22
    # (3,1) -> 16 | (part1(3)|part1(1)<<1) = 16 | 7 = 23
    
    h_codes = quadtree.block_morton_codes.numpy()
    idx_5 = -1
    for i in range(4):
        if h_codes[i] == 5:
            idx_5 = i
            break
    
    quadtree.refine_blocks([idx_5], state, basis)
    
    neighbors = state.neighbors.numpy()
    codes = quadtree.block_morton_codes.numpy()
    levels = quadtree.block_levels.numpy()
    active = state.active_block_indices.numpy()[:quadtree.num_blocks]
    
    # Identify the new L2 blocks
    idx_21 = -1 # Right-most, bottom-most L2 block (3,0)
    for idx in active:
        if codes[idx] == 21:
            idx_21 = idx
            break
            
    assert idx_21 != -1
    assert levels[idx_21] == 2
    
    # Block 21 is at ix=3, iy=0 (Level 2).
    # Its RIGHT neighbor in periodic X should be ix=0, iy=0 (Level 1 or 2).
    # In our case, ix=0, iy=0 is Block 4 (Level 1).
    # Since it's L2 looking at L1, it should be a MORTAR interface (-2).
    
    assert neighbors[idx_21, 1] == -2, f"Expected MORTAR_FLAG (-2) for periodic right neighbor of L2 block, got {neighbors[idx_21, 1]}"
    
    # Also check the L1 block 4 (0,0) looking LEFT.
    # Its LEFT neighbor is ix=3, iy=0 (Level 2).
    # Should also be MORTAR_FLAG (-2) because the neighbor is finer.
    idx_4 = -1
    for idx in active:
        if codes[idx] == 4:
            idx_4 = idx
            break
    assert idx_4 != -1
    assert neighbors[idx_4, 0] == -2, f"Expected MORTAR_FLAG (-2) for periodic left neighbor of L1 block, got {neighbors[idx_4, 0]}"

if __name__ == "__main__":
    test_periodic_connectivity_2x2()
