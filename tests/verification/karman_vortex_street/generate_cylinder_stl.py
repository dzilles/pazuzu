import numpy as np
import meshio
import os

def generate_cylinder_stl(filename="cylinder.stl", radius=0.5, height=2.0, num_segments=64):
    """
    Generates a simple triangulated cylinder centered at (0,0) extending along Z.
    """
    
    # 1. Generate Vertices
    # Top circle (z = +h/2), Bottom circle (z = -h/2)
    z_top = height / 2.0
    z_bot = -height / 2.0
    
    theta = np.linspace(0, 2*np.pi, num_segments, endpoint=False)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    
    # Vertices: [Top Circle ..., Bottom Circle ...]
    # Top: 0 to N-1
    # Bot: N to 2N-1
    
    top_verts = np.column_stack((x, y, np.full_like(x, z_top)))
    bot_verts = np.column_stack((x, y, np.full_like(x, z_bot)))
    
    vertices = np.vstack((top_verts, bot_verts))
    
    # 2. Generate Faces (Side walls)
    # 2 triangles per segment
    # Quad: Top[i], Top[i+1], Bot[i+1], Bot[i]
    
    faces = []
    N = num_segments
    
    for i in range(N):
        i_next = (i + 1) % N
        
        # Indices
        t1 = i
        t2 = i_next
        b1 = i + N
        b2 = i_next + N
        
        # Triangle 1: t1, b1, t2
        faces.append([t1, b1, t2])
        
        # Triangle 2: t2, b1, b2
        faces.append([t2, b1, b2])
        
    # Optional: Caps (Top and Bottom) - Not strictly needed for 2D slice if it's open, 
    # but good for completeness. The SDF kernel might rely on "inside/outside" which requires a closed manifold.
    # Let's add center points for caps to make it a closed volume.
    
    # Add Center Top (Index 2N) and Center Bot (Index 2N+1)
    vertices = np.vstack((vertices, [[0,0,z_top], [0,0,z_bot]]))
    c_top = 2*N
    c_bot = 2*N + 1
    
    for i in range(N):
        i_next = (i + 1) % N
        
        # Top Cap: Center, i, i_next (CCW viewed from top)
        faces.append([c_top, i, i_next])
        
        # Bottom Cap: Center, i_next, i (CW viewed from bottom -> CCW normal out)
        # Indices for bottom circle are i+N
        faces.append([c_bot, i_next+N, i+N])

    faces = np.array(faces)
    
    # 3. Write STL
    cells = [("triangle", faces)]
    mesh = meshio.Mesh(vertices, cells)
    
    print(f"Writing {filename} with {len(vertices)} vertices and {len(faces)} faces...")
    mesh.write(filename)
    print("Done.")

if __name__ == "__main__":
    output_path = os.path.join(os.path.dirname(__file__), "cylinder.stl")
    generate_cylinder_stl(output_path, radius=0.5, num_segments=100)
