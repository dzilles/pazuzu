import gmsh
import sys
import os

def view_mesh(filename):
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    print(f"Initializing Gmsh to view '{filename}'...")
    gmsh.initialize()
    
    try:
        gmsh.open(filename)
        print("Mesh loaded. Launching GUI...")
        print("interaction: Use mouse to rotate/pan. Close window to exit.")
        gmsh.fltk.run()
    except Exception as e:
        print(f"Error: {e}")
    finally:
        gmsh.finalize()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        target_file = sys.argv[1]
    else:
        # Default fallback for convenience
        target_file = os.path.join(
            os.path.dirname(__file__), 
            "../tests/verification/cylinder_flow/cylinder.msh"
        )
        target_file = os.path.normpath(target_file)

    view_mesh(target_file)
