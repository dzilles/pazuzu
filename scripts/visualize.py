import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mc
import matplotlib.animation as animation
import argparse
import os
import glob
import sys

# Add project root to path so we can import from src
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.physics.equations import conservative_to_primitive

def calculate_scalar_field(q, quantity):
    """
    Helper to extract the scalar field (rho, p, etc.) from the Q state.
    """
    num_elements, Np, _ = q.shape
    
    # Transpose q to shape (4, num_elements, Np) then reshape
    q_reshaped = q.transpose(2, 0, 1).reshape(4, -1)
    prim_reshaped = conservative_to_primitive(q_reshaped)
    # Reshape back to (num_elements, Np, 4)
    prim = prim_reshaped.reshape(4, num_elements, Np).transpose(1, 2, 0)

    if quantity == 'rho':
        return np.mean(prim[:, :, 0], axis=1), "Density (rho)"
    elif quantity == 'p':
        return np.mean(prim[:, :, 3], axis=1), "Pressure (p)"
    elif quantity == 'u':
        return np.mean(prim[:, :, 1], axis=1), "Velocity X (u)"
    elif quantity == 'v':
        return np.mean(prim[:, :, 2], axis=1), "Velocity Y (v)"
    elif quantity == 'vorticity':
        # Placeholder (würde Geometrie-Faktoren benötigen, die nicht im Output sind)
        raise NotImplementedError("Vorticity plotting from file is not fully implemented yet (missing metric factors).")
    else:
        raise ValueError(f"Unknown quantity: {quantity}")

def create_video(args):
    """
    Creates an MP4 video or GIF from the step files.
    """
    # 1. Dateien finden und sortieren
    search_path = os.path.join(args.steps_dir, "step_*.npz")
    files = sorted(glob.glob(search_path))
    
    if not files:
        print(f"Keine Dateien gefunden in {search_path}")
        return

    print(f"Erstelle Video aus {len(files)} Dateien...")

    # 2. Setup Plot mit der ersten Datei
    data0 = np.load(files[0])
    q0 = data0['q']
    vertices = data0['vertices']
    num_elements = q0.shape[0]
    
    vals, label = calculate_scalar_field(q0, args.quantity)
    
    # Bestimme globale Min/Max Werte für stabile Farbskala über die Zeit
    # Optional: Man könnte auch über alle Files iterieren, um das globale Min/Max zu finden.
    # Hier nehmen wir dynamische Skalierung oder feste Werte.
    vmin, vmax = vals.min(), vals.max()
    
    fig, ax = plt.subplots(figsize=(10, 8))
    polygons = [vertices[i] for i in range(num_elements)]
    collection = mc.PolyCollection(polygons, array=vals, cmap='viridis', edgecolor='k', lw=0.1)
    
    # Start-Skalierung (wird im Update angepasst, wenn autoscale gewünscht)
    collection.set_clim(vmin, vmax)
    
    ax.add_collection(collection)
    ax.autoscale_view()
    ax.set_aspect('equal')
    
    cbar = plt.colorbar(collection, ax=ax, label=label)
    title = ax.set_title(f"Step 0 - {label}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # 3. Update Funktion
    def update(frame_idx):
        filename = files[frame_idx]
        try:
            data = np.load(filename)
            q = data['q']
            # Wir nehmen an, das Gitter bewegt sich nicht (Euler), daher vertices nicht neu laden
            # Wenn sich das Gitter bewegt (Lagrange/ALE), müssten wir vertices hier updaten.
            
            new_vals, _ = calculate_scalar_field(q, args.quantity)
            
            collection.set_array(new_vals)
            
            # Dynamische Farbskala (optional, sieht oft besser aus bei Explosionen/Wellen)
            # collection.set_clim(new_vals.min(), new_vals.max())
            
            step_num = os.path.basename(filename).split('_')[1].split('.')[0]
            title.set_text(f"Step {step_num} - {label}")
            
            if frame_idx % 10 == 0:
                print(f"Render frame {frame_idx}/{len(files)}")
                
        except Exception as e:
            print(f"Fehler bei Frame {frame_idx}: {e}")
        
        return collection, title

    # 4. Animation erstellen
    ani = animation.FuncAnimation(fig, update, frames=len(files), blit=False)
    
    output_file = args.save if args.save else "simulation_video.mp4"
    
    print(f"Speichere Video als {output_file} ...")
    
    if output_file.endswith('.gif'):
        ani.save(output_file, writer='pillow', fps=args.fps)
    else:
        # Benötigt ffmpeg installiert!
        try:
            ani.save(output_file, writer='ffmpeg', fps=args.fps, dpi=200)
        except Exception as e:
            print(f"Fehler beim Speichern als MP4 (ist ffmpeg installiert?): {e}")
            print("Versuche als GIF zu speichern...")
            ani.save(output_file.replace(".mp4", ".gif"), writer='pillow', fps=args.fps)

    print("Fertig.")

def plot_single_frame(args):
    """
    Plots a single result file.
    """
    try:
        data = np.load(args.file)
    except FileNotFoundError:
        print(f"Error: Output file not found at {args.file}")
        return

    q = data['q']
    vertices = data['vertices']
    num_elements = q.shape[0]
    
    vals, label = calculate_scalar_field(q, args.quantity)

    polygons = [vertices[i] for i in range(num_elements)]
    
    fig, ax = plt.subplots(figsize=(10, 8))
    collection = mc.PolyCollection(polygons, array=vals, cmap='viridis', edgecolor='k', lw=0.1)
    
    ax.add_collection(collection)
    ax.autoscale_view()
    ax.set_aspect('equal')
    
    plt.colorbar(collection, ax=ax, label=label)
    ax.set_title(f"2D DG Solver Result: {label}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    if args.save:
        plt.savefig(args.save, dpi=300)
        print(f"Plot saved to {args.save}")
    else:
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize 2D DG solver results.")
    
    # Modus Auswahl
    parser.add_argument("--video", action="store_true", help="Create a video from all steps in the steps directory.")
    
    # Datei/Pfad Argumente
    parser.add_argument("file", type=str, nargs='?', default="data/output.npz", help="Path to a single .npz file (for single frame mode).")
    parser.add_argument("--steps_dir", type=str, default="data/steps", help="Directory containing step_*.npz files (for video mode).")
    
    # Plot Optionen
    parser.add_argument("-q", "--quantity", type=str, default="rho", choices=['rho', 'p', 'u', 'v', 'vorticity'], help="Quantity to plot.")
    parser.add_argument("-s", "--save", type=str, help="Output filename (e.g., plot.png or video.mp4).")
    parser.add_argument("--fps", type=int, default=10, help="Frames per second for video.")
    
    args = parser.parse_args()
    
    if args.video:
        create_video(args)
    else:
        plot_single_frame(args)