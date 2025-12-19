import numpy as np
import warp as wp

class Mesh:
    def __init__(self, nx, ny, x_min=0.0, x_max=1.0, y_min=0.0, y_max=1.0, device="cuda"):
        self.nx = nx
        self.ny = ny
        self.num_elements = nx * ny
        self.device = device
        
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

        self.dx = (x_max - x_min) / nx
        self.dy = (y_max - y_min) / ny

        # --- Mesh Data (Host) ---
        self.vertices_host = np.zeros((self.num_elements, 4, 2), dtype=np.float32)
        self.connectivity_host = -np.ones((self.num_elements, 4, 2), dtype=np.int32)
        
        self._create_mesh()
        self._create_connectivity()

        # --- GEOMETRIE FAKTOREN (WICHTIG!) ---
        # Für kartesische Gitter sind diese konstant, aber wir speichern sie pro Element
        # um später auch verzerrte Gitter unterstützen zu können.
        
        # Ableitungs-Skalierung: dr/dx und ds/dy
        # Wenn r von [-1, 1] geht und x von [0, dx], ist der Faktor 2/dx
        self.rx_host = np.full(self.num_elements, 2.0 / self.dx, dtype=np.float32)
        self.sy_host = np.full(self.num_elements, 2.0 / self.dy, dtype=np.float32)
        
        # Volumen-Skalierung (Determinante der Jacobi-Matrix dx/dr * dy/ds)
        # J = (dx/2) * (dy/2)
        self.J_host  = np.full(self.num_elements, (self.dx * self.dy) / 4.0, dtype=np.float32)
        
        # Oberflächen-Skalierung (J_surf) für die Flux-Terme
        # Für Flux in x-Richtung (linke/rechte Kante) ist die Länge dy -> Faktor dy/2
        # Für Flux in y-Richtung (untere/obere Kante) ist die Länge dx -> Faktor dx/2
        self.Js_x_host = np.full(self.num_elements, self.dy / 2.0, dtype=np.float32) 
        self.Js_y_host = np.full(self.num_elements, self.dx / 2.0, dtype=np.float32)

        # --- Transfer to Warp ---

# --- Transfer to Warp ---
        self.vertices = wp.array(self.vertices_host, dtype=wp.vec2, device=self.device)
        
        # KORREKTUR: Der Kernel erwartet connectivity mit ndim=2 (nur Neighbor-IDs).
        # Wir nehmen daher nur den ersten Teil des Tupels [NeighborID, FaceID] -> [:, :, 0]
        # self.connectivity_host hat Shape (NumElements, 4, 2) -> Slice gibt (NumElements, 4)
        self.connectivity = wp.array(self.connectivity_host[:, :, 0], dtype=wp.int32, device=self.device)
        
        self.rx = wp.array(self.rx_host, dtype=wp.float32, device=self.device)
        self.sy = wp.array(self.sy_host, dtype=wp.float32, device=self.device)
        self.J  = wp.array(self.J_host, dtype=wp.float32, device=self.device)
        self.Js_x = wp.array(self.Js_x_host, dtype=wp.float32, device=self.device)
        self.Js_y = wp.array(self.Js_y_host, dtype=wp.float32, device=self.device)

    def _create_mesh(self):
        # (Dein Code war hier korrekt)
        for j in range(self.ny):
            for i in range(self.nx):
                element_id = j * self.nx + i
                x0 = self.x_min + i * self.dx
                y0 = self.y_min + j * self.dy
                
                self.vertices_host[element_id, 0, :] = [x0, y0]
                self.vertices_host[element_id, 1, :] = [x0 + self.dx, y0]
                self.vertices_host[element_id, 2, :] = [x0 + self.dx, y0 + self.dy]
                self.vertices_host[element_id, 3, :] = [x0, y0 + self.dy]

    def _create_connectivity(self):
        # (Dein Code war hier korrekt - die Logik passt zur Basis-Klasse)
        # Face 0: Bottom, Face 1: Right, Face 2: Top, Face 3: Left
        # Die Face-Indizes passen zur Standard CCW Nummerierung.
        for j in range(self.ny):
            for i in range(self.nx):
                element_id = j * self.nx + i

                # Left (Face 3) -> Neighbor Right (Face 1)
                ni, nj = (i - 1) if i > 0 else (self.nx - 1), j
                self.connectivity_host[element_id, 3] = [nj * self.nx + ni, 1]

                # Right (Face 1) -> Neighbor Left (Face 3)
                ni, nj = (i + 1) if i < self.nx - 1 else 0, j
                self.connectivity_host[element_id, 1] = [nj * self.nx + ni, 3]

                # Bottom (Face 0) -> Neighbor Top (Face 2)
                ni, nj = i, (j - 1) if j > 0 else (self.ny - 1)
                self.connectivity_host[element_id, 0] = [nj * self.nx + ni, 2]

                # Top (Face 2) -> Neighbor Bottom (Face 0)
                ni, nj = i, (j + 1) if j < self.ny - 1 else 0
                self.connectivity_host[element_id, 2] = [nj * self.nx + ni, 0]
