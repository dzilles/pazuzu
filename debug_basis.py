import numpy as np
from scipy.special import legendre

def gauss_lobatto_quadrature(N):
    if N == 0: return np.array([0.0]), np.array([2.0])
    if N == 1: roots = np.array([])
    else: roots = np.roots(legendre(N).deriv(1))
    
    nodes = np.concatenate(([-1.0], np.sort(roots), [1.0]))
    weights = 2 / (N * (N+1) * legendre(N)(nodes)**2)
    return nodes, weights

def gauss_quadrature(Nq):
    from scipy.special import roots_legendre
    nodes, weights = roots_legendre(Nq)
    return nodes, weights

def interpolation_matrix_1d(nodes_from, nodes_to):
    ni = len(nodes_to)
    nj = len(nodes_from)
    Interp = np.zeros((ni, nj))
    for i in range(ni):
        for j in range(nj):
            val = 1.0
            for k in range(nj):
                if j != k:
                    val *= (nodes_to[i] - nodes_from[k]) / (nodes_from[j] - nodes_from[k])
            Interp[i, j] = val
    return Interp

def debug():
    N = 2
    nodes_1d, weights_1d = gauss_lobatto_quadrature(N)
    
    # 1. Consistent Mass Matrix Setup
    # Need to integrate degree 2N = 4. Gauss N+1=3 points is exact for 5.
    nq = N + 1
    nq_nodes, nq_weights = gauss_quadrature(nq)
    
    Vq = interpolation_matrix_1d(nodes_1d, nq_nodes)
    M_cons = Vq.T @ np.diag(nq_weights) @ Vq
    
    # 2. RHS for Left/Right
    nodes_q_left = 0.5 * (nq_nodes - 1.0)
    nodes_q_right = 0.5 * (nq_nodes + 1.0)
    Vq_left = interpolation_matrix_1d(nodes_1d, nodes_q_left)
    Vq_right = interpolation_matrix_1d(nodes_1d, nodes_q_right)
    
    # RHS_L_ik = 0.5 * Integral( phi_i(parent) * phi_k(child) )
    # phi_i(parent) at child Gauss nodes is Vq_left[m, i]
    # phi_k(child) at child Gauss nodes is Vq[m, k]
    RHS_L = 0.5 * Vq_left.T @ np.diag(nq_weights) @ Vq
    RHS_R = 0.5 * Vq_right.T @ np.diag(nq_weights) @ Vq
    
    M_inv = np.linalg.inv(M_cons)
    R_left = M_inv @ RHS_L
    R_right = M_inv @ RHS_R
    
    # 3. Prolongation (same as before)
    nodes_child_left = 0.5 * (nodes_1d - 1.0)
    nodes_child_right = 0.5 * (nodes_1d + 1.0)
    P_left = interpolation_matrix_1d(nodes_1d, nodes_child_left)
    P_right = interpolation_matrix_1d(nodes_1d, nodes_child_right)
    
    # Test Reversibility
    u_parent = nodes_1d**2 + nodes_1d + 1.0
    u_child_left = P_left @ u_parent
    u_child_right = P_right @ u_parent
    
    u_recovered = R_left @ u_child_left + R_right @ u_child_right
    
    print(f"Original: {u_parent}")
    print(f"Recovered: {u_recovered}")
    err = np.max(np.abs(u_parent - u_recovered))
    print(f"Error: {err}")
    
    # Check Identity
    Identity = R_left @ P_left + R_right @ P_right
    print("Identity Check:")
    print(Identity)
    print(f"Is Identity? {np.allclose(Identity, np.eye(N+1))}")
    
    # Check Conservation
    # sum(w_p * u_p) == 0.25 * (sum(w_c * u_c_L) + sum(w_c * u_c_R))?
    # No, Jacobian for 1D is 0.5.
    # sum(w_p * u_p) == 0.5 * (sum(w_c * u_c_L) + sum(w_c * u_c_R))
    int_p = np.sum(u_recovered * weights_1d)
    int_c = 0.5 * (np.sum(u_child_left * weights_1d) + np.sum(u_child_right * weights_1d))
    print(f"Int Parent: {int_p}, Int Children: {int_c}")
    print(f"Conservation Error: {abs(int_p - int_c)}")

if __name__ == "__main__":
    debug()