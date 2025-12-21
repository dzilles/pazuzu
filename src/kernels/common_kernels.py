import warp as wp

@wp.kernel
def rk_stage_1(
    q: wp.array(dtype=wp.vec4, ndim=2),
    rhs: wp.array(dtype=wp.vec4, ndim=2),
    dt: wp.float32,
    q_out: wp.array(dtype=wp.vec4, ndim=2)
):
    """
    Performs the first stage of the Low-Storage SSP-RK3 time integration scheme.
    
    Q(1) = Q_n + dt * RHS(Q_n)

    Args:
        q (wp.array): The state vector at time n (Q_n).
        rhs (wp.array): The Right-Hand Side evaluated at Q_n.
        dt (float): Time step size.
        q_out (wp.array): Output buffer for the first intermediate stage Q(1).
    """
    e, i = wp.tid()
    q_out[e, i] = q[e, i] + dt * rhs[e, i]

@wp.kernel
def rk_stage_2(
    q: wp.array(dtype=wp.vec4, ndim=2),      # Q_n (Initial state)
    q_1: wp.array(dtype=wp.vec4, ndim=2),    # Q(1) from Stage 1
    rhs: wp.array(dtype=wp.vec4, ndim=2),    # RHS(Q(1))
    dt: wp.float32,
    q_out: wp.array(dtype=wp.vec4, ndim=2)   # Destination for Q(2)
):
    """
    Performs the second stage of the Low-Storage SSP-RK3 time integration scheme.
    
    Q(2) = 3/4 * Q_n + 1/4 * (Q(1) + dt * RHS(Q(1)))

    Args:
        q (wp.array): The state vector at time n (Q_n).
        q_1 (wp.array): The state vector from the first stage (Q(1)).
        rhs (wp.array): The Right-Hand Side evaluated at Q(1).
        dt (float): Time step size.
        q_out (wp.array): Output buffer for the second intermediate stage Q(2).
    """
    e, i = wp.tid()
    # 0.75 * Q_n + 0.25 * (Q_1 + dt * RHS)
    q_out[e, i] = 0.75 * q[e, i] + 0.25 * (q_1[e, i] + dt * rhs[e, i])

@wp.kernel
def rk_stage_3(
    q: wp.array(dtype=wp.vec4, ndim=2),      # Q_n
    q_2: wp.array(dtype=wp.vec4, ndim=2),    # Q(2) from Stage 2
    rhs: wp.array(dtype=wp.vec4, ndim=2),    # RHS(Q(2))
    dt: wp.float32,
    q_out: wp.array(dtype=wp.vec4, ndim=2)   # Destination for Q(n+1) (often overwrites self.Q)
    ):
    """
    Performs the third and final stage of the Low-Storage SSP-RK3 time integration scheme.
    
    Q(n+1) = 1/3 * Q_n + 2/3 * (Q(2) + dt * RHS(Q(2)))

    Args:
        q (wp.array): The state vector at time n (Q_n).
        q_2 (wp.array): The state vector from the second stage (Q(2)).
        rhs (wp.array): The Right-Hand Side evaluated at Q(2).
        dt (float): Time step size.
        q_out (wp.array): Output buffer for the next time step Q(n+1).
    """
    e, i = wp.tid()
    # 1/3 * Q_n + 2/3 * (Q_2 + dt * RHS)
    q_out[e, i] = (1.0/3.0) * q[e, i] + (2.0/3.0) * (q_2[e, i] + dt* rhs[e, i])
