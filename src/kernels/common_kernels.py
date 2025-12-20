import warp as wp

@wp.kernel
def rk_stage_1(
    q: wp.array(dtype=wp.vec4, ndim=2),
    rhs: wp.array(dtype=wp.vec4, ndim=2),
    dt: wp.float32,
    q_out: wp.array(dtype=wp.vec4, ndim=2)
):
    """
    Stage 1: Q(1) = Q_n + dt * RHS(Q_n)
    """
    e, i = wp.tid()
    q_out[e, i] = q[e, i] + dt * rhs[e, i]

@wp.kernel
def rk_stage_2(
    q: wp.array(dtype=wp.vec4, ndim=2),      # Q_n (Startzustand)
    q_1: wp.array(dtype=wp.vec4, ndim=2),    # Q(1) aus Stage 1
    rhs: wp.array(dtype=wp.vec4, ndim=2),    # RHS(Q(1))
    dt: wp.float32,
    q_out: wp.array(dtype=wp.vec4, ndim=2)   # Ziel für Q(2)
):
    """
    Stage 2: Q(2) = 3/4 * Q_n + 1/4 * (Q(1) + dt * RHS(Q(1)))
    """
    e, i = wp.tid()
    # 0.75 * Q_n + 0.25 * (Q_1 + dt * RHS)
    q_out[e, i] = 0.75 * q[e, i] + 0.25 * (q_1[e, i] + dt * rhs[e, i])

@wp.kernel
def rk_stage_3(
    q: wp.array(dtype=wp.vec4, ndim=2),      # Q_n
    q_2: wp.array(dtype=wp.vec4, ndim=2),    # Q(2) aus Stage 2
    rhs: wp.array(dtype=wp.vec4, ndim=2),    # RHS(Q(2))
    dt: wp.float32,
    q_out: wp.array(dtype=wp.vec4, ndim=2)   # Ziel für Q(n+1) (überschreibt oft self.Q)
    ):
    """
    Stage 3: Q(n+1) = 1/3 * Q_n + 2/3 * (Q(2) + dt * RHS(Q(2)))
    """
    e, i = wp.tid()
    # 1/3 * Q_n + 2/3 * (Q_2 + dt * RHS)
    q_out[e, i] = (1.0/3.0) * q[e, i] + (2.0/3.0) * (q_2[e, i] + dt* rhs[e, i])
