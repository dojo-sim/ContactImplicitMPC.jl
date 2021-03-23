
# model = get_model("quadruped")
κ = 1e-4
ref_traj0 = get_trajectory("quadruped", "gait1")
ref_traj0.κ .= κ
H = ref_traj0.H
h = 0.1
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq+nc+nb
nr = nq+nu+nc+nb+nd

# Test Jacobian!
cost0 = CostFunction(H, model.dim,
    Qq=fill(Diagonal(1e-2*SizedVector{nq}([0.02,0.02,1,.15,.15,.15,.15,.15,.15,.15,.15,])), H),
    Qu=fill(Diagonal(3e-2*ones(SizedVector{nu})), H),
    Qγ=fill(Diagonal(1e-6*ones(SizedVector{nc})), H),
    Qb=fill(Diagonal(1e-6*ones(SizedVector{nb})), H),
    )
n_opts0 = NewtonOptions(r_tol=1e-5, solver_inner_iter=1)
core0 = Newton(H, h, model, cost=cost0, n_opts=n_opts0)
impl0 = ImplicitTraj(H, model)
linearization!(model, ref_traj0, impl0)

# The solve initialized with the correct tra, is solved in 1 step.
newton_solve!(model, core0, impl0, ref_traj0)
implicit_dynamics!(model, core0.traj, impl0; κ=core0.traj.κ)
residual!(model, core0, core0.r, core0.ν, impl0, core0.traj, ref_traj0)
@test norm(core0.r.r, 1) / length(core0.r.r) < 1e-5


# The solver can recover from initial disturbances
off = [0.02; 0.02; zeros(nq-2)]
core0.n_opts.solver_inner_iter = 20
newton_solve!(model, core0, impl0, ref_traj0, q0=ref_traj0.q[1]+off, q1=ref_traj0.q[2]+off)
implicit_dynamics!(model, core0.traj, impl0; κ=core0.traj.κ)
residual!(model, core0, core0.r, core0.ν, impl0, core0.traj, ref_traj0)
@test norm(core0.r.r, 1) / length(core0.r.r) < 1e-5
