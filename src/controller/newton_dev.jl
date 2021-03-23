model = get_model("quadruped")
κ = 1.0e-4
ref_traj = get_trajectory("quadruped", "gait1")
ref_traj.κ .= κ
H = ref_traj.H
h = 0.1
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# Test Jacobian!
cost = cost_function(H, model.dim,
    q = [Diagonal(1.0e-2 * ([0.02,0.02,1,.15,.15,.15,.15,.15,.15,.15,.15])) for t = 1:H],
    u = [Diagonal(3.0e-2 * ones(nu)) for t = 1:H],
    γ = [Diagonal(1.0e-6 * ones(nc)) for t = 1:H],
    b = [Diagonal(1.0e-6 * ones(nb)) for t = 1:H])
model.dim
opts = NewtonOptions(r_tol = 1.0e-5, solver_inner_iter = 1)
core = Newton(H, h, model, cost = cost, opts = opts)
im_traj = implicit_trajectory(H, model)
linearization!(im_traj, model, ref_traj)

# The solve initialized with the correct traj, is solved in 1 step.
newton_solve!(core, model, im_traj, ref_traj)
implicit_dynamics!(im_traj, model, core.traj, κ = core.traj.κ)
residual!(core.res, model, core, core.ν, im_traj, core.traj, ref_traj)
@test norm(core.res.r, 1) / length(core.res.r) < 1.0e-5

# The solver can recover from initial disturbances
off = [0.02; 0.02; zeros(nq - 2)]
core.opts.solver_inner_iter = 20
newton_solve!(core, model, im_traj, ref_traj,
    q0 = ref_traj.q[1] + off, q1 = ref_traj.q[2] + off)
implicit_dynamics!(im_traj, model, core.traj, κ = core.traj.κ)
residual!(core.res, model, core, core.ν, im_traj, core.traj, ref_traj)
@test norm(core.res.r, 1) / length(core.res.r) < 1.0e-5
