include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
model = get_model("hopper_2D")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
ref_traj = get_trajectory("hopper_2D", "gait_forward", load_type=:joint_traj)
H = ref_traj.H
h = ref_traj.h
κ = 1.0e-3

# Cost function
cost = CostFunction(H, model.dim,
    q = [[Diagonal(1.0e-5 * ones(nq)) for t = 1:H-1]; [Diagonal(1.0e3 * ones(nq))]],
    u = [Diagonal(3.0e-2 * [0.1, 1.0]) for t = 1:H],
    γ = [Diagonal(1.0e-100 * ones(nc)) for t = 1:H],
    b = [Diagonal(1.0e-100 * ones(nb)) for t = 1:H])
n_opts = NewtonOptions(r_tol=3e-6, κ_init=κ, κ_tol=2κ, solver_inner_iter=100)
core = Newton(H, h, model, cost = cost, opts = n_opts)
# im_traj = ImplicitTraj(ref_traj, model, κ=κ)
# q0_dist = deepcopy(ref_traj.q[1] + [-0.1,0.0,0,0])
# q1_dist = deepcopy(ref_traj.q[2] + [-0.1,0.0,0,0])
# newton_solve!(core, model, im_traj, ref_traj, q0=q0_dist, q1=q1_dist, verbose=true)

visualize!(vis, model, ref_traj.q, Δt=ref_traj.h)

plot(hcat(Vector.(ref_traj.u)...)')
# Check ~perfect loop
@show round.(ref_traj.q[end-2] - ref_traj.q[1], digits=3)
@show round.(ref_traj.q[end-1] - ref_traj.q[2], digits=3)
