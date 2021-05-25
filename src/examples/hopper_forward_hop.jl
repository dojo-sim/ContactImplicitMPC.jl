include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
nq = s.model.dim.q
nu = s.model.dim.u
nc = s.model.dim.c
nb = nc * friction_dim(s.env)
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
ref_traj = get_trajectory(s.model, s.env,
    joinpath(pwd(), "src/dynamics/hopper_2D/gaits/gait_in_place.jld2"),
    load_type=:joint_traj)

H = ref_traj.H
h = ref_traj.h
κ = 1.0e-8

# Cost function
obj = TrackingObjective(s.model, s.env, H,
    q = [[Diagonal(1.0e-5 * ones(nq)) for t = 1:H-1]; [Diagonal(1.0e9 * ones(nq))]],
    u = [Diagonal(1.0e-5 * [0.1, 1.0]) for t = 1:H],
    γ = [Diagonal(1.0e-100 * ones(nc)) for t = 1:H],
    b = [Diagonal(1.0e-100 * ones(nb)) for t = 1:H])

n_opts = NewtonOptions(r_tol=3e-8, max_iter=100)

im_traj = ImplicitTraj(ref_traj, s, κ=κ)

core = Newton(s, H, h, ref_traj, im_traj, obj = obj, opts = n_opts)

q0_dist = deepcopy(ref_traj.q[1] + [-0.1,0.0,0,0])
q1_dist = deepcopy(ref_traj.q[2] + [-0.1,0.0,0,0])
newton_solve!(core, s, im_traj, ref_traj, q0=q0_dist, q1=q1_dist, verbose=true)

visualize_robot!(vis, s.model, core.traj.q)

plot(hcat(Vector.(ref_traj.u)...)')
plot(hcat(Vector.(core.traj.u)...)')
# Check ~perfect loop
@show round.(core.traj.q[end-2] - ref_traj.q[1], digits=3)
@show round.(core.traj.q[end-1] - ref_traj.q[2], digits=3)

# Save trajectory
traj = deepcopy(core.traj)
gait_path = joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "gaits", "gait_forward.jld2")
@save gait_path traj

# Reload trajectory
res = JLD2.jldopen(gait_path)
loaded_traj = res["traj"]

traj = get_trajectory("hopper_2D", "gait_forward", load_type=:joint_traj)
plot(hcat(Vector.(traj.q)...)')
