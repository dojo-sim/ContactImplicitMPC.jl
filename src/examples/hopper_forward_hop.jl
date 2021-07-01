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
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_in_place.jld2"),
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

n_opts = NewtonOptions(r_tol=3e-8, max_iter=100, verbose = true)

im_traj_ = ImplicitTraj(ref_traj, s, κ=κ)

core = Newton(s, H, h, ref_traj, im_traj_, obj = obj, opts = n_opts)

q0_dist = deepcopy(ref_traj.q[1] + [-0.1,0.0,0,0])
q1_dist = deepcopy(ref_traj.q[2] + [-0.1,0.0,0,0])
newton_solve!(core, s, im_traj_, ref_traj, q0=q0_dist, q1=q1_dist)


im_traj_.d[1]
im_traj_.ip[1].θ
core.traj.θ[1]

l = 1
lu = 1
plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.([q[1:end] for q in core.traj.q])...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][lu:lu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[3,1], hcat(Vector.(vcat([fill(ref_traj.γ[i][1:nc], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:end] for u in core.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)


implicit_dynamics!(im_traj, s, core.traj, κ = im_traj.ip[1].κ)
norm.(im_traj.d)
im_traj.d
core.traj.θ[1]
core.traj.θ[2]



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

traj = get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_forward.jld2"),
    load_type=:joint_traj)
