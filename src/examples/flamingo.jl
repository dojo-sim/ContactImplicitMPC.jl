include(joinpath(@__DIR__, "..", "dynamics", "flamingo", "visuals.jl"))
T = Float64
vis = Visualizer()
render(vis)
open(vis)

# get hopper model
# model_sim = get_model("flamingo", surf="sinusoidal")
model_sim = get_model("flamingo", surf="flat")
model = get_model("flamingo", surf="flat")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd
nz = num_var(model)
nθ = num_data(model)

# get trajectory
# ref_traj = get_trajectory("flamingo", "gait1", load_type=:split_traj_alt, model=model)
# ref_traj = get_trajectory("flamingo", "gait_forward_36_3", load_type=:split_traj_alt, model=model)
ref_traj = get_trajectory("flamingo", "gait_forward_36_4", load_type=:split_traj_alt, model=model)


H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 15
h_sim = h / N_sample
H_sim = 1000#35000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingVelocityObjective(H_mpc, model.dim,
    v = [Diagonal(1e-3 * [1e0,1,1e4,1,1,1,1,1e4,1e4]) for t = 1:H_mpc],
    q = [Diagonal(1e-1 * [3e2, 1e-6, 3e2, 1, 1, 1, 1, 0.1, 0.1]) for t = 1:H_mpc],
    u = [Diagonal(3e-1 * [0.1; 0.1; 0.3; 0.3; ones(nu-6); 2; 2]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, model, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        # altitude_update = true,
        altitude_impact_threshold = 0.02,
        altitude_verbose = true,
        )
    )

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = SimulatorOptions(warmstart = true)
    )

@time status = simulate!(sim)

# save trajectory
@save joinpath(pwd(), "src/dynamics/flamingo/simulations/flat.jld2") sim
@load joinpath(pwd(), "src/dynamics/flamingo/simulations/flat.jld2") sim

l = 9
lu = 1
plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.([q[l:l] for q in sim.traj.q])...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][lu:lu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[3,1], hcat(Vector.(vcat([fill(ref_traj.γ[i][1:nc], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[lu:lu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
# plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)
# plot!(plt[3,1], hcat(Vector.([b[1:nb] for b in sim.traj.b]*N_sample)...)', color=:red, linewidth=1.0)

plot_surface!(vis, sim.model.env, xlims=[-1, 7.5], ylims = [-0.5, 0.5])
anim = visualize_robot!(vis, sim.model, sim.traj, sample=10)
anim = visualize_force!(vis, model_sim, sim.traj, anim=anim, h=h_sim, sample=10)

convert_config(model_sim, sim.traj.q[1])

filename = "flamingo_100_steps"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)


settransform!(vis["/Cameras/default"],
		compose(Translation(0.0, 0.5, -1.0),LinearMap(RotZ(-pi / 2.0))))

# Ghost
flamingo_ghost!(vis, sim, x -> 0.0)

# Animation
flamingo_animation!(vis, sim, x -> 0.0)
