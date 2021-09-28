include(joinpath(@__DIR__, "..", "src/dynamics", "flamingo", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)


# get hopper model
s_sim = get_simulation("flamingo", "slope_smooth_2D_lc", "slope")
s = get_simulation("flamingo", "flat_2D_lc", "flat")
model = s.model
env = s.env
const ContactImplicitMPC = Main
ref_traj_ = deepcopy(ContactImplicitMPC.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))
ref_traj = deepcopy(ref_traj_)

H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 15
h_sim = h / N_sample
H_sim = 10000 #35000

# barrier parameter
κ_mpc = 2.0e-4

obj = TrackingVelocityObjective(model, env, H_mpc,
    v = [Diagonal(1e-3 * [1e0,1,1e4,1,1,1,1,1e4,1e4]) for t = 1:H_mpc],
    q = [Diagonal(1e-1 * [3e2, 1e-6, 3e2, 1, 1, 1, 1, 0.1, 0.1]) for t = 1:H_mpc],
    u = [Diagonal(3e-1 * [0.1; 0.1; 0.3; 0.3; ones(model.dim.u-6); 2; 2]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	ip_type = :interior_point,
	mode = :configuration,
	ip_opts = InteriorPointOptions(
					undercut = 5.0,
					κ_tol = κ_mpc,
					r_tol = 1.0e-8,
					diff_sol = true,
					solver = :empty_solver,
					max_time = 1e5,),
    n_opts = NewtonOptions(
        r_tol = 3e-4,
		# verbose = true,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        altitude_update = true,
        altitude_impact_threshold = 0.02,
        # altitude_verbose = true,
        )
    )


q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = simulator(s_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
	ip_opts = InteriorPointOptions(
		undercut = Inf,
		γ_reg = 0.0,
		# verbose = true,
        r_tol = 1.0e-8,
        κ_tol = 1.0e-8),
    sim_opts = SimulatorOptions(warmstart = true),
    )

@time status = simulate!(sim, verbose = true)
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

plot_surface!(vis, s_sim.env, xlims=[-0.5, 10.5], ylims = [-0.5, 0.5])
plot_lines!(vis, s_sim.model, sim.traj.q)
# anim = visualize_robot!(vis, model, sim.traj, sample=10)
anim = visualize_meshrobot!(vis, s_sim.model, sim.traj, sample=10)
anim = visualize_force!(vis, s_sim.model, s_sim.env, sim.traj, anim=anim, h=h_sim, sample=10)

(s_sim.env.surf(5.0) - s_sim.env.surf(4.0)) / (5.0 - 4.0) * 360 / (2π)


# filename = "flamingo_40_steps_sine"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
#
# flamingo_ghost!(vis, sim, model_sim.env.surf)
# flamingo_animation!(vis, sim, model_sim.env.surf)
