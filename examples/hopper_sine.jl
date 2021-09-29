const ContactImplicitMPC = Main
include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
s_sim = get_simulation("hopper_2D", "sine2_2D_lc", "sinusoidal")
s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
model = s.model
env = s.env

# get trajectory
ref_traj = deepcopy(ContactImplicitMPC.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_forward.jld2"),
    load_type = :joint_traj))
# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 100*H*N_sample #500*H*N_sample

# barrier parameter
κ_mpc = 2.0e-4

obj = TrackingObjective(model, env, H_mpc,
    q = [Diagonal(1.0e-1 * [0.1,3,1,3])   for t = 1:H_mpc],
    u = [Diagonal(1.0e-0 * [1e-3, 1e0]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
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
        altitude_impact_threshold = 0.05,
        # altitude_verbose = true,
        )
    )

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = ContactImplicitMPC.simulator(s_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactImplicitMPC.InteriorPointOptions(
		γ_reg = 0.0,
		undercut = Inf,
		r_tol = 1.0e-8,
		κ_tol = 1.0e-8,),
    sim_opts = ContactImplicitMPC.SimulatorOptions(warmstart = true))


status = ContactImplicitMPC.simulate!(sim, verbose = true)

################################################################################
# Timing result
################################################################################
process!(sim)
# Time budget
ref_traj.h
# Time used on average
sim.stats.μ_dt
# Speed ratio
H_sim * h_sim / sum(sim.stats.dt)



plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)

plot_lines!(vis, model, sim.traj.q[1:10:end])
plot_surface!(vis, s_sim.env, n=200)
anim = visualize_robot!(vis, model, sim.traj, sample=5)
anim = visualize_force!(vis, model, s_sim.env, sim.traj, anim=anim, h=h_sim, sample = 5)

plot(hcat([u[[1,2]] for u in Vector.(sim.traj.u)[1:end]]...)')


# filename = "hopper_ipopt"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)

# const ContactImplicitMPC = Main
