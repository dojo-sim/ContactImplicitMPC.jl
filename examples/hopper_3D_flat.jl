const ContactImplicitMPC = Main
include(joinpath(@__DIR__, "..", "src/dynamics", "hopper_3D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
s = get_simulation("hopper_3D", "flat_3D_lc", "flat")
model = s.model
env = s.env

nq = model.dim.q
nu = model.dim.u
nc = model.dim.c

# get trajectory
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_3D/gaits/gait_forward.jld2"),
    load_type = :joint_traj))

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 10
H_mpc = 20
h_sim = h / N_sample
H_sim = 1000 # 10000 # 12000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingObjective(model, env, H_mpc,
    q = [Diagonal(1.0e-1 * [3,3,0.1,5e+1,5e+1,5e+1,10])   for t = 1:H_mpc],
    u = [Diagonal(1.0e-0 * [1e-1, 1e-1, 1e1]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
        # altitude_update = false,
        altitude_update = true,
        altitude_impact_threshold = 0.05,
        altitude_verbose = true,
        # live_plotting = true,
        )
    )

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	# mode = :configurationforce,
	mode = :configuration,
	ip_type = :mehrotra,
    n_opts = NewtonOptions(
		r_tol = 3e-4,
		max_iter = 5,
		max_time = ref_traj.h, # HARD REAL TIME
		),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        # altitude_update = true,
        # altitude_impact_threshold = 0.05,
        # altitude_verbose = true,
        ),
	ip_opts = InteriorPointOptions(
		max_iter = 100,
		verbose = false,
		r_tol = 1.0e-4,
		κ_tol = 1.0e-4,
		diff_sol = true,
		# κ_reg = 1e-3,
		# γ_reg = 1e-1,
		solver = :empty_solver,
		),
    )

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_tol = 1.0e-6),
    sim_opts = SimulatorOptions(warmstart = true))

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



plot_surface!(vis, env, n=100, xlims=[-0.3, 1.5], ylims=[-1.1, 0.3])
plot_lines!(vis, model, sim.traj.q, offset=0.0)
visualize_robot!(vis, model, sim.traj)


plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)

# visualize_robot!(vis, model, sim.traj, sample=10)



# filename = "hopper_3d_sine"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
