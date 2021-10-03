include(joinpath(module_dir(), "src", "dynamics", "quadruped", "visuals.jl"))
include(joinpath(module_dir(), "src/controller/newton_structure_solver/methods.jl"))

vis = Visualizer()
open(vis)
const ContactControl = Main

s = get_simulation("quadruped", "flat_2D_lc", "flat")
model = s.model
env = s.env

ref_traj_ = deepcopy(ContactControl.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
    load_type = :split_traj_alt))

ref_traj = deepcopy(ref_traj_)

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 10000 #3000

# barrier parameter
κ_mpc = 1.0e-4

obj_mpc = quadratic_objective(model, H_mpc,
    q = [Diagonal(5e-2 * [1.0; 0.02; 0.25; 0.25 * ones(model.dim.q-3)]) for t = 1:H_mpc+2],
    v = [Diagonal(1.0e-4 * ones(model.dim.q)) for t = 1:H_mpc],
    u = [Diagonal(10e-2 * ones(model.dim.u)) for t = 1:H_mpc-1])
# obj = TrackingObjective(model, env, H_mpc,
#     q = [Diagonal(1e-2 * [1.0; 0.02; 0.25; 0.25 * ones(model.dim.q-3)]) for t = 1:H_mpc],
#     u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
#     γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
#     b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj_mpc,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	mode = :configuration,
	newton_mode = :structure,
    n_opts = NewtonOptions(
        r_tol = 1.0e-5,
		β_init = 1.0e-5,
        # solver = :ldl_solver,
        verbose=true,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions())

p = linearized_mpc_policy(ref_traj, s, obj_mpc,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	# mode = :configurationforce,
	mode = :configuration,
	newton_mode = :structure,
	ip_type = :mehrotra,
    n_opts = NewtonOptions(
		verbose = true,
		solver = :lu_solver,
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
	ip_opts = MehrotraOptions(
		max_iter_inner = 100,
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

sim = ContactControl.simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6,
        diff_sol = true),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

@time status = ContactControl.simulate!(sim, verbose = true)

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


# @elapsed status = ContactControl.simulate!(sim)
# @profiler status = ContactControl.simulate!(sim)

plot_lines!(vis, model, sim.traj.q[1:1:end])
plot_surface!(vis, env, ylims=[0.3, -0.05])
anim = visualize_meshrobot!(vis, model, sim.traj, sample=5)
# anim = visualize_robot!(vis, model, sim.traj, anim=anim)
anim = visualize_force!(vis, model, env, sim.traj, anim=anim, h=h_sim, sample = 5)

# Display ghosts
t_ghosts = [1, 1333, 2666]
mvis_ghosts = []
for (i,t) in enumerate(t_ghosts)
    α = i/(length(t_ghosts)+1)
    name = Symbol("ghost$i")
    mvis = build_meshrobot!(vis, model, name=name, α=α)
    push!(mvis_ghosts, mvis)
end

for (i,t) in enumerate(t_ghosts)
    name = Symbol("ghost$i")
    set_meshrobot!(vis, mvis_ghosts[i], model, sim.traj.q[t], name=name)
end



filename = "quadruped_struct_fail"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)
