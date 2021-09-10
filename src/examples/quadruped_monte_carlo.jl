include(joinpath(@__DIR__, "..", "dynamics", "quadruped", "visuals.jl"))
vis = Visualizer()
open(vis)
const ContactControl = Main

s = get_simulation("quadruped", "flat_2D_lc", "flat")
model = s.model
env = s.env


function run_policy(s::Simulation; H_sim::Int = 2000, verbose = false,
		qinit::AbstractVector = zeros(s.model.dim.q))
	model = s.model
	env = s.env

	ref_traj = deepcopy(ContactControl.get_trajectory(model, env,
	    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
	    load_type = :split_traj_alt))

	# time
	H = ref_traj.H
	h = ref_traj.h
	N_sample = 5
	H_mpc = 10
	h_sim = h / N_sample
	H_sim = H_sim

	# barrier parameter
	κ_mpc = 1.0e-4

	obj = TrackingObjective(model, env, H_mpc,
	    q = [Diagonal(1e-2 * [1.0; 0.02; 0.25; 0.75 * ones(model.dim.q-3)]) for t = 1:H_mpc],
	    u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
		γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
	    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

	p = linearized_mpc_policy(ref_traj, s, obj,
	    H_mpc = H_mpc,
	    N_sample = N_sample,
	    κ_mpc = κ_mpc,
		mode = :configuration,
		ip_type = :mehrotra,
	    n_opts = NewtonOptions(
			r_tol = 3e-4,
			max_iter = 5,
			max_time = ref_traj.h, # HARD REAL TIME
			),
	    mpc_opts = LinearizedMPCOptions(
	        ),
		ip_opts = MehrotraOptions(
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

	q1_ref = copy(qinit)
	q0_ref = copy(qinit)
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

	status = ContactControl.simulate!(sim, verbose = false)
	return deepcopy(sim.traj)
end

function collect_runs(s::Simulation; n::Int = 1, H_sim::Int = 2000, verbose::Bool = false)
	nq = s.model.dim.q
	trajs = []
	Random.seed!(100)
	conf_max = [0.05, 0.8, 0.8, 0.8, +0.2,  0.10]
	conf_min = [0.00, 0.6, 0.6, 0.6, -0.2, -0.30]
	conf_Δ = conf_max .- conf_min
	for i = 1:n
		verbose && println("sample = $i/$n")
		conf = conf_min .+ conf_Δ .* rand(length(conf_min))
		conf[end] = max(conf[end], 0.0)
		qinit = initial_configuration(s.model, conf...)
		traj = run_policy(s, H_sim = H_sim, qinit = qinit, verbose = verbose)
		push!(trajs, traj)
	end
	return trajs
end


function initial_configuration(model::Quadruped, θ0, θ1, θ2, θ3, x, Δz)
    q = zeros(model.dim.q)
	q[1] = x
    q[3] = pi / 2.0
    q[4] = -θ1
    q[5] = θ2

    q[8] = -θ1
    q[9] = θ2

    q[2] = model.l_thigh1 * cos(q[4]) + model.l_calf1 * cos(q[5])

    q[10] = -θ3
    q[11] = acos((q[2] - model.l_thigh2 * cos(q[10])) / model.l_calf2)

    q[6] = -θ3
    q[7] = acos((q[2] - model.l_thigh2 * cos(q[6])) / model.l_calf2)

	q[2] += Δz
	q[3] += θ0
    return q
end

function visualize_runs!(vis::Visualizer, model::ContactModel, trajs::AbstractVector;
		sample=max(1, Int(floor(trajs[1].H / 100))), h=trajs[1].h*sample,  α=1.0,
		anim::MeshCat.Animation=MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=model_name(model))

	for (i, traj) in enumerate(trajs)
		name_i = Symbol(string(name) * string(i))
		anim = visualize_meshrobot!(vis, model, traj, sample = sample, h = h, α = α, anim = anim, name = name_i)
	end
	return anim
end

# monte carlo
H_sim = 1000
trajs = collect_runs(s, n = 100, H_sim = H_sim, verbose = true)
anim = visualize_runs!(vis, s.model, trajs, α=0.3, sample=5)
rep_ref_traj = repeat_ref_traj(ref_traj, Int(ceil(H_sim/N_sample/H)), idx_shift=[1])
anim = visualize_meshrobot!(vis, model, rep_ref_traj, name=:ref, anim=anim, sample=1, α=1.0)
settransform!(vis[:ref], Translation(0.0, -0.4, 0.0))
plot_surface!(vis, s.env, xlims=[-0.5,4], ylims=[-0.4,0.4])

# high drop
conf = [-0.05, 0.6, 0.6, 0.6, 0.0, 0.30]
qinit = initial_configuration(s.model, conf...)
traj = run_policy(s, H_sim = 300, qinit=qinit)
anim = visualize_meshrobot!(vis, model, traj, name=:drop, anim=anim, sample=5, α=1.0)

# filename = "quadruped_monte_carlo_clean"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
