# PREAMBLE

# PKG_SETUP

# ## Setup

using ContactImplicitMPC 
using LinearAlgebra 
using StaticArrays
using Random

# ## Simulation
s = get_simulation("quadruped", "flat_2D_lc", "flat")
model = s.model
env = s.env

# ## Run policy
function run_policy(s::Simulation; H_sim::Int = 2000, verbose = false,
		qinit::AbstractVector = zeros(s.model.dim.q))
	model = s.model
	env = s.env

	ref_traj = deepcopy(ContactImplicitMPC.get_trajectory(model, env,
	    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
	    load_type = :split_traj_alt))

	H = ref_traj.H
	h = ref_traj.h
	N_sample = 5
	H_mpc = 10
	h_sim = h / N_sample
	H_sim = H_sim
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
	    n_opts = NewtonOptions(
			r_tol = 3e-4,
			max_iter = 5,
			max_time = ref_traj.h, # HARD REAL TIME
			),
	    mpc_opts = LinearizedMPCOptions(
	        ),
		ip_opts = InteriorPointOptions(
			max_iter = 100,
			verbose = false,
			r_tol = 1.0e-4,
			κ_tol = 1.0e-4,
			diff_sol = true,
			solver = :empty_solver,
			),
	    )

	q1_ref = copy(qinit)
	q0_ref = copy(qinit)
	q1_sim = SVector{model.dim.q}(q1_ref)
	q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
	@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

	sim = ContactImplicitMPC.simulator(s, q0_sim, q1_sim, h_sim, H_sim,
	    p = p,
		ip_opts = InteriorPointOptions(
			γ_reg = 0.0,
			undercut = Inf,
			r_tol = 1.0e-8,
			κ_tol = 1.0e-8,),
	    sim_opts = ContactImplicitMPC.SimulatorOptions(warmstart = true))

	status = ContactImplicitMPC.simulate!(sim, verbose = false)
	return deepcopy(sim.traj)
end

# ## Collect runs
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

# ## Initial configurations
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

# ## Visualize runs
function visualize_runs!(vis::ContactImplicitMPC.Visualizer, model::ContactModel, trajs::AbstractVector;
		sample=max(1, Int(floor(trajs[1].H / 100))), h=trajs[1].h*sample,  α=1.0,
		anim::ContactImplicitMPC.MeshCat.Animation=ContactImplicitMPC.MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=ContactImplicitMPC.model_name(model))

	for (i, traj) in enumerate(trajs)
		name_i = Symbol(string(name) * string(i))
		anim = visualize_meshrobot!(vis, model, traj, sample = sample, h = h, α = α, anim = anim, name = name_i)
	end
	return anim
end

# ## Experiments
H_sim = 1000
trajs = collect_runs(s, n = 100, H_sim = H_sim, verbose = true)

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
open(vis)

# ## Visualize
anim = visualize_runs!(vis, model, trajs, α=0.3, sample=5)

# ## High drop
conf = [-0.05, 0.6, 0.6, 0.6, 0.0, 0.30]
qinit = initial_configuration(s.model, conf...)
traj = run_policy(s, H_sim = 300, qinit=qinit)
anim = visualize_meshrobot!(vis, model, traj, name=:drop, anim=anim, sample=5, α=1.0)

