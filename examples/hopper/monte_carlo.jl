# PREAMBLE

# PKG_SETUP

# ## Setup

using ContactImplicitMPC 
using LinearAlgebra 
using StaticArrays
using Random

# ## Simulation
s = get_simulation("hopper_2D", "flat_2D_lc", "flat");
model = s.model
env = s.env
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c

# ## Run Policy
function run_policy(s::Simulation; H_sim::Int = 4000, verbose = verbose,
		offset::AbstractVector = zeros(s.model.dim.q))
	model = s.model
	env = s.env
	ref_traj = deepcopy(get_trajectory(model, env,
	    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_in_place.jld2"),
	    load_type = :joint_traj))
	
	H = ref_traj.H
	h = ref_traj.h
	N_sample = 5
	H_mpc = 10
	h_sim = h / N_sample
	κ_mpc = 1.0e-4

	obj = TrackingVelocityObjective(model, env, H_mpc,
		q = [[Diagonal(1.0e-2 * [0.1,3,1,3])   for t = 1:H_mpc-2]; [Diagonal(1.0e-1 * [0.1,3,1,3])   for t = 1:2]],
		v = [[Diagonal(1.0e+1 * [0.1,3,1,3])   for t = 1:2]; [Diagonal(1.0e-3 * [0.1,3,1,3])   for t = 1:H_mpc-2]],
	    u = [Diagonal(1.0e-0 * [3e-3, 1e0]) for t = 1:H_mpc],
	    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
	    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

	p = ci_mpc_policy(ref_traj, s, obj,
	    H_mpc = H_mpc,
	    N_sample = N_sample,
	    κ_mpc = κ_mpc,
		mode = :configuration,
	    n_opts = NewtonOptions(
			r_tol = 3e-4,
			max_iter = 5,
			max_time = ref_traj.h, # HARD REAL TIME
			),
	    mpc_opts = CIMPCOptions(
	        ),
		ip_opts = InteriorPointOptions(
			undercut = 5.0,
			κ_tol = κ_mpc,
			r_tol = 1.0e-8,
			diff_sol = true,
			solver = :empty_solver,
			max_time = 1e5),
	    )

	q1_ref = copy(ref_traj.q[2]) + offset
	q0_ref = copy(ref_traj.q[1]) + offset
	q1_sim = SVector{model.dim.q}(q1_ref)
	q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
	@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

	sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
	    p = p,
	    ip_opts = InteriorPointOptions(
			γ_reg = 0.0,
			undercut = Inf,
			r_tol = 1.0e-8,
			κ_tol = 1.0e-8,),
	    sim_opts = SimulatorOptions(warmstart = true))

	status = simulate!(sim, verbose = false)
	return deepcopy(sim.traj)
end

# ## Collect runs
function collect_runs(s::Simulation; n::Int = 1, H_sim::Int = 4000, verbose::Bool = false)
	trajs = []
	Random.seed!(100)
	offset_max = [+0.50, 0.30, +0.30, +0.00]
	offset_min = [-0.50, 0.00, -0.30, -0.30]
	offset_Δ = offset_max .- offset_min
	for i = 1:n
		verbose && println("sample = $i/$n")
		offset = offset_min .+ offset_Δ .* rand(s.model.dim.q)
		traj = run_policy(s, H_sim = H_sim, offset = offset, verbose = verbose)
		push!(trajs, traj)
	end
	return trajs
end

# ## Visualize runs
function visualize_runs!(vis::ContactImplicitMPC.Visualizer, model::ContactModel, trajs::AbstractVector;
		sample=max(1, Int(floor(trajs[1].H / 100))), h=trajs[1].h*sample,  α=1.0,
		anim::ContactImplicitMPC.MeshCat.Animation=ContactImplicitMPC.MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=ContactImplicitMPC.model_name(model))

	for (i, traj) in enumerate(trajs)
		name_i = Symbol(string(name) * string(i))
		anim = visualize_robot!(vis, model, traj, sample = sample, h = h, α = α, anim = anim, name = name_i)
	end
	return anim
end

# ## Experiment
H_sim = 1000
trajs = collect_runs(s, n = 100, H_sim = H_sim, verbose = true);

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
anim = visualize_runs!(vis, s.model, trajs, α=0.1, sample=5)



