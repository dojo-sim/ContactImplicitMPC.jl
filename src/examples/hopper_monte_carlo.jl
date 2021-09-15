const ContactControl = Main
include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)


# get hopper model
s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
model = s.model
env = s.env
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c


function run_policy(s::Simulation; H_sim::Int = 4000, verbose = verbose,
		offset::AbstractVector = zeros(s.model.dim.q))
	model = s.model
	env = s.env
	ref_traj = deepcopy(get_trajectory(model, env,
	    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_in_place.jld2"),
	    load_type = :joint_traj))
	# time
	H = ref_traj.H
	h = ref_traj.h
	N_sample = 5
	H_mpc = 10
	h_sim = h / N_sample

	# barrier parameter
	κ_mpc = 1.0e-4

	obj = TrackingVelocityObjective(model, env, H_mpc,
		q = [[Diagonal(1.0e-2 * [0.1,3,1,3])   for t = 1:H_mpc-2]; [Diagonal(1.0e-1 * [0.1,3,1,3])   for t = 1:2]],
		v = [[Diagonal(1.0e+1 * [0.1,3,1,3])   for t = 1:2]; [Diagonal(1.0e-3 * [0.1,3,1,3])   for t = 1:H_mpc-2]],
	    u = [Diagonal(1.0e-0 * [3e-3, 1e0]) for t = 1:H_mpc],
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

	q1_ref = copy(ref_traj.q[2]) + offset
	q0_ref = copy(ref_traj.q[1]) + offset
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

function visualize_runs!(vis::Visualizer, model::ContactModel, trajs::AbstractVector;
		sample=max(1, Int(floor(trajs[1].H / 100))), h=trajs[1].h*sample,  α=1.0,
		anim::MeshCat.Animation=MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=model_name(model))

	for (i, traj) in enumerate(trajs)
		name_i = Symbol(string(name) * string(i))
		anim = visualize_robot!(vis, model, traj, sample = sample, h = h, α = α, anim = anim, name = name_i)
	end
	return anim
end


H_sim = 5000
trajs = collect_runs(s, n = 100, H_sim = H_sim, verbose = true)
anim = visualize_runs!(vis, s.model, trajs, α=0.1, sample=5)
rep_ref_traj = repeat_ref_traj(ref_traj, Int(ceil(H_sim/N_sample/H)))
anim = visualize_robot!(vis, model, rep_ref_traj, name=:ref, anim=anim, sample=1, α=1.0)
plot_surface!(vis, s.env, xlims=[-2,2], ylims=[-0.2,0.2])



plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)

rep_ref_traj = repeat_ref_traj(ref_traj, Int(ceil(H_sim/N_sample/H)))
plot_surface!(vis, env, n=200)
plot_lines!(vis, model, sim.traj.q[1:10:end])
anim = visualize_robot!(vis, model, sim.traj, sample=5)
anim = visualize_robot!(vis, model, rep_ref_traj, name=:ref, anim=anim, sample=1, α=0.3)

plot(hcat([u[[1,2]] for u in Vector.(sim.traj.u)[1:end]]...)')


filename = "quadruped_slow_drop"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)

# const ContactControl = Main
