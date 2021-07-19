using Random
using LinearAlgebra
const ContactControl = Main

################################################################################
# Benchmark Utils
################################################################################

function get_initialization(ref_traj::ContactTraj{T, nq}, t::Int; z_dis = 1e-2, θ_dis = 1e-2) where {T,nq}
	Random.seed!(10)
	z = 1e-1ones(length(ref_traj.z[t]))
	z[1:nq] = deepcopy(ref_traj.θ[t][nq+1:2nq])

	θ = deepcopy(ref_traj.θ[t])
	θ[1:end - 2] += θ_dis * ones(length(θ) - 2)
	return z, θ
end

function get_mehrotra_expanded_solver(s::Simulation, z, θ;
		r_tol = 1e-8, κ_tol = 1e-8, linear = false)
	model = s.model
	env = s.env
	lin = LinearizedStep(s, z, θ, 0.0) # /!| here we do not use κ from ref_traj
	@warn "wrong mehrotra solver"
	if linear
		ip = mehrotra(z, θ,
			ix = linearization_var_index(model, env)[1],
			iy1 = linearization_var_index(model, env)[2],
			iy2 = linearization_var_index(model, env)[3],
			idyn = linearization_term_index(model, env)[1],
			irst = linearization_term_index(model, env)[2],
			ibil = linearization_term_index(model, env)[3],
			idx_ineq = inequality_indices(model, env),
			idx_soc = soc_indices(model, env),
			r!  = r!,
			rm! = rm!,
			rz! = rz!,
			rθ! = rθ!,
			r  = RLin(s, lin.z, lin.θ, lin.r, lin.rz, lin.rθ),
			rz = RZLin(s, lin.rz),
			rθ = RθLin(s, lin.rθ),
			v_pr = view(zeros(1,1), 1,1),
			v_du = view(zeros(1,1), 1,1),
			# opts = Mehrotra12Options(
			opts = MehrotraOptions(
				max_iter_inner = 100,
				r_tol = r_tol,
				κ_tol = κ_tol,
				solver = :empty_solver,
				))
	else
		ip = mehrotra(z, θ,
			ix = linearization_var_index(model, env)[1],
			iy1 = linearization_var_index(model, env)[2],
			iy2 = linearization_var_index(model, env)[3],
			idyn = linearization_term_index(model, env)[1],
			irst = linearization_term_index(model, env)[2],
			ibil = linearization_term_index(model, env)[3],
		    idx_ineq = inequality_indices(model, env),
		    idx_soc = soc_indices(model, env),
			r! = s.res.r!,
		    rm! = s.res.rm!,
		    rz! = s.res.rz!,
		    rθ! = s.res.rθ!,
		    rz = s.rz,
		    rθ = s.rθ,
			# opts = Mehrotra12Options(
		    opts = MehrotraOptions(
		        max_iter_inner = 100,
		        r_tol = r_tol,
		        κ_tol = κ_tol,
				))
	end
	return ip
end

function get_interior_point_solver(s::Simulation, z, θ;
		r_tol = 1e-8, κ_tol = 1e-8, linear = false)
	model = s.model
	env = s.env
	lin = LinearizedStep(s, z, θ, 0.0) # /!| here we do not use κ from ref_traj
	if linear
		ip = interior_point(z, θ,
			ix = linearization_var_index(model, env)[1],
			iy1 = linearization_var_index(model, env)[2],
			iy2 = linearization_var_index(model, env)[3],
			idyn = linearization_term_index(model, env)[1],
			irst = linearization_term_index(model, env)[2],
			ibil = linearization_term_index(model, env)[3],
			idx_ineq = inequality_indices(model, env),
			idx_soc = soc_indices(model, env),
			r!  = r!,
			rm! = rm!,
			rz! = rz!,
			rθ! = rθ!,
			r  = RLin(s, lin.z, lin.θ, lin.r, lin.rz, lin.rθ),
			rz = RZLin(s, lin.rz),
			rθ = RθLin(s, lin.rθ),
			v_pr = view(zeros(1,1), 1,1),
			v_du = view(zeros(1,1), 1,1),
			opts = InteriorPointOptions(
				max_iter_inner = 100,
				r_tol = r_tol,
				κ_tol = κ_tol,
				solver = :empty_solver,
				))
	else
		ip = interior_point(z, θ,
		    idx_ineq = inequality_indices(model, env),
		    idx_soc = soc_indices(model, env),
		    r! = s.res.r!,
		    rm! = s.res.rm!,
		    rz! = s.res.rz!,
		    rθ! = s.res.rθ!,
		    rz = s.rz,
		    rθ = s.rθ,
		    opts = InteriorPointOptions(
				max_iter_inner = 100,
				r_tol = r_tol,
				κ_init = κ_tol,
				κ_tol = 2κ_tol,
				))
	end
	return ip
end

################################################################################
# Benchmark Structures & Methods
################################################################################

# Benchmark
@with_kw mutable struct BenchmarkOptions12{T}
	r_tol::T = 1.0e-8            # residual tolerance
	κ_tol::T = 1.0e-8            # bilinear tolerance
	z_dis::T = 1.0e-4            # disturbance on initial z
    θ_dis::T = 1.0e-4            # disturbance on initial θ
    linear::Bool = false         # whether to use or not the linearized residuals
end

mutable struct BenchmarkStatistics13{T}
	failure::AbstractVector{Bool}
	iter::AbstractVector{Int}
	r_vio::AbstractVector{T}
	κ_vio::AbstractVector{T}
end

function BenchmarkStatistics13()
	failure = trues(0)
	iter = zeros(Int, 0)
	r_vio = zeros(0)
	κ_vio = zeros(0)
	return BenchmarkStatistics13(failure, iter, r_vio, κ_vio)
end

function record!(stats::BenchmarkStatistics13, s, i, r, κ)
	push!(stats.failure, s)
	push!(stats.iter, i)
	push!(stats.r_vio, r)
	push!(stats.κ_vio, κ)
	return nothing
end

function benchmark_mehrotra_expanded!(stats::BenchmarkStatistics13, s::Simulation,
		ref_traj::ContactTraj; opts::BenchmarkOptions12)

	nz = num_var(s.model, s.env)
	for t = 1:ref_traj.H
		z, θ = get_initialization(ref_traj, t,
			z_dis = opts.z_dis,
			θ_dis = opts.θ_dis)
		ip = get_mehrotra_expanded_solver(s, deepcopy(ref_traj.z[t]), deepcopy(ref_traj.θ[t]),
			r_tol = opts.r_tol,
			κ_tol = opts.κ_tol,
			linear = opts.linear)
		interior_point_solve!(ip, z, θ)
		rv = residual_violation(ip, ip.r)
		κv = bilinear_violation(ip, ip.r)
		fl = !((rv < opts.r_tol) && (κv < opts.κ_tol))
		record!(stats, fl, ip.iterations, rv, κv)
	end
	return nothing
end

function benchmark_interior_point!(stats::BenchmarkStatistics13, s::Simulation,
		ref_traj::ContactTraj; opts::BenchmarkOptions12)

	nz = num_var(s.model, s.env)
	for t = 1:ref_traj.H
		z, θ = get_initialization(ref_traj, t,
			z_dis = opts.z_dis,
			θ_dis = opts.θ_dis)
		ip = get_interior_point_solver(s, deepcopy(ref_traj.z[t]), deepcopy(ref_traj.θ[t]),
			r_tol = opts.r_tol,
			κ_tol = opts.κ_tol,
			linear = opts.linear)
		interior_point_solve!(ip, z, θ)
		rv = norm(ip.r, Inf)
		fl = !(rv < opts.r_tol)
		record!(stats, fl, ip.iterations, rv, NaN)
	end
	return nothing
end

function benchmark!(stats::BenchmarkStatistics13, s::Simulation,
		ref_traj::ContactTraj; opts::BenchmarkOptions12, solver::Symbol)
	if solver == :interior_point
		benchmark_interior_point!(stats, s, ref_traj; opts = opts)
	elseif solver == :mehrotra_expanded
		benchmark_mehrotra_expanded!(stats, s, ref_traj; opts = opts)
	else
		@warn "Unknown solver"
	end
	return nothing
end

function benchmark(s::AbstractVector{Simulation},
		ref_traj::AbstractVector{<:ContactTraj}, dist::AbstractVector{<:Real};
		opts = BenchmarkOptions12(), solver::Symbol)

	stats = Vector{BenchmarkStatistics13}()
	nsim = length(s)
	ndis = length(dist)
	@assert length(ref_traj) == nsim

	for i = 1:ndis
		sta = BenchmarkStatistics13()
		opts.z_dis = dist[i]
		opts.θ_dis = dist[i]
		for j = 1:nsim
			benchmark!(sta, s[j], ref_traj[j]; opts = opts, solver = solver)
		end
		push!(stats, sta)
	end
	return stats
end


################################################################################
# Loading Problems
################################################################################

# s_quadruped = ContactControl.get_simulation("quadruped", "flat_2D_lc", "flat")
# s_flamingo = ContactControl.get_simulation("flamingo", "flat_2D_lc", "flat")
# s_hopper = ContactControl.get_simulation("hopper_2D", "flat_2D_lc", "flat")
# ref_traj_quadruped = deepcopy(ContactControl.get_trajectory(s_quadruped.model, s_quadruped.env,
#     joinpath(ContactControl.module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
#     load_type = :split_traj_alt))
# ref_traj_flamingo = deepcopy(ContactControl.get_trajectory(s_flamingo.model, s_flamingo.env,
#     joinpath(ContactControl.module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
#     load_type = :split_traj_alt))
# ref_traj_hopper = deepcopy(ContactControl.get_trajectory(s_hopper.model, s_hopper.env,
#     joinpath(ContactControl.module_dir(), "src/dynamics/hopper_2D/gaits/gait_forward.jld2"),
#     load_type = :joint_traj))


################################################################################
# Running Benchmark
################################################################################

s = [s_quadruped, s_flamingo, s_hopper]
ref_traj = [ref_traj_quadruped, ref_traj_flamingo, ref_traj_hopper]
dist = [1e-6, 1e-5, 1e-4, 1e-4, 1e-3, 1e-2]
tol = 1e-8
opts_nonlin = BenchmarkOptions12(r_tol = tol, κ_tol = tol, linear = false)
opts_lin = BenchmarkOptions12(r_tol = tol, κ_tol = tol, linear = true)
stats = Dict{Symbol, AbstractVector{BenchmarkStatistics13}}()
keys = [
	:interior_point_lin,
	:interior_point_nonlin,
	:mehrotra_expanded_lin,
	:mehrotra_expanded_nonlin,
	]

stats[:mehrotra_expanded_nonlin] = benchmark(
	s,
	deepcopy.(ref_traj),
	dist,
	opts = opts_nonlin,
	solver = :mehrotra_expanded,
	)

stats[:interior_point_nonlin] = benchmark(
	s,
	deepcopy.(ref_traj),
	dist,
	opts = opts_nonlin,
	solver = :interior_point,
	)

stats[:mehrotra_expanded_lin] = benchmark(
	s,
	deepcopy.(ref_traj),
	dist,
	opts = opts_lin,
	solver = :mehrotra_expanded,
	)

stats[:interior_point_lin] = benchmark(
	s,
	deepcopy.(ref_traj),
	dist,
	opts = opts_lin,
	solver = :interior_point,
	)

function display_statistics!(plt, stats::AbstractVector{BenchmarkStatistics13}, dist;
		field::Symbol = :iter, solver::String = "solver")
	data = [mean(getfield(s, field)) for s in stats]
	linestyle = occursin("nonlin", solver) ? :solid : :dot
	linewidth = occursin("nonlin", solver) ? 3.0 : 5.0
	marker = occursin("nonlin", solver) ? :circ : :utriangle
	color = occursin("mehrotra_expanded", solver) ? :blue : :red

	plot!(plt, dist, data,
		title = "solver comparison",
		xlabel = "disturbances",
		ylabel = String(field),
		label = solver,
		xlims = [minimum(dist), maximum(dist)*100],
		ylims = [0, Inf],
		xaxis = :log,
		legend = :topright,
		linestyle = linestyle,
		color = color,
		linewidth = linewidth)
	scatter!(plt, dist, data,
		label = false,
		marker = marker,
		markercolor = color,
		markeralpha = 0.6,
		markerstrokecolor = :black,
		markerstrokewidth = 2.0,
		color = color,
		markersize = 8.0,
		)
	display(plt)
	return nothing
end


################################################################################
# Benchmark Results
################################################################################

plt = plot()
for k in keys
	display_statistics!(plt, stats[k], dist, field = :iter, solver = String(k))
end
plt = plot()
for k in keys
	display_statistics!(plt, stats[k], dist, field = :failure,  solver = String(k))
end
