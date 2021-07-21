using Random
using LinearAlgebra
const ContactControl = Main

include(joinpath(module_dir(), "src", "solver", "mehrotra.jl"))
include(joinpath(module_dir(), "src", "solver", "mehrotra_expanded.jl"))
include(joinpath(module_dir(), "src", "solver", "mehrotra_latest.jl"))

################################################################################
# Benchmark Utils
################################################################################

function get_benchmark_initialization(ref_traj::ContactTraj{T, nq}, t::Int; z_dis = 1e-2, θ_dis = 1e-2) where {T,nq}
	Random.seed!(10)
	z = 1e-1ones(length(ref_traj.z[t]))
	z[1:nq] = deepcopy(ref_traj.θ[t][nq+1:2nq])

	θ = deepcopy(ref_traj.θ[t])
	# θ[1:end - 2] += θ_dis * ones(length(θ) - 2)
	θ[1:nq] += θ_dis * ones(nq)
	return z, θ
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

function get_mehrotra_expanded_solver(s::Simulation, z, θ;
		r_tol = 1e-8, κ_tol = 1e-8, linear = false)
	model = s.model
	env = s.env
	lin = LinearizedStep(s, z, θ, 0.0) # /!| here we do not use κ from ref_traj
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
			opts = Mehrotra12Options(
				max_iter_inner = 100,
				r_tol = r_tol,
				κ_tol = κ_tol,
				solver = :empty_solver,
				κ_reg = 1e-3,
				γ_reg = 1.0,
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
			opts = Mehrotra12Options(
		        max_iter_inner = 100,
		        r_tol = r_tol,
		        κ_tol = κ_tol,
				κ_reg = 1e-3,
				γ_reg = 1.0,
				))
	end
	return ip
end

function get_mehrotra_latest_solver(s::Simulation, z, θ;
		r_tol = 1e-8, κ_tol = 1e-8, linear = false)
	model = s.model
	env = s.env
	lin = LinearizedStep(s, z, θ, 0.0) # /!| here we do not use κ from ref_traj
	if linear
		ip = mehrotra_latest(z, θ,
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
			opts = Mehrotra113Options(
				max_iter_inner = 100,
				r_tol = r_tol,
				κ_tol = κ_tol,
				solver = :empty_solver,
				κ_reg = 1e-3,
				γ_reg = 1.0,
				ls_relaxation = 1.0,
				))
	else
		ip = mehrotra_latest(z, θ,
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
			opts = Mehrotra113Options(
		        max_iter_inner = 100,
		        r_tol = r_tol,
		        κ_tol = κ_tol,
				κ_reg = 1e-3,
				γ_reg = 1.0,
				ls_relaxation = 1.0,
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

mutable struct BenchmarkStatistics15{T}
	failure::AbstractVector{Bool}
	iter::AbstractVector{Int}
	r_vio::AbstractVector{T}
	κ_vio::AbstractVector{T}
	cond_schur::AbstractVector{T}
end

function BenchmarkStatistics15()
	failure = trues(0)
	iter = zeros(Int, 0)
	r_vio = zeros(0)
	κ_vio = zeros(0)
	cond_schur = zeros(0)
	return BenchmarkStatistics15(failure, iter, r_vio, κ_vio, cond_schur)
end

function record!(stats::BenchmarkStatistics15, s, i, r, κ, c)
	push!(stats.failure, s)
	push!(stats.iter, i)
	push!(stats.r_vio, r)
	push!(stats.κ_vio, κ)
	push!(stats.cond_schur, c)
	return nothing
end

function benchmark_interior_point!(stats::BenchmarkStatistics15, s::Simulation,
	ref_traj::ContactTraj; opts::BenchmarkOptions12)

	nz = num_var(s.model, s.env)
	for t = 1:ref_traj.H
		z, θ = get_benchmark_initialization(ref_traj, t,
		z_dis = opts.z_dis,
		θ_dis = opts.θ_dis)
		ip = get_interior_point_solver(s, deepcopy(ref_traj.z[t]), deepcopy(ref_traj.θ[t]),
		r_tol = opts.r_tol,
		κ_tol = opts.κ_tol,
		linear = opts.linear)
		interior_point_solve!(ip, z, θ)
		rv = norm(ip.r, Inf)
		fl = !(rv < opts.r_tol)
		record!(stats, fl, ip.iterations, rv, NaN, NaN)
	end
	return nothing
end

function benchmark_mehrotra_expanded!(stats::BenchmarkStatistics15, s::Simulation,
		ref_traj::ContactTraj; opts::BenchmarkOptions12)

	nz = num_var(s.model, s.env)
	for t = 1:ref_traj.H
		z, θ = get_benchmark_initialization(ref_traj, t,
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
		cond_schur = opts.linear ? log(10, cond(ip.rz.D)) : NaN
		record!(stats, fl, ip.iterations, rv, κv, cond_schur)
	end
	return nothing
end

function benchmark_mehrotra_latest!(stats::BenchmarkStatistics15, s::Simulation,
		ref_traj::ContactTraj; opts::BenchmarkOptions12)

	nz = num_var(s.model, s.env)
	for t = 1:ref_traj.H
		z, θ = get_benchmark_initialization(ref_traj, t,
			z_dis = opts.z_dis,
			θ_dis = opts.θ_dis)
		ip = get_mehrotra_latest_solver(s, deepcopy(ref_traj.z[t]), deepcopy(ref_traj.θ[t]),
			r_tol = opts.r_tol,
			κ_tol = opts.κ_tol,
			linear = opts.linear)
		interior_point_solve!(ip, z, θ)
		rv = residual_violation(ip, ip.r)
		κv = bilinear_violation(ip, ip.r)
		fl = !((rv < opts.r_tol) && (κv < opts.κ_tol))
		cond_schur = opts.linear ? log(10, cond(ip.rz.D)) : NaN
		record!(stats, fl, ip.iterations, rv, κv, cond_schur)
	end
	return nothing
end




function benchmark!(stats::BenchmarkStatistics15, s::Simulation,
		ref_traj::ContactTraj; opts::BenchmarkOptions12, solver::Symbol)
	if solver == :interior_point
		benchmark_interior_point!(stats, s, ref_traj; opts = opts)
	elseif solver == :mehrotra_expanded
		benchmark_mehrotra_expanded!(stats, s, ref_traj; opts = opts)
	elseif solver == :mehrotra_latest
		benchmark_mehrotra_latest!(stats, s, ref_traj; opts = opts)
	else
		@warn "Unknown solver"
	end
	return nothing
end

function benchmark(s::AbstractVector{Simulation},
		ref_traj::AbstractVector{<:ContactTraj}, dist::AbstractVector{<:Real};
		opts = BenchmarkOptions12(), solver::Symbol)

	stats = Vector{BenchmarkStatistics15}()
	nsim = length(s)
	ndis = length(dist)
	@assert length(ref_traj) == nsim

	for i = 1:ndis
		sta = BenchmarkStatistics15()
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
dist = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
tol = 1e-8
opts_nonlin = BenchmarkOptions12(r_tol = tol, κ_tol = tol, linear = false)
opts_lin = BenchmarkOptions12(r_tol = tol, κ_tol = tol, linear = true)
stats = Dict{Symbol, AbstractVector{BenchmarkStatistics15}}()
keys = [
	:interior_point_lin,
	:interior_point_nonlin,
	:mehrotra_expanded_lin,
	:mehrotra_expanded_nonlin,
	:mehrotra_latest_lin,
	:mehrotra_latest_nonlin,
	]

stats[:interior_point_nonlin] = benchmark(
	s,
	deepcopy.(ref_traj),
	dist,
	opts = opts_nonlin,
	solver = :interior_point,
	)

stats[:interior_point_lin] = benchmark(
	s,
	deepcopy.(ref_traj),
	dist,
	opts = opts_lin,
	solver = :interior_point,
	)

stats[:mehrotra_expanded_nonlin] = benchmark(
	s,
	deepcopy.(ref_traj),
	dist,
	opts = opts_nonlin,
	solver = :mehrotra_expanded,
	)

stats[:mehrotra_expanded_lin] = benchmark(
	s,
	deepcopy.(ref_traj),
	dist,
	opts = opts_lin,
	solver = :mehrotra_expanded,
	)

stats[:mehrotra_latest_nonlin] = benchmark(
	s,
	deepcopy.(ref_traj),
	dist,
	opts = opts_nonlin,
	solver = :mehrotra_latest,
	)

stats[:mehrotra_latest_lin] = benchmark(
	s,
	deepcopy.(ref_traj),
	dist,
	opts = opts_lin,
	solver = :mehrotra_latest,
	)

function display_statistics!(plt, plt_ind, stats::AbstractVector{BenchmarkStatistics15}, dist;
		field::Symbol = :iter, solver::String = "solver")
	data = [mean(getfield(s, field)) for s in stats]
	linestyle = occursin("nonlin", solver) ? :solid : :dot
	linewidth = occursin("nonlin", solver) ? 3.0 : 5.0
	marker    = occursin("nonlin", solver) ? :circ : :utriangle
	if occursin("interior_point", solver)
		color = :orange
	elseif occursin("mehrotra_expanded", solver)
		color = :cornflowerblue
	elseif occursin("mehrotra_latest", solver)
		color = :red
	end

	plot!(plt[plt_ind...], dist, data,
		xlabel = "disturbances",
		ylabel = String(field),
		label = solver,
		xlims = [minimum(dist), maximum(dist)*100],
		ylims = [0, Inf],
		xaxis = :log,
		legend = :right,
		linestyle = linestyle,
		color = color,
		linewidth = linewidth)
	scatter!(plt[plt_ind...], dist, data,
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

plt = plot(size = (800, 600), layout = (3, 1))
for k in keys
	display_statistics!(plt, [1,1], stats[k], dist, field = :iter, solver = String(k))
end
for k in keys
	display_statistics!(plt, [2,1], stats[k], dist, field = :failure,  solver = String(k))
end
for k in keys
	display_statistics!(plt, [3,1], stats[k], dist, field = :cond_schur,  solver = String(k))
end



# Dev


z2, θ2 = get_benchmark_initialization(ref_traj[1], 10)
ip2 = mehrotra_latest(z2, θ2,
	ix = linearization_var_index(s[1].model, s[1].env)[1],
	iy1 = linearization_var_index(s[1].model, s[1].env)[2],
	iy2 = linearization_var_index(s[1].model, s[1].env)[3],
	idyn = linearization_term_index(s[1].model, s[1].env)[1],
	irst = linearization_term_index(s[1].model, s[1].env)[2],
	ibil = linearization_term_index(s[1].model, s[1].env)[3],
    idx_ineq = inequality_indices(s[1].model, s[1].env),
    idx_soc = soc_indices(s[1].model, s[1].env),
	r! = s[1].res.r!,
    rm! = s[1].res.rm!,
    rz! = s[1].res.rz!,
    rθ! = s[1].res.rθ!,
    rz = s[1].rz,
    rθ = s[1].rθ,
    opts = Mehrotra113Options(
        max_iter_inner=30,
        r_tol=1e-8,
        κ_tol=1e-8,
		# verbose=true
		))
interior_point_solve!(ip2, z2, θ2)
r2 = zeros(num_var(s[1].model, s[1].env))
ip2.methods.r!(r2, ip2.z, ip2.θ, 0.0)
@testset "Mehrotra Nonlinear" begin
	@test norm(r2, Inf) < 1e-8
	@test abs(norm(r2, Inf) - 7.84e-9) < 1e-11
	@test ip2.iterations == 9
end
