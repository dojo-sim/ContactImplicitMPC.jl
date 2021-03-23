#
# res_path = joinpath(@__DIR__, "../../src/dynamics")
# include(joinpath(res_path, "particle/model.jl"))
# model = particle
# instantiate_base!(model, joinpath(res_path, "particle/base.jld2"))
# instantiate_dynamics!(model, joinpath(res_path, "particle/dynamics.jld2"))
# instantiate_residual!(model, joinpath(res_path, "particle/residual.jld2"))
# @load joinpath(res_path, "particle/sparse_jacobians.jld2") rz_sp rθ_sp
# model.spa.rz_sp = rz_sp
#
# # time
# h = 0.01
# T = 100
#
# ## DROP
# # initial conditions
# q0 = @SVector [0.0, 0.0, 1.0]
# q1 = @SVector [0.0, 0.0, 1.0]
#
# # simulator
# sim = simulator(model, q0, q1, h, T,
#     r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
#     rz = rz_sp,
#     rθ = rθ_sp,
#     ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
#     sim_opts = SimulatorOptions(warmstart = false))
#
# # simulate
# status = simulate!(sim)
# @test status
# @test all(isapprox.(sim.traj.q[end], 0.0, atol = 1.0e-6))
#
#
#
#
#
#
# function easy_Rz(model, bil_terms, bil_vars,  Rz0, z)
# 	# bil_terms, bil_vars = get_bilinear_indices(model)
# 	_Rz = deepcopy(Rz0)
# 	for i = 1:length(bil_terms)
# 		t = bil_terms[i]
# 		v1 = bil_vars[i][1]
# 		v2 = bil_vars[i][2]
# 		_Rz[t,v1] = Diagonal(z[v2])
# 		_Rz[t,v2] = Diagonal(z[v1])
# 	end
# 	return _Rz
# end
#
#
# function triv_lin_rzθ(model::ContactDynamicsModel, lin::LinearizedStep13, z, θ, κ)
# 	@assert norm(κ - lin.κ)/κ < 1e-10
# 	r = lin.r0 + lin.Rz0 * (z-lin.z0) + lin.Rθ0 * (θ-lin.θ0)
# 	# Bilinearities
# 	# γ2 s1 - κ
# 	# ψ s2 - κ
# 	# b2 η - κ
# 	for i = 1:length(lin.bil_terms)
# 		t = lin.bil_terms[i]
# 		v1 = lin.bil_vars[i][1]
# 		v2 = lin.bil_vars[i][2]
# 		# r[t] = z[v1].*z[v2] .- κ
# 		bil_addition!(r, t, z[v1], z[v2], κ)
# 	end
# 	return r
# end
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# r1! = model.res.r
# rz1! = model.res.rz
# rθ1! = model.res.rθ
# nz = num_var(model)
# nθ = num_data(model)
# z_ = rand(SizedVector{nz,Float64})
# θ_ = rand(SizedVector{nθ,Float64})
# κ_ = 1e-5
# r_ = rand(SizedVector{nz,Float64})
# rz_ = spzeros(nz,nz)
# rz_ = similar(rz_sp, Float64)
# rθ_ = zeros(nz,nθ)
#
# r1!(r_, z_, θ_, κ_)
# rz1!(rz_, z_, θ_, κ_)
# rθ1!(rθ_, z_, θ_, κ_)
#
# LinearizedStep(model, z_, θ_, κ_)
#
#
#
#
#
# r_
# rz_
# rθ_
#
# a = 10
# a = 10
# a = 10
# a = 10
#
#
# a = [1 0 ; 0 0]
# sa = sparse(a)
# b = spzeros(2,2)
# sb = similar(sa, Float64)
# sb
#
#
#
# # check that inequality constraints are satisfied
# inequality_check(x, idx_ineq) = any(view(x, idx_ineq) .<= 0.0) ? true : false
# inequality_check(x) = any(x .<= 0.0) ? true : false
#
#
# # residual
# function r!(r, z, θ, κ)
#     @warn "residual not defined"
#     nothing
# end
#
# # residual Jacobian wrt z
# function rz!(rz, z, θ, κ)
#     @warn "residual Jacobian wrt z not defined"
#     nothing
# end
#
# # residual Jacobian wrt θ
# function rθ!(rθ, z, θ, κ)
#     @warn "residual Jacobian wrt θ not defined"
#     nothing
# end
#
# struct ResidualMethods
#     r!
#     rz!
#     rθ!
# end
#
# struct InteriorPoint{T}
#     methods::ResidualMethods
#     z::Vector{T}               # current point
#     z̄::Vector{T}               # candidate point
#     r::Vector{T}               # residual
#     r_norm::T                  # residual norm
#     r̄::Vector{T}               # candidate residual
#     r̄_norm::T                  # candidate residual norm
#     rz#::SparseMatrixCSC{T,Int} # residual Jacobian wrt z
#     rθ#::SparseMatrixCSC{T,Int} # residual Jacobian wrt θ
#     Δ::Vector{T}               # search direction
#     idx_ineq::Vector{Int}      # indices for inequality constraints
#     z̄_ineq                     # variables subject to inequality constraints
#     δz::SparseMatrixCSC{T,Int} # solution gradients
#     θ::Vector{T}               # problem data
#     κ::Vector{T}               # barrier parameter
#     num_var::Int
#     num_data::Int
# end
#
# function interior_point(num_var::Int, num_data::Int, idx_ineq::Vector{Int};
#         r! = r!, rz! = rz!, rθ! = rθ!,
#         rz = spzeros(num_var, num_var),
#         rθ = spzeros(num_var, num_data)) where T
#
#     InteriorPoint(
#         ResidualMethods(r!, rz!, rθ!),
#         zeros(num_var),
#         zeros(num_var),
#         zeros(num_var),
#         0.0,
#         zeros(num_var),
#         0.0,
#         rz,
#         rθ,
#         zeros(num_var),
#         idx_ineq,
#         view(zeros(num_var), idx_ineq),
#         spzeros(num_var, num_data),
#         zeros(num_data),
#         zeros(1),
#         num_var,
#         num_data)
# end
#
# # interior-point solver options
# @with_kw mutable struct InteriorPointOptions{T}
#     r_tol::T = 1.0e-5
#     κ_tol::T = 1.0e-5
#     κ_init::T = 1.0
#     κ_scale::T = 0.1
#     max_iter::Int = 100
#     max_ls::Int = 50
#     diff_sol::Bool = false
# end
#
# # interior point solver
# function interior_point!(ip::InteriorPoint{T};
#     opts = InteriorPointOptions{T}()) where T
#
#     # methods
#     r! = ip.methods.r!
#     rz! = ip.methods.rz!
#     rθ! = ip.methods.rθ!
#
#     # options
#     r_tol = opts.r_tol
#     κ_tol = opts.κ_tol
#     κ_init = opts.κ_init
#     κ_scale = opts.κ_scale
#     max_iter = opts.max_iter
#     max_ls = opts.max_ls
#     diff_sol = opts.diff_sol
#
#     # unpack pre-allocated data
#     z = ip.z
#     z̄ = ip.z̄
#     r = ip.r
#     r_norm = ip.r_norm
#     r̄ = ip.r̄
#     r̄_norm = ip.r̄_norm
#     rz = ip.rz
#     Δ = ip.Δ
#     idx_ineq = ip.idx_ineq
#     ip.z̄_ineq .= view(ip.z̄, ip.idx_ineq)
#     θ = ip.θ
#     κ = ip.κ
#
#     # initialize barrier parameter
#     κ[1] = κ_init
#
#     # compute residual, residual Jacobian
#     r!(r, z, θ, κ[1])
#     r_norm = norm(r, Inf)
#
#     while true
#         for i = 1:max_iter
#             # check for converged residual
#             if r_norm < r_tol
#                 continue
#             end
#
#             # compute residual Jacobian
#             rz!(rz, z, θ, κ[1])
#
#             # compute step
#             # Δ .= r
#             # info = LAPACK.getrf!(Array(rz))
#             # LAPACK.getrs!('N', info[1], info[2], Δ)
#             Δ .= rz \ r
#             # Δ .= lu(rz) \ r
#
#             # initialize step length
#             α = 1.0
#
#             # candidate point
#             z̄ .= z - α * Δ
#
#             # check inequality constraints
#             iter = 0
#             while inequality_check(view(z̄, idx_ineq))
#                 α = 0.5 * α
#                 z̄ .= z - α * Δ
#                 iter += 1
#                 if iter > max_ls
#                     @error "backtracking line search fail"
#                     return false
#                 end
#             end
#
#             # reduce norm of residual
#             r!(r̄, z̄, θ, κ[1])
#             r̄_norm = norm(r̄, Inf)
#
#             while r̄_norm^2.0 >= (1.0 - 0.001 * α) * r_norm^2.0
#                 α = 0.5 * α
#                 z̄ .= z - α * Δ
#                 r!(r̄, z̄, θ, κ[1])
#                 r̄_norm = norm(r̄, Inf)
#
#                 iter += 1
#                 if iter > max_ls
#                     @error "line search fail"
#                     return false
#                 end
#             end
#
#             # update
#             z .= z̄
#             r .= r̄
#             r_norm = r̄_norm
#         end
#
#         if κ[1] < κ_tol
#             # differentiate solution
#             diff_sol && differentiate_solution!(ip)
#             return true
#         else
#             # update barrier parameter
#             κ[1] *= κ_scale
#
#             # update residual
#             r!(r, z, θ, κ[1])
#             r_norm = norm(r, Inf)
#         end
#     end
# end
#
# function interior_point!(ip::InteriorPoint{T}, z::Vector{T}, θ::Vector{T};
#     opts = InteriorPointOptions{T}()) where T
#     ip.z .= z
#     ip.θ .= θ
#     interior_point!(ip, opts = opts)
# end
#
# function differentiate_solution!(ip::InteriorPoint)
#     z = ip.z
#     θ = ip.θ
#     rz = ip.rz
#     rθ = ip.rθ
#     δz = ip.δz
#     κ = ip.κ
#
#     ip.methods.rz!(rz, z, θ, κ[1]) # maybe not needed
#     ip.methods.rθ!(rθ, z, θ, κ[1])
#
#     δz .= -1.0 * rz \ Array(rθ) # TODO: fix
#     nothing
# end
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# function easy_lin_step(lin::LinearizedStep, q1_ref, q2_ref,
# 	model::ContactDynamicsModel; u2_ref=zeros(SVector{model.dim.u,Float64}),
# 	ρ0=1e0, outer_iter=5, inner_iter=100, tol=1e-8, step_print::Bool=false, z_init=3e-2)
# 	nq = model.dim.q
# 	nu = model.dim.u
# 	nγ = model.dim.γ
# 	nb = model.dim.b
# 	nz = nq + 4nγ + 4nb
#
# 	function rz(z::AbstractVector)
# 		θ = [q1_ref; q2_ref; u2_ref]
# 		return triv_lin_rzθ(model, lin, z, θ, ρ)
# 	end
#
# 	# Jacobian
# 	function Rz(z::AbstractVector)
# 		# differentiate r
# 		return easy_Rz(model, lin.bil_terms, lin.bil_vars, lin.Rz0, z)
# 	end
#
#     # initialize
# 	z = z_init * (ones(eltype(q1_ref), nz) + ones(eltype(q2_ref), nz) +  ones(eltype(u2_ref), nz))
# 	z[1:nq] = copy(q2_ref)
#
#     ρ = ρ0 # barrier parameter
#     flag = false
#
#     for k = 1:outer_iter
#         for i = 1:inner_iter
#             # compute residual, residual Jacobian
#             res = rz(z)
#             if norm(res) < tol
#                 step_print ? println("   iter ($i) - norm: $(scn(norm(res)))") : nothing
#                 # return z, true
#                 flag = true
#                 # continue
#                 break
#             end
#
#             jac = Rz(z)
#             Δ = jac \ res
#
#             # line search the step direction
#             α = 1.0
#             iter = 0
#             while check_variables(model, z - α * Δ) # backtrack inequalities
#                 α = 0.5 * α
#                 iter += 1
#                 if iter > 50
#                     @error "backtracking line search fail"
#                     flag = false
#                     return z, false
#                 end
#             end
#
#             while norm(rz(z - α * Δ))^2.0 >= (1.0 - 0.001 * α) * norm(res)^2.0
#                 α = 0.5 * α
#                 # println("   α = $α")
#                 iter += 1
#                 if iter > 50
#                     @error "line search fail"
#                     flag = false
#                     return z, false
#                 end
#             end
#
# 			# update
#             z = z - α * Δ
#         end
#
#         ρ = 0.1 * ρ
#         step_print ? println("ρ: $(scn(ρ))") : nothing
#     end
#
# 	_Rz = easy_Rz(model, lin.bil_terms, lin.bil_vars, lin.Rz0, z)
# 	_Rθ = lin.Rθ0
# 	∇ = easy_∇(model, lin.bil_terms, lin.bil_vars, lin.Rz0, lin.Rθ0, z)
#     return z, ∇, flag, _Rz, _Rθ
# end
#
#
#
#
#
#
# function easy_lin_step!(sim, t)
#     # unpack
#     model = sim.model
#     q = sim.traj.q
#     u = sim.u
#     w = sim.w
#     h = sim.h
#     ip = sim.ip
#     z = ip.z
#     θ = ip.θ
#
#     # initialize
#     if sim.sim_opts.warmstart
#         z .+= sim.sim_opts.z_warmstart * rand(ip.num_var)
#         sim.ip_opts.κ_init = sim.sim_opts.κ_warmstart
#     else
#         z_initialize!(z, model, q[t+1])
#     end
#     θ_initialize!(θ, model, q[t], q[t+1], u[t], w[t], h)
#
#     # solve
#     status = interior_point!(ip, opts = sim.ip_opts)
#
#     if status
#         # parse result
#         q2, γ, b, _ = unpack_z(model, z)
#         sim.traj.q[t+2] = copy(q2)
#         sim.γ[t] = γ
#         sim.b[t] = b
#
#         if sim.ip_opts.diff_sol
#             nq = model.dim.q
#             nu = model.dim.u
#             nc = model.dim.c
#             nb = model.dim.b
#
#             sim.dq2dq0[t] = view(ip_data.δz, 1:nq, 1:nq)
#             sim.dq2dq1[t] = view(ip_data.δz, 1:nq, nq .+ (1:nq))
#             sim.dq2du[t] = view(ip_data.δz, 1:nq, 2 * nq .+ (1:nu))
#             sim.dγdq0[t] = view(ip_data.δz, nq .+ (1:nc), 1:nq)
#             sim.dγdq1[t] = view(ip_data.δz, nq .+ (1:nc), nq .+ (1:nq))
#             sim.dγdu[t] = view(ip_data.δz, nq .+ (1:nc), 2 * nq .+ (1:nu))
#             sim.dbdq0[t] = view(ip_data.δz, nq + nc .+ (1:nb), 1:nq)
#             sim.dbdq1[t] = view(ip_data.δz, nq + nc .+ (1:nb), nq .+ (1:nq))
#             sim.dbdu[t] = view(ip_data.δz, nq + nc .+ (1:nb), 2 * nq .+ (1:nu))
#         end
#     end
#     return status
# end
#
#
#
#
#
# @with_kw struct SimulatorOptions{T}
#     warmstart::Bool = true
#     z_warmstart::T = 0.001
#     κ_warmstart::T = 0.001
# end
#
# struct Simulator{S,nq,nu,nc,nb,nw}
#     model
#
#     T::Int
#     h::S
#
#     q::Vector{SArray{Tuple{nq},S,1,nq}}
#     u::Vector{SArray{Tuple{nu},S,1,nu}}
#     γ::Vector{SArray{Tuple{nc},S,1,nc}}
#     b::Vector{SArray{Tuple{nb},S,1,nb}}
#     w::Vector{SArray{Tuple{nw},S,1,nw}}
#
#     dq2dq0::Vector{SizedArray{Tuple{nq,nq},S,2,2}}
#     dq2dq1::Vector{SizedArray{Tuple{nq,nq},S,2,2}}
#     dq2du::Vector{SizedArray{Tuple{nq,nu},S,2,2}}
#     dγdq0::Vector{SizedArray{Tuple{nc,nq},S,2,2}}
#     dγdq1::Vector{SizedArray{Tuple{nc,nq},S,2,2}}
#     dγdu::Vector{SizedArray{Tuple{nc,nu},S,2,2}}
#     dbdq0::Vector{SizedArray{Tuple{nb,nq},S,2,2}}
#     dbdq1::Vector{SizedArray{Tuple{nb,nq},S,2,2}}
#     dbdu::Vector{SizedArray{Tuple{nb,nu},S,2,2}}
#
#     ip::InteriorPoint{S}
#     ip_opts::InteriorPointOptions{S}
#
#     sim_opts::SimulatorOptions{S}
# end
#
# function simulator(model, q0::SVector, q1::SVector, h::S, T::Int;
#     u = [@SVector zeros(model.dim.u) for t = 1:T],
#     w = [@SVector zeros(model.dim.w) for t = 1:T],
#     ip_opts = InteriorPointOptions{S}(),
#     r! = r!, rz! = rz!, rθ! = rθ!,
#     rz = spzeros(num_var(model), num_var(model)),
#     rθ = spzeros(num_var(model), num_data(model)),
#     sim_opts = SimulatorOptions{S}()) where S
#
#     nq = model.dim.q
#     nu = model.dim.u
#     nw = model.dim.w
#     nc = model.dim.c
#     nb = model.dim.b
#
#     q = [q0, q1, [@SVector zeros(nq) for t = 1:T]...]
#     γ = [@SVector zeros(nc) for t = 1:T]
#     b = [@SVector zeros(nb) for t = 1:T]
#
#     dq2dq0 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:T]
#     dq2dq1 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:T]
#     dq2du = [SizedMatrix{nq,nu}(zeros(nq, nu)) for t = 1:T]
#     dγdq0 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:T]
#     dγdq1 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:T]
#     dγdu = [SizedMatrix{nc,nu}(zeros(nc, nu)) for t = 1:T]
#     dbdq0 = [SizedMatrix{nb,nq}(zeros(nb, nq)) for t = 1:T]
#     dbdq1 = [SizedMatrix{nb,nq}(zeros(nb, nq)) for t = 1:T]
#     dbdu = [SizedMatrix{nb,nu}(zeros(nb, nu)) for t = 1:T]
#
#     ip = interior_point(
#         num_var(model),
#         num_data(model),
#         inequality_indices(model),
#         r! = r!, rz! = rz!, rθ! = rθ!,
#         rz = rz,
#         rθ = rθ)
#
#     Simulator(
#         model,
#         T, h,
#         q,
#         u,
#         γ,
#         b,
#         w,
#         dq2dq0,
#         dq2dq1,
#         dq2du,
#         dγdq0,
#         dγdq1,
#         dγdu,
#         dbdq0,
#         dbdq1,
#         dbdu,
#         ip,
#         ip_opts,
#         sim_opts)
# end
#
#
#
#
#
#
# T = Float64
# model = get_model("quadruped")
# model = get_model("particle")
# nq = model.dim.q
# nc = model.dim.c
# nu = model.dim.u
# nw = model.dim.w
# nz = num_var(model)
# nθ = num_data(model)
#
# ip_opts = InteriorPointOptions(κ_init=1e-3, κ_tol=1e-3)
# ip = interior_point(
# 	num_var(model),
# 	num_data(model),
# 	inequality_indices(model),
# 	r! = model.res.r,
# 	rz! = model.res.rz,
# 	rθ! = model.res.rθ,
# 	rz = model.spa.rz_sp,
# 	rθ = model.spa.rθ_sp)
#
# q0 = zeros(SizedVector{nq,T})
# q1 = zeros(SizedVector{nq,T})
# u1 = rand(SizedVector{nu,T})
# w1 = rand(SizedVector{nw,T})
# h = 0.01
# θ_initialize!(ip.θ, model, q0, q1, u1, w1, h)
# z_initialize!(ip.z, model, q1)
#
# @btime status = interior_point!(ip, opts = ip_opts)
#
# z0 = deepcopy(ip.z)
# θ0 = deepcopy(ip.θ)
# κ0 = deepcopy(ip.κ)
# lin = LinearizedStep(model, z0, θ0, κ0[1])
#
# # residual
# function r!(r, z, θ, κ)
#     @warn "linearized"
# 	r_linearized!(lin, r, z, θ, κ)
# end
#
# # residual Jacobian wrt z
# function rz!(rz, z, θ, κ)
# 	@warn "linearized"
# 	rz_linearized!(lin, rz, z, θ, κ)
# end
#
# # residual Jacobian wrt θ
# function rθ!(rθ, z, θ, κ)
# 	@warn "linearized"
# 	rθ_linearized!(lin, rθ, z, θ, κ)
# 	nothing
# end
#
#
# ip_linearized = interior_point(
# 	num_var(model),
# 	num_data(model),
# 	inequality_indices(model),
# 	r! = r!, rz! = rz!, rθ! = rθ!,
# 	rz = model.spa.rz_sp,
# 	rθ = model.spa.rθ_sp)
#
# θ_initialize!(ip_linearized.θ, model, q0, q1, u1, w1, h)
# z_initialize!(ip_linearized.z, model, q1)
#
# status = interior_point!(ip_linearized, opts = ip_opts)
# @test norm(ip_linearized.z - z0, Inf) < 1e-8
# @test norm(ip_linearized.θ - θ0, Inf) < 1e-8
# @test norm(ip_linearized.r, Inf) < 1e-5
#
#
#
#
#
#
# function inner(b; lin=100, pol=200)
#     return lin+pol*b
# end
#
# function outer(a, b; kwargs...)
#     @show kwargs
#     return a + inner(b; kwargs...)
#     return nothing
# end
#
# outer(10,10; lin=100, pol=300)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # vis = Visualizer()
# # open(vis)
#
# T = Float64
# model = get_model("particle")
# H = 10
# h = 0.03
# nq = model.dim.q
# nu = model.dim.u
# nw = model.dim.w
# # nγ = model.dim.c
# # nb = model.dim.b
# # nz = num_var(model)
# # nθ = num_data(model)
#
# q0 = SVector{nq,T}([0.0, 0.0, 0.2])
# q1 = SVector{nq,T}([0.1, 0.1, 0.2])
# u = [rand(SizedVector{nu,T}) for k = 1:H]
# w = [rand(SizedVector{nw,T}) for k = 1:H]
# ip_opts = InteriorPointOptions(κ_init=1e-3, κ_tol=1e-3)
# sim0 = simulator2_base(model, q0, q1, h, H;
#     u = [@SVector zeros(model.dim.u) for t = 1:H],
#     w = [@SVector zeros(model.dim.w) for t = 1:H],
#     ip_opts = ip_opts,
#     sim_opts = SimulatorOptions{T}())
#
# simulate!(sim0; verbose = false)
# visualize!(vis, model, sim0.traj.q)
# ref_traj0 = deepcopy(sim0.traj)
# traj0 = deepcopy(ref_traj0)
#
# impl0 = ImplicitTraj(H, model)
# linearization!(model, ref_traj0, impl0)
# implicit_dynamics!(model, traj0, impl0, κ=κ)
#
# impl0
#
#
#
# ################################################################################
# # Test
# ################################################################################
# T = Float64
# model = get_model("particle")
# H = 10
# h = 0.03
# κ = 1e-3
# nq = model.dim.q
#
# q0 = SVector{nq,T}([0.0, 0.0, 0.2])
# q1 = SVector{nq,T}([0.1, 0.1, 0.2])
# ip_opts = InteriorPointOptions(κ_init=κ, κ_tol=κ*2, r_tol=1e-5)
# sim0 = simulator2_base(model, q0, q1, h, H;
#     u = [@SVector zeros(model.dim.u) for t = 1:H],
#     w = [@SVector zeros(model.dim.w) for t = 1:H],
#     ip_opts = ip_opts,
#     sim_opts = SimulatorOptions{T}())
#
# simulate!(sim0; verbose = false)
# ref_traj0 = deepcopy(sim0.traj)
# traj0 = deepcopy(ref_traj0)
#
# impl0 = ImplicitTraj(H, model)
# linearization!(model, ref_traj0, impl0)
# implicit_dynamics!(model, traj0, impl0, κ=κ)
#
# # Implicit dynamics contraint violation is ≈ 0
# for k = 1:H
# 	@test norm(impl0.d[k], Inf) < 1e-3
# end
#
#
#
# impl0.δz
#
# impl0
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # Trajectory Optimizer
# function easy_lin_trajopt(probsize, model, q1_init, q2_init, q1_ref, q2_ref, lint0;
#     ρ0=1e-3,#############################################
#     β=1e1,#############################################
#     outer_iter=1,#############################################
#     res_tol=1e-8,#############################################
#     solver_outer_iter::Int=3,#############################################
#     solver_inner_iter::Int=10,#############################################
#     utr_ref=[zeros(probsize.nu)  for k = 1:probsize.N-1],######################################################
#     qtr_ref=[zeros(probsize.nq)  for k = 1:probsize.N-1],######################################################
#     γtr_ref=[zeros(probsize.nγ)  for k = 1:probsize.N-1],######################################################
#     btr_ref=[zeros(2probsize.nb) for k = 1:probsize.N-1],######################################################
#     utr_wrm=utr_ref,############################################################################
#     qtr_wrm=qtr_ref,############################################################################
#     γtr_wrm=γtr_ref,############################################################################
#     btr_wrm=btr_ref,############################################################################
#     Qu=fill(Diagonal(1e-1*ones(probsize.nu)), probsize.N-1),##########################
#     Qq=fill(Diagonal(1e+1*ones(probsize.nq)), probsize.N-1),##########################
#     Qγ=fill(Diagonal(1e-7*ones(probsize.nγ)), probsize.N-1),##########################
#     Qb=fill(Diagonal(1e-7*ones(2probsize.nb)), probsize.N-1),#########################
#     u_amp=1e-1,
#     live_plot::Bool=false,
#     z_init=z_init,
#     )
#     N = probsize.N
#     nq = probsize.nq
#     nu = probsize.nu
#     nγ = probsize.nγ
#     nb = probsize.nb
#     nz = nq+4nγ+4nb # residual of the step and variables of the step
#     nθ = 2nq+nu # parameters of the step
#     nλ = nq+nγ+2nb # dynamics lagrange multiplier
#     m = nz+nu+nλ # primal dual vars in a step
#     M = (N-1)*m # primal dual vars
#     nr = 2nq+nu+2nγ+4nb # outer problem resiudal size
#
#     traj = zeros(M)
#     # Initilaization
#     set_var_traj!(N,model, traj, utr_wrm, :u)
#     set_var_traj!(N,model, traj, qtr_wrm, :q)
#     set_var_traj!(N,model, traj, γtr_wrm, :γ)
#     set_var_traj!(N,model, traj, btr_wrm, :b)
#
#     function res(traj)
#         r = zeros(eltype(traj), (N-1)*nr)
#         rtr = [view(r, (k-1)*nr .+ (1:nr)) for k=1:N-1]
#         # we organise the residual as follows [u,q,γ,b,constraint]
#         indu = Vector(1:nu)
#         indq = Vector(nu .+ (1:nq))
#         indγ = Vector(nu+nq .+ (1:nγ))
#         indb = Vector(nu+nq+nγ .+ (1:2nb))
#         indc = Vector(nu+nq+nγ+2nb .+ (1:nλ))
#         # indθ = vcat(indq-) # q1,q2,u2
#         inds = vcat(indq, indγ, indb) # q3,γ2,b2
#
#         utr = get_var_traj(N,model,traj,:u)
#         qtr = get_var_traj(N,model,traj,:q)
#         γtr = get_var_traj(N,model,traj,:γ)
#         btr = get_var_traj(N,model,traj,:b)
#         λdtr = get_var_traj(N,model,traj,:λd)
#         for k = 1:N-1
#
#             # optimality and constraints
#             q1, q2 = get_init_config(N, [q1_init, q2_init], traj, k)
#             zk, ∇k, flag = easy_lin_step(lint0[k],
#                 q1, q2, model, u2_ref=utr[k], ρ0=ρ0, outer_iter=outer_iter, z_init=z_init)
#             ck = zk[1:nλ] - [qtr[k]; γtr[k]; btr[k]]
#
#             rtr[k][indu] .+= Qu[k]*(utr[k] - utr_ref[k])
#             rtr[k][indq] .+= Qq[k]*(qtr[k] - qtr_ref[k])
#             rtr[k][indγ] .+= Qγ[k]*(γtr[k] - γtr_ref[k])
#             rtr[k][indb] .+= Qb[k]*(btr[k] - btr_ref[k])
#             rtr[k][indc] .+= ck
#
#             # Constraint derivative
#             # Minus Identity term #∇qk1, ∇γk, ∇bk
#             rtr[k][nu .+ (1:nλ)] += -λdtr[k]
#             # Implicit function theorem part #∇qk_1, ∇qk, ∇uk
#             off = 0
#             ∇q1 = ∇k[1:nλ, off .+ (1:nq)] # ∈ nλ×nθ
#             off += nq
#             ∇q2 = ∇k[1:nλ, off .+ (1:nq)] # ∈ nλ×nθ
#             off += nq
#             ∇u2 = ∇k[1:nλ, off .+ (1:nu)] # ∈ nλ×nθ
#             off += nu
#
#             k >= 3 ? rtr[k-2][indq] += ∇q1'*λdtr[k] : nothing
#             k >= 2 ? rtr[k-1][indq] += ∇q2'*λdtr[k] : nothing
#             rtr[k][indu] += ∇u2'*λdtr[k]
#
#             # # AL terms
#             # if k >= 3
#             #     rtr[k][inds] += 0.5/ρ0 * -I * ck
#             #     k >= 3 ? rtr[k-2][indq] += 0.5/ρ0 * ∇q1'*ck : nothing
#             #     k >= 2 ? rtr[k-1][indq] += 0.5/ρ0 * ∇q2'*ck : nothing
#             #     rtr[k][indu] += 0.5/ρ0 * ∇u2'*ck
#             # end
#         end
#         return r
#     end
#
#     function jac(traj)
#         jcb = spzeros((N-1)*nr,(N-1)*nr)
#         indu = Vector(1:nu)
#         indq = Vector(nu .+ (1:nq))
#         indγ = Vector(nu+nq .+ (1:nγ))
#         indb = Vector(nu+nq+nγ .+ (1:2nb))
#         indc = Vector(nu+nq+nγ+2nb .+ (1:nλ))
#         # indθ = vcat(indq-) # q1,q2,u2
#         inds = vcat(indq, indγ, indb) # q3,γ2,b2
#         indd = vcat(indq .- 2nr, indq .- nr, indu, indq, indγ, indb) # q1,q2,u2,q3,γ2,b2
#
#         utr = get_var_traj(N,model,traj,:u)
#         qtr = get_var_traj(N,model,traj,:q)
#         γtr = get_var_traj(N,model,traj,:γ)
#         btr = get_var_traj(N,model,traj,:b)
#         λdtr = get_var_traj(N,model,traj,:λd)
#
#         jcbtr = [view(jcb, (k-1)*nr .+ (1:nr), (k-1)*nr .+ (1:nr)) for k=1:N-1]
#         for k = 1:N-1
#             jcbtr[k][indu,indu] .+= Qu[k]
#             jcbtr[k][indq,indq] .+= Qq[k]
#             jcbtr[k][indγ,indγ] .+= Qγ[k]
#             jcbtr[k][indb,indb] .+= Qb[k]
#             jcbtr[k][indc,inds] += -I
#             jcbtr[k][inds,indc] += -I
#
#             # Linearization
#             q1, q2 = get_init_config(N, [q1_init, q2_init], traj, k)
#             zk, ∇k, flag = easy_lin_step(lint0[k],
#                 q1, q2, model, u2_ref=utr[k], ρ0=ρ0, outer_iter=outer_iter, z_init=z_init)
#
#             # Implicit function theorem part #∇qk_1, ∇qk, ∇uk
#             off = 0
#             ∇q1 = ∇k[1:nλ, off .+ (1:nq)] # ∈ nλ×nθ
#             off += nq
#             ∇q2 = ∇k[1:nλ, off .+ (1:nq)] # ∈ nλ×nθ
#             off += nq
#             ∇u2 = ∇k[1:nλ, off .+ (1:nu)] # ∈ nλ×nθ
#             off += nu
#
#             ind_k_2 = (k-3)*nr .+ (1:nr)
#             ind_k_1 = (k-2)*nr .+ (1:nr)
#             ind_k   = (k-1)*nr .+ (1:nr)
#
#             k >= 3 ? jcb_2  = view(jcb, ind_k, ind_k_2) : nothing
#             k >= 3 ? jcb_2T = view(jcb, ind_k_2, ind_k) : nothing
#             k >= 2 ? jcb_1  = view(jcb, ind_k, ind_k_1) : nothing
#             k >= 2 ? jcb_1T = view(jcb, ind_k_1, ind_k) : nothing
#             jcb_0 = view(jcb, ind_k, ind_k)
#
#             k >= 3 ? jcb_2[indc, indq] .+= ∇q1 : nothing
#             k >= 2 ? jcb_1[indc, indq] .+= ∇q2 : nothing
#             jcb_0[indc, indu] .+= ∇u2
#
#             k >= 3 ? jcb_2T[indq, indc] .+= ∇q1' : nothing
#             k >= 2 ? jcb_1T[indq, indc] .+= ∇q2' : nothing
#             jcb_0[indu, indc] .+= ∇u2'
#
#             # Regularization
#             jcbtr[k][indc, indc] += -β*ρ0*I
#
#             # AL terms
#             # ∇c = sparse([∇k[1:nλ, 1:nθ] -I(nλ)])
#             # # @show size(∇c)
#             # # @show size(∇c')
#             # # @show size(∇c'*∇c)
#             # k >= 3 ? (@show size(jcb[(k-1)*nr .+ indd, (k-1)*nr .+ indd] )) : nothing
#             # k >= 3 ? jcb[(k-1)*nr .+ indd, (k-1)*nr .+ indd] += 1/ρ0 * ∇c'*∇c : nothing
#         end
#         return jcb
#     end
#
#     function easy_jac(traj; ϵ=1e-5)
#         M = length(traj)
#         indu = get_var_ind(N,model,traj,:u)
#         indq = get_var_ind(N,model,traj,:q)
#         indγ = get_var_ind(N,model,traj,:γ)
#         indb = get_var_ind(N,model,traj,:b)
#         indλd = get_var_ind(N,model,traj,:λd)
#         masktr = [[indu[k]; indq[k]; indγ[k]; indb[k]; indλd[k]] for k=1:N-1]
#         mask = vcat(masktr...)
#
#         Nr = (N-1)*nr
#         jcb = zeros(Nr, Nr)
#         ep = zeros(M)
#         em = zeros(M)
#         for (i,j) in enumerate(mask)
#             @show i/length(mask)
#             ep[j] += ϵ
#             em[j] -= ϵ
#             jcb[:,i] = (res(traj+ep) - res(traj+em))./(2ϵ)
#             ep *= 0.0
#             em *= 0.0
#         end
#         return jcb
#     end
#
#     function solver(traj)
#         indu = get_var_ind(N,model,traj,:u)
#         indq = get_var_ind(N,model,traj,:q)
#         indγ = get_var_ind(N,model,traj,:γ)
#         indb = get_var_ind(N,model,traj,:b)
#         indλd = get_var_ind(N,model,traj,:λd)
#         masktr = [[indu[k]; indq[k]; indγ[k]; indb[k]; indλd[k]] for k=1:N-1]
#         mask = vcat(masktr...)
#         unmask = setdiff(1:length(traj), mask)
#
#         traj_back = deepcopy(traj)
#         traj_back[mask] .= 0.0
#
#         for i = 1:solver_outer_iter
#             live_plot ? visualize!(vis, model, get_var_traj(N,model,traj,:q), Δt=model.dt) : nothing
#             # sleep(1.0)
#             for l = 1:solver_inner_iter
#                 j = jac(traj)
#                 r = res(traj)
#                 Δ = j \ r
#                 if norm(r,1)/length(r) < res_tol
#                     break
#                 end
#
#
#                 # line search the step direction
#                 α = 1.0
#                 iter = 0
#                 traj_trial = deepcopy(traj)
#                 traj_trial[mask] -= α * Δ
#
#                 # plt = plot()
#                 # plot!(reshape(vcat(get_var_traj(N,model,traj,:u)...), (nu,N-1))', label="u")
#                 # display(plt)
#
#                 while norm(res(traj_trial))^2.0 >= (1.0 - 0.001 * α) * norm(r)^2.0
#                     α = 0.5 * α
#                     # println("   α = $α")
#                     iter += 1
#                     if iter > 6
#                         break
#                     end
#                     # if iter > 50
#                     #     @error "line search fail"
#                     #     flag = false
#                     #     return traj, false
#                     # end
#                     traj_trial = deepcopy(traj)
#                     traj_trial[mask] -= α * Δ
#                 end
#
#                 # update
#                 if iter > 6
#                     β = min(β*1.3, 1e2)
#                 else
#                     β = max(1e1, β/1.3)
#                 end
#                 println("ρ0: ", scn(ρ0, digits=0) ,
#                     "     r: ", scn(norm(r,1)/length(r), digits=0),
#                     "     Δ: ", scn(norm(Δ,1)/length(Δ), digits=0),
#                     "     α: ", -Int(round(log(α))))
#                 traj[mask] -= α * Δ
#             end
#             ρ0 /= 10.0
#         end
#         return traj, flag
#     end
#
#     set_var_traj!(N,model, traj, [zeros(nλ) for k=1:N-1], :λd) ##?????????????????????? why
#     traj, flag = solver(traj)
#     return traj#, res(traj), jac(traj)#, easy_jac(traj)
# end
#
#
# function control!(model::ContactDynamicsModel, impl::ImplicitTraj,
#     ref_traj::ContactTraj, cost::CostFunction, s_opts::NewtonOptions{T}) where {T}
#
#     # we have an initial state q0 q1
#
#
#     return nothing
# end
#
# function jac(traj)
#     jcb = spzeros((N-1)*nr,(N-1)*nr)
#     indu = Vector(1:nu)
#     indq = Vector(nu .+ (1:nq))
#     indγ = Vector(nu+nq .+ (1:nγ))
#     indb = Vector(nu+nq+nγ .+ (1:2nb))
#     indc = Vector(nu+nq+nγ+2nb .+ (1:nλ))
#     # indθ = vcat(indq-) # q1,q2,u2
#     inds = vcat(indq, indγ, indb) # q3,γ2,b2
#     indd = vcat(indq .- 2nr, indq .- nr, indu, indq, indγ, indb) # q1,q2,u2,q3,γ2,b2
#
#     utr = get_var_traj(N,model,traj,:u)
#     qtr = get_var_traj(N,model,traj,:q)
#     γtr = get_var_traj(N,model,traj,:γ)
#     btr = get_var_traj(N,model,traj,:b)
#     λdtr = get_var_traj(N,model,traj,:λd)
#
#     jcbtr = [view(jcb, (k-1)*nr .+ (1:nr), (k-1)*nr .+ (1:nr)) for k=1:N-1]
#     for k = 1:N-1
#         jcbtr[k][indu,indu] .+= Qu[k]
#         jcbtr[k][indq,indq] .+= Qq[k]
#         jcbtr[k][indγ,indγ] .+= Qγ[k]
#         jcbtr[k][indb,indb] .+= Qb[k]
#         jcbtr[k][indc,inds] += -I
#         jcbtr[k][inds,indc] += -I
#
#         # Linearization
#         q1, q2 = get_init_config(N, [q1_init, q2_init], traj, k)
#         zk, ∇k, flag = easy_lin_step(lint0[k],
#             q1, q2, model, u2_ref=utr[k], ρ0=ρ0, outer_iter=outer_iter, z_init=z_init)
#
#         # Implicit function theorem part #∇qk_1, ∇qk, ∇uk
#         off = 0
#         ∇q1 = ∇k[1:nλ, off .+ (1:nq)] # ∈ nλ×nθ
#         off += nq
#         ∇q2 = ∇k[1:nλ, off .+ (1:nq)] # ∈ nλ×nθ
#         off += nq
#         ∇u2 = ∇k[1:nλ, off .+ (1:nu)] # ∈ nλ×nθ
#         off += nu
#
#         ind_k_2 = (k-3)*nr .+ (1:nr)
#         ind_k_1 = (k-2)*nr .+ (1:nr)
#         ind_k   = (k-1)*nr .+ (1:nr)
#
#         k >= 3 ? jcb_2  = view(jcb, ind_k, ind_k_2) : nothing
#         k >= 3 ? jcb_2T = view(jcb, ind_k_2, ind_k) : nothing
#         k >= 2 ? jcb_1  = view(jcb, ind_k, ind_k_1) : nothing
#         k >= 2 ? jcb_1T = view(jcb, ind_k_1, ind_k) : nothing
#         jcb_0 = view(jcb, ind_k, ind_k)
#
#         k >= 3 ? jcb_2[indc, indq] .+= ∇q1 : nothing
#         k >= 2 ? jcb_1[indc, indq] .+= ∇q2 : nothing
#         jcb_0[indc, indu] .+= ∇u2
#
#         k >= 3 ? jcb_2T[indq, indc] .+= ∇q1' : nothing
#         k >= 2 ? jcb_1T[indq, indc] .+= ∇q2' : nothing
#         jcb_0[indu, indc] .+= ∇u2'
#
#         # Regularization
#         jcbtr[k][indc, indc] += -β*ρ0*I
#     end
#     return jcb
# end
#
# function res(traj)
#     r = zeros(eltype(traj), (N-1)*nr)
#     rtr = [view(r, (k-1)*nr .+ (1:nr)) for k=1:N-1]
#     # we organise the residual as follows [u,q,γ,b,constraint]
#     indu = Vector(1:nu)
#     indq = Vector(nu .+ (1:nq))
#     indγ = Vector(nu+nq .+ (1:nγ))
#     indb = Vector(nu+nq+nγ .+ (1:2nb))
#     indc = Vector(nu+nq+nγ+2nb .+ (1:nλ))
#     # indθ = vcat(indq-) # q1,q2,u2
#     inds = vcat(indq, indγ, indb) # q3,γ2,b2
#
#     utr = get_var_traj(N,model,traj,:u)
#     qtr = get_var_traj(N,model,traj,:q)
#     γtr = get_var_traj(N,model,traj,:γ)
#     btr = get_var_traj(N,model,traj,:b)
#     λdtr = get_var_traj(N,model,traj,:λd)
#     for k = 1:N-1
#
#         # optimality and constraints
#         q1, q2 = get_init_config(N, [q1_init, q2_init], traj, k)
#         zk, ∇k, flag = easy_lin_step(lint0[k],
#             q1, q2, model, u2_ref=utr[k], ρ0=ρ0, outer_iter=outer_iter, z_init=z_init)
#         ck = zk[1:nλ] - [qtr[k]; γtr[k]; btr[k]]
#
#         rtr[k][indu] .+= Qu[k]*(utr[k] - utr_ref[k])
#         rtr[k][indq] .+= Qq[k]*(qtr[k] - qtr_ref[k])
#         rtr[k][indγ] .+= Qγ[k]*(γtr[k] - γtr_ref[k])
#         rtr[k][indb] .+= Qb[k]*(btr[k] - btr_ref[k])
#         rtr[k][indc] .+= ck
#
#         # Constraint derivative
#         # Minus Identity term #∇qk1, ∇γk, ∇bk
#         rtr[k][nu .+ (1:nλ)] += -λdtr[k]
#         # Implicit function theorem part #∇qk_1, ∇qk, ∇uk
#         off = 0
#         ∇q1 = ∇k[1:nλ, off .+ (1:nq)] # ∈ nλ×nθ
#         off += nq
#         ∇q2 = ∇k[1:nλ, off .+ (1:nq)] # ∈ nλ×nθ
#         off += nq
#         ∇u2 = ∇k[1:nλ, off .+ (1:nu)] # ∈ nλ×nθ
#         off += nu
#
#         k >= 3 ? rtr[k-2][indq] += ∇q1'*λdtr[k] : nothing
#         k >= 2 ? rtr[k-1][indq] += ∇q2'*λdtr[k] : nothing
#         rtr[k][indu] += ∇u2'*λdtr[k]
#
#         # # AL terms
#         # if k >= 3
#         #     rtr[k][inds] += 0.5/ρ0 * -I * ck
#         #     k >= 3 ? rtr[k-2][indq] += 0.5/ρ0 * ∇q1'*ck : nothing
#         #     k >= 2 ? rtr[k-1][indq] += 0.5/ρ0 * ∇q2'*ck : nothing
#         #     rtr[k][indu] += 0.5/ρ0 * ∇u2'*ck
#         # end
#     end
#     return r
# end



#
#
# nz = num_var(model)
# nθ = num_data(model)
# r0 = zeros(nz)
# # z0 = zeros(nz)
# # θ0 = zeros(nθ)
# r!(r0, z0, θ0, κ) = model.linearized.r(r0, z0, θ0, κ, impl.lin[k].z0, impl.lin[k].θ0, impl.lin[k].r0, impl.lin[k].rz0, impl.lin[k].rθ0)
# k = 1
# model.linearized.r(r0, ref_traj0.z[k], ref_traj0.θ[k], κ, impl0.lin[k].z0, impl0.lin[k].θ0, impl0.lin[k].r0, impl0.lin[k].rz0, impl0.lin[k].rθ0)
# @test norm(r0) < 1e-8
#
#
# q0 = Vector(ref_traj0.q[1])
# q1 = Vector(ref_traj0.q[2])
# u1 = Vector(ref_traj0.u[1])
# w1 = Vector(ref_traj0.w[1])
# γ1 = Vector(ref_traj0.γ[1])
# b1 = Vector(ref_traj0.b[1])
# q2 = Vector(ref_traj0.q[3])
# dynamics(model, h, q0, q1, u1, w1, γ1, b1, q2)


include("newton.jl")



# vis = Visualizer()
# open(vis)

n_opts.r_tol = 1e-6
core1 = Newton(H, h, model)
linearization!(model, ref_traj0, impl0)
@time newton_solve!(model, core1, impl0, cost0, ref_traj0, n_opts, initial_offset=false)
@time newton_solve!(model, core1, impl0, cost0, ref_traj0, n_opts, warm_start=true, initial_offset=true)
















T = Float64
κ = 1e-4
# model = get_model("quadruped")
@load joinpath(pwd(), "src/dynamics/quadruped/gaits/gait1.jld2") z̄ x̄ ū h̄ q u γ b

# time
h = h̄
H = length(u)
# H = 15

nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
impl0 = ImplicitTraj(H, model)

# initial conditions
q0 = SVector{model.dim.q}(q[1])
q1 = SVector{model.dim.q}(q[2])

function z_initialize!(z, model::Quadruped, q1)
	nq = model.dim.q
    z .= 1.0e-1
    z[1:nq] = q1
end

sim0 = simulator(model, q0, q1, h, H,
	p = open_loop_policy([SVector{model.dim.u}(h * u[i]) for i=1:H], h),
    r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
    rz = model.spa.rz_sp,
    rθ = model.spa.rθ_sp,
    ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 2κ, κ_init = κ),
    sim_opts = SimulatorOptions(warmstart = true))

simulate!(sim0; verbose = false)
ref_traj0 = deepcopy(sim0.traj)
traj0 = deepcopy(sim0.traj)

norm(q[1]   - traj0.q[1])/norm(traj0.q[1])
norm(q[2]   - traj0.q[2])/norm(traj0.q[2])
norm(q[3]   - traj0.q[3])/norm(traj0.q[3])
norm(h*u[1] - traj0.u[1])/norm(traj0.u[1])
norm(h*γ[1] - traj0.γ[1])/norm(traj0.γ[1])
norm(h*b[1] - traj0.b[1])/norm(traj0.b[1])
norm(h*b[2] - traj0.b[2])/norm(traj0.b[2])
norm(h*b[3] - traj0.b[3])/norm(traj0.b[3])


nz = num_var(model)
r0 = zeros(nz)
model.res.r(r0, traj0.z[1], traj0.θ[1], κ)
@test norm(r0) < 1e-8

# vis = Visualizer()
# open(vis)
visualize!(vis, model, traj0.q)

linearization!(model, ref_traj0, impl0)
@time implicit_dynamics!(model, ref_traj0, impl0, κ=κ)
@test mean([norm(d) for d in impl0.d]) < 1e-6
mean([norm(d) for d in impl0.d])
plot([norm(d,Inf) for d in impl0.d])







δz_ = [impl0.δq0[1] impl0.δq1[1] impl0.δu1[1]]
@test norm(impl0.δz[1][1:size(δz_)[1], 1:size(δz_)[2]] - δz_) == 0.0


cost0 = CostFunction(H, model.dim,
    Qq=fill(Diagonal(1e-2*SizedVector{nq}([0.02,0.02,1,.15,.15,.15,.15,.15,.15,.15,.15,])), H),
    Qu=fill(Diagonal(3e-2*ones(SizedVector{nu})), H),
    Qγ=fill(Diagonal(1e-6*ones(SizedVector{nc})), H),
    Qb=fill(Diagonal(1e-6*ones(SizedVector{nb})), H),
    )
n_opts = NewtonOptions()

core0 = Newton(H, h, model)
core0.r
traj1 = deepcopy(traj0)
for t = 1:H
    traj1.q[t+2] .+= ones(nq)
    traj1.u[t] .+= ones(nu)
    traj1.w[t] .+= ones(nw)
    traj1.γ[t] .+= ones(nc)
    traj1.b[t] .+= ones(nb)
end
@time residual!(model, core0, core0.r, core0.ν, impl0, cost0, traj0, ref_traj0, n_opts)
norm(core0.r.r) == 0.0
@time residual!(model, core0, core0.r, core0.ν, impl0, cost0, traj1, ref_traj0, n_opts)
norm(core0.r.r) == 0.0
#
# off = 0
# core0.r.r[off .+ (1:nq)]
# off += nq
# core0.r.r[off .+ (1:nu)]
# off += nu
# core0.r.r[off .+ (1:nc)]
# off += nc
# core0.r.r[off .+ (1:nb)]
# off += nb
# core0.r.r[off:end]
@time residual!(model, core0, core0.r, core0.ν, impl0, cost0, ref_traj0, ref_traj0, n_opts)
# @allocated residual!(model, core0, core0.r, impl0, cost0, traj0, ref_traj0, n_opts)
# @code_warntype residual!(model, core0, core0.r, impl0, cost0, traj0, ref_traj0, n_opts)
# @benchmark residual!(model, core0, core0.r, impl0, cost0, traj0, ref_traj0, n_opts)

@time jacobian!(model, core0, core0.j, impl0, cost0, n_opts)
@time jacobian!(model, core0, core0.j, impl0, cost0, n_opts)
# @allocated jacobian!(model, core0, core0.j, impl0, cost0, n_opts)
# @code_warntype jacobian!(model, core0, core0.j, impl0, cost0, n_opts)
# @benchmark jacobian!(model, core0, core0.j, impl0, cost0, n_opts)

# core1 = Newton(H, h, model)
# @time newton_solve!(model, core1, impl0, cost0, ref_traj0, n_opts)
visualize!(vis, model, ref_traj0.q, Δt=5*h)

scn(norm(core0.r.r, 1)/length(core0.r.r))
