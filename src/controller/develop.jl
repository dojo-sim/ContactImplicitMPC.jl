
res_path = joinpath(@__DIR__, "../../src/dynamics")
include(joinpath(res_path, "particle/model.jl"))
model = particle
instantiate_base!(model, joinpath(res_path, "particle/base.jld2"))
instantiate_dynamics!(model, joinpath(res_path, "particle/dynamics.jld2"))
instantiate_residual!(model, joinpath(res_path, "particle/residual.jld2"))
@load joinpath(res_path, "particle/sparse_jacobians.jld2") rz_sp rθ_sp
model.spa.rz_sp = rz_sp

# time
h = 0.01
T = 100

## DROP
# initial conditions
q0 = @SVector [0.0, 0.0, 1.0]
q1 = @SVector [0.0, 0.0, 1.0]

# simulator
sim = simulator(model, q0, q1, h, T,
    r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
    rz = rz_sp,
    rθ = rθ_sp,
    ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 1.0e-8),
    sim_opts = SimulatorOptions(warmstart = false))

# simulate
status = simulate!(sim)
@test status
@test all(isapprox.(sim.q[end], 0.0, atol = 1.0e-6))






function easy_Rz(model, bil_terms, bil_vars,  Rz0, z)
	# bil_terms, bil_vars = get_bilinear_indices(model)
	_Rz = deepcopy(Rz0)
	for i = 1:length(bil_terms)
		t = bil_terms[i]
		v1 = bil_vars[i][1]
		v2 = bil_vars[i][2]
		_Rz[t,v1] = Diagonal(z[v2])
		_Rz[t,v2] = Diagonal(z[v1])
	end
	return _Rz
end


function triv_lin_rzθ(model::ContactDynamicsModel, lin::LinStep13, z, θ, κ)
	@assert norm(κ - lin.κ)/κ < 1e-10
	r = lin.r0 + lin.Rz0 * (z-lin.z0) + lin.Rθ0 * (θ-lin.θ0)
	# Bilinearities
	# γ2 s1 - κ
	# ψ s2 - κ
	# b2 η - κ
	for i = 1:length(lin.bil_terms)
		t = lin.bil_terms[i]
		v1 = lin.bil_vars[i][1]
		v2 = lin.bil_vars[i][2]
		# r[t] = z[v1].*z[v2] .- κ
		bil_addition!(r, t, z[v1], z[v2], κ)
	end
	return r
end


function r_approx!(model::ContactDynamicsModel, lin::LinStep14, r::AbstractVector{T},
	z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where {T}
	@assert norm(κ - lin.κ0)/κ < 1e-10
	r .= lin.r0 + lin.rz0 * (z-lin.z0) + lin.rθ0 * (θ-lin.θ0)
	for i = 1:length(lin.bil_terms)
		t = lin.bil_terms[i]
		v1 = lin.bil_vars[i][1]
		v2 = lin.bil_vars[i][2]
		# r[t] = z[v1].*z[v2] .- κ
		bil_addition!(r, t, z[v1], z[v2], κ)
	end
    return nothing
end

function rz_approx!(model::ContactDynamicsModel, lin::LinStep14, rz::AbstractMatrix{T},
	z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where {T}
	rz .= lin.rz0
	for i = 1:length(lin.bil_terms)
		t = lin.bil_terms[i]
		v1 = lin.bil_vars[i][1]
		v2 = lin.bil_vars[i][2]
		rz[t,v1] .= Diagonal(z[v2])
		rz[t,v2] .= Diagonal(z[v1])
	end
    return nothing
end



################################################################################
# Test
################################################################################

model = get_model("quadruped")
model = get_model("particle")

nz = num_var(model)
nθ = num_data(model)
z_ = rand(SizedVector{nz,Float64})
θ_ = rand(SizedVector{nθ,Float64})
κ_ = 1e-5
r0_ = rand(SizedVector{nz,Float64})
r1_ = rand(SizedVector{nz,Float64})
rz0_ = spzeros(nz,nz)
rz0_ = similar(model.spa.rz_sp, Float64)
rz1_ = deepcopy(rz0_)
rθ_ = zeros(nz,nθ)
lin = LinStep14(model, z_, θ_, κ_)

# Test r!
model.res.r(r0_, z_, θ_, κ_)
r_approx!(model, lin, r1_, z_, θ_, κ_)
@test norm(r0_ - r1_, Inf) < 1e-8

α = 1.1
r_approx!(model, lin, r1_, α*z_, α*θ_, κ_)
@test r0_ != r1_

# Test rz!
model.res.rz(rz0_, z_, θ_, κ_)
rz_approx!(model, lin, rz1_, z_, θ_, κ_)
@test norm(rz0_ - rz1_, Inf) < 1e-8

α = 1.1
rz_approx!(model, lin, rz1_, α*z_, α*θ_, κ_)
@test r0_ != r1_
















r1! = model.res.r
rz1! = model.res.rz
rθ1! = model.res.rθ
nz = num_var(model)
nθ = num_data(model)
z_ = rand(SizedVector{nz,Float64})
θ_ = rand(SizedVector{nθ,Float64})
κ_ = 1e-5
r_ = rand(SizedVector{nz,Float64})
rz_ = spzeros(nz,nz)
rz_ = similar(rz_sp, Float64)
rθ_ = zeros(nz,nθ)

r1!(r_, z_, θ_, κ_)
rz1!(rz_, z_, θ_, κ_)
rθ1!(rθ_, z_, θ_, κ_)

LinStep14(model, z_, θ_, κ_)





r_
rz_
rθ_

a = 10
a = 10
a = 10
a = 10


a = [1 0 ; 0 0]
sa = sparse(a)
b = spzeros(2,2)
sb = similar(sa, Float64)
sb



# check that inequality constraints are satisfied
inequality_check(x, idx_ineq) = any(view(x, idx_ineq) .<= 0.0) ? true : false
inequality_check(x) = any(x .<= 0.0) ? true : false


# residual
function r!(r, z, θ, κ)
    @warn "residual not defined"
    nothing
end

# residual Jacobian wrt z
function rz!(rz, z, θ, κ)
    @warn "residual Jacobian wrt z not defined"
    nothing
end

# residual Jacobian wrt θ
function rθ!(rθ, z, θ, κ)
    @warn "residual Jacobian wrt θ not defined"
    nothing
end

struct InteriorPointMethods
    r!
    rz!
    rθ!
end

struct InteriorPoint{T}
    methods::InteriorPointMethods
    z::Vector{T}               # current point
    z̄::Vector{T}               # candidate point
    r::Vector{T}               # residual
    r_norm::T                  # residual norm
    r̄::Vector{T}               # candidate residual
    r̄_norm::T                  # candidate residual norm
    rz#::SparseMatrixCSC{T,Int} # residual Jacobian wrt z
    rθ#::SparseMatrixCSC{T,Int} # residual Jacobian wrt θ
    Δ::Vector{T}               # search direction
    idx_ineq::Vector{Int}      # indices for inequality constraints
    z̄_ineq                     # variables subject to inequality constraints
    δz::SparseMatrixCSC{T,Int} # solution gradients
    θ::Vector{T}               # problem data
    κ::Vector{T}               # barrier parameter
    num_var::Int
    num_data::Int
end

function interior_point(num_var::Int, num_data::Int, idx_ineq::Vector{Int};
        r! = r!, rz! = rz!, rθ! = rθ!,
        rz = spzeros(num_var, num_var),
        rθ = spzeros(num_var, num_data)) where T

    InteriorPoint(
        InteriorPointMethods(r!, rz!, rθ!),
        zeros(num_var),
        zeros(num_var),
        zeros(num_var),
        0.0,
        zeros(num_var),
        0.0,
        rz,
        rθ,
        zeros(num_var),
        idx_ineq,
        view(zeros(num_var), idx_ineq),
        spzeros(num_var, num_data),
        zeros(num_data),
        zeros(1),
        num_var,
        num_data)
end

# interior-point solver options
@with_kw mutable struct InteriorPointOptions{T}
    r_tol::T = 1.0e-5
    κ_tol::T = 1.0e-5
    κ_init::T = 1.0
    κ_scale::T = 0.1
    max_iter::Int = 100
    max_ls::Int = 50
    diff_sol::Bool = false
end

# interior point solver
function interior_point!(ip::InteriorPoint{T};
    opts = InteriorPointOptions{T}()) where T

    # methods
    r! = ip.methods.r!
    rz! = ip.methods.rz!
    rθ! = ip.methods.rθ!

    # options
    r_tol = opts.r_tol
    κ_tol = opts.κ_tol
    κ_init = opts.κ_init
    κ_scale = opts.κ_scale
    max_iter = opts.max_iter
    max_ls = opts.max_ls
    diff_sol = opts.diff_sol

    # unpack pre-allocated data
    z = ip.z
    z̄ = ip.z̄
    r = ip.r
    r_norm = ip.r_norm
    r̄ = ip.r̄
    r̄_norm = ip.r̄_norm
    rz = ip.rz
    Δ = ip.Δ
    idx_ineq = ip.idx_ineq
    ip.z̄_ineq .= view(ip.z̄, ip.idx_ineq)
    θ = ip.θ
    κ = ip.κ

    # initialize barrier parameter
    κ[1] = κ_init

    # compute residual, residual Jacobian
    r!(r, z, θ, κ[1])
    r_norm = norm(r, Inf)

    while true
        for i = 1:max_iter
            # check for converged residual
            if r_norm < r_tol
                continue
            end

            # compute residual Jacobian
            rz!(rz, z, θ, κ[1])

            # compute step
            # Δ .= r
            # info = LAPACK.getrf!(Array(rz))
            # LAPACK.getrs!('N', info[1], info[2], Δ)
            Δ .= rz \ r
            # Δ .= lu(rz) \ r

            # initialize step length
            α = 1.0

            # candidate point
            z̄ .= z - α * Δ

            # check inequality constraints
            iter = 0
            while inequality_check(view(z̄, idx_ineq))
                α = 0.5 * α
                z̄ .= z - α * Δ
                iter += 1
                if iter > max_ls
                    @error "backtracking line search fail"
                    return false
                end
            end

            # reduce norm of residual
            r!(r̄, z̄, θ, κ[1])
            r̄_norm = norm(r̄, Inf)

            while r̄_norm^2.0 >= (1.0 - 0.001 * α) * r_norm^2.0
                α = 0.5 * α
                z̄ .= z - α * Δ
                r!(r̄, z̄, θ, κ[1])
                r̄_norm = norm(r̄, Inf)

                iter += 1
                if iter > max_ls
                    @error "line search fail"
                    return false
                end
            end

            # update
            z .= z̄
            r .= r̄
            r_norm = r̄_norm
        end

        if κ[1] < κ_tol
            # differentiate solution
            diff_sol && differentiate_solution!(ip)
            return true
        else
            # update barrier parameter
            κ[1] *= κ_scale

            # update residual
            r!(r, z, θ, κ[1])
            r_norm = norm(r, Inf)
        end
    end
end

function interior_point!(ip::InteriorPoint{T}, z::Vector{T}, θ::Vector{T};
    opts = InteriorPointOptions{T}()) where T
    ip.z .= z
    ip.θ .= θ
    interior_point!(ip, opts = opts)
end

function differentiate_solution!(ip::InteriorPoint)
    z = ip.z
    θ = ip.θ
    rz = ip.rz
    rθ = ip.rθ
    δz = ip.δz
    κ = ip.κ

    ip.methods.rz!(rz, z, θ, κ[1]) # maybe not needed
    ip.methods.rθ!(rθ, z, θ, κ[1])

    δz .= -1.0 * rz \ Array(rθ) # TODO: fix
    nothing
end
