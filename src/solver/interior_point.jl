abstract type LinearSolver end

# check that inequality constraints are satisfied
function inequality_check(x, idx_ineq)
    for i in idx_ineq
        if x[i] <= 0.0
            return true
        end
    end
    return false
end

# residual
function r!(r, z, θ, κ)
    @warn "residual not defined"
    nothing
end

# residual Jacobian wrt z
function rz!(rz, z, θ)
    @warn "residual Jacobian wrt z not defined"
    nothing
end

# residual Jacobian wrt θ
function rθ!(rθ, z, θ)
    @warn "residual Jacobian wrt θ not defined"
    nothing
end

# interior-point solver options
@with_kw mutable struct InteriorPointOptions{T}
    r_tol::T = 1.0e-5
    κ_tol::T = 1.0e-5
    κ_init::T = 1.0
    κ_scale::T = 0.1
    max_iter_inner::Int = 100
    max_iter_outer::Int = 10
    max_ls::Int = 50
    max_time::T = 60.0
    diff_sol::Bool = false
    solver::Symbol = :lu_solver
end

mutable struct InteriorPointMethods
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
    δz::SparseMatrixCSC{T,Int} # solution gradients
    θ::Vector{T}               # problem data
    κ::Vector{T}               # barrier parameter
    num_var::Int
    num_data::Int
    solver::LinearSolver
    opts::InteriorPointOptions
end

function interior_point(x, θ;
        num_var = length(x),
        num_data = length(θ),
        idx_ineq = collect(1:0),
        r! = r!, rz! = rz!, rθ! = rθ!,
        rz = spzeros(num_var, num_var),
        rθ = spzeros(num_var, num_data),
        opts = InteriorPointOptions()) where T

    rz!(rz, x, θ) # compute Jacobian for pre-factorization

    InteriorPoint(
        InteriorPointMethods(r!, rz!, rθ!),
        x,
        zeros(num_var),
        zeros(num_var),
        0.0,
        zeros(num_var),
        0.0,
        rz,
        rθ,
        zeros(num_var),
        idx_ineq,
        spzeros(num_var, num_data),
        θ,
        zeros(1),
        num_var,
        num_data,
        eval(opts.solver)(rz),
        opts)
end

# interior point solver
function interior_point!(ip::InteriorPoint{T}) where T

    # methods
    r! = ip.methods.r!
    rz! = ip.methods.rz!
    rθ! = ip.methods.rθ!

    # options
    opts = ip.opts
    r_tol = opts.r_tol
    κ_tol = opts.κ_tol
    κ_init = opts.κ_init
    κ_scale = opts.κ_scale
    max_iter_inner = opts.max_iter_inner
    max_iter_outer = opts.max_iter_outer
    max_time = opts.max_time
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
    θ = ip.θ
    κ = ip.κ
    solver = ip.solver

    # initialize barrier parameter
    κ[1] = κ_init

    # compute residual, residual Jacobian
    r!(r, z, θ, κ[1])
    r_norm = norm(r, Inf)

    for i = 1:max_iter_outer
        for j = 1:max_iter_inner
            # check for converged residual
            if r_norm < r_tol
                break
            end

            # compute residual Jacobian
            rz!(rz, z, θ)

            # compute step
            linear_solve!(solver, Δ, rz, r)

            # initialize step length
            α = 1.0

            # candidate point
            z̄ .= z - α * Δ

            # check inequality constraints
            iter = 0
            while inequality_check(z̄, idx_ineq)
                α *= 0.5
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

            while r̄_norm >= (1.0 - 0.001 * α) * r_norm
                α *= 0.5
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

        if κ[1] <= κ_tol
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

function interior_point!(ip::InteriorPoint{T}, z::AbstractVector{T}, θ::AbstractVector{T}) where T
    ip.z .= z
    ip.θ .= θ
    interior_point!(ip)
end

function differentiate_solution!(ip::InteriorPoint)
    z = ip.z
    θ = ip.θ
    rz = ip.rz
    rθ = ip.rθ
    δz = ip.δz
    κ = ip.κ

    ip.methods.rz!(rz, z, θ) # maybe not needed
    ip.methods.rθ!(rθ, z, θ)

    δz .= -1.0 * rz \ Array(rθ) # TODO: fix
    nothing
end

linear_solve!(solver::LinearSolver, x::Vector{T}, A::Array{T, 2}, b::Vector{T}) where T = linear_solve!(solver, x, sparse(A), b)
