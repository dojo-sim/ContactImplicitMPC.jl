# check that inequality constraints are satisfied
function inequality_check(z, idx_ineq)
    if any(view(z, idx_ineq) .<= 0.0)
        return true
    else
        return false
    end
end

# residual
# function r!(r, z, θ)
#     @warn "residual not defined"
#     nothing
# end
#
# # residual Jacobian wrt z
# function rz!(rz, z, θ)
#     @warn "residual Jacobian wrt z not defined"
#     nothing
# end
#
# # residual Jacobian wrt θ
# function rθ!(rθ, z, θ)
#     @warn "residual Jacobian wrt θ not defined"
#     nothing
# end

mutable struct ResidualData{T}
    r::Vector{T} # pre-allocated memory for residual
    θ::Vector{T} # pre-allocated memory for problem data
    κ::T         # barrier parameter
    info::Dict   # additional info
end

function residual_data(num_var, num_data)
    ResidualData(zeros(num_var), zeros(num_data), 0.0, Dict())
end

# data structure containing pre-allocated memory for interior-point solver
struct InteriorPointData{T}
    z::Vector{T}           # current point
    z̄::Vector{T}           # candidate point
    r::Vector{T}           # residual
    r_norm::T              # residual norm
    r̄::Vector{T}           # candidate residual
    r̄_norm::T              # candidate residual norm
    rz::SparseMatrixCSC{T,Int}         # residual Jacobian wrt z
    rθ::SparseMatrixCSC{T,Int}         # residual Jacobian wrt θ
    Δ::Vector{T}           # search direction
    idx_ineq::Vector{Int}  # indices for inequality constraints
    δz::SparseMatrixCSC{T,Int}         # solution gradients
    data::ResidualData{T}  # residual data
end

function interior_point_data(num_var, num_data, idx_ineq;
        rz = spzeros(num_var, num_var),
        rθ = spzeros(num_var, num_data))

    r_data = residual_data(num_var, num_data)

    InteriorPointData(
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
        spzeros(num_var, num_data),
        r_data)
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
function interior_point!(data::InteriorPointData;
    opts = InteriorPointOptions())

    # options
    r_tol = opts.r_tol
    κ_tol = opts.κ_tol
    κ_init = opts.κ_init
    κ_scale = opts.κ_scale
    max_iter = opts.max_iter
    max_ls = opts.max_ls
    diff_sol = opts.diff_sol

    # unpack pre-allocated data
    z = data.z
    z̄ = data.z̄
    θ = data.data
    r = data.r
    r_norm = data.r_norm
    r̄ = data.r̄
    r̄_norm = data.r̄_norm
    rz = data.rz
    Δ = data.Δ
    idx_ineq = data.idx_ineq

    # initialize barrier parameter
    θ.κ = κ_init

    # compute residual, residual Jacobian
    r!(r, z, θ)
    r_norm = norm(r, Inf)

    for k = 1:10
        for i = 1:max_iter
            # check for converged residual
            if r_norm < r_tol
                continue
            end

            # compute residual Jacobian
            rz!(rz, z, θ)

            # compute step
            # Δ .= r
            # info = LAPACK.getrf!(rz)
            # LAPACK.getrs!('N', info[1], info[2], Δ)
            Δ .= rz \ r

            # initialize step length
            α = 1.0

            # candidate point
            z̄ .= z - α * Δ

            # check inequality constraints
            iter = 0
            while inequality_check(z̄, idx_ineq)
                α = 0.5 * α
                z̄ .= z - α * Δ
                iter += 1
                if iter > max_ls
                    @error "backtracking line search fail"
                    return false
                end
            end

            # reduce norm of residual
            r!(r̄, z̄, θ)
            r̄_norm = norm(r̄, Inf)

            while r̄_norm^2.0 >= (1.0 - 0.001 * α) * r_norm^2.0
                α = 0.5 * α
                z̄ .= z - α * Δ
                r!(r̄, z̄, θ)
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

        if θ.κ < κ_tol
            return true
        else
            # update barrier parameter
            θ.κ *= κ_scale

            # update residual
            r!(r, z, θ)
            r_norm = norm(r, Inf)
        end
    end

    diff_sol && differentiate_solution!(data)

    return true
end

function differentiate_solution!(data)
    z = data.z
    θ = data.θ
    rz = data.rz
    rθ = data.rθ
    δz = data.δz

    rz!(rz, z, θ) # maybe not needed
    rθ!(rθ, z, θ)

    δz .= -1.0 * rz \ rθ
    nothing
end
