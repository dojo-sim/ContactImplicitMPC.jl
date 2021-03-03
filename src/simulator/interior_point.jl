# check that inequality constraints are satisfied
function inequality_check(z, idx_ineq)
    if any(view(z, idx_ineq) .<= 0.0)
        return true
    else
        return false
    end
end

# residual
function r!(r, z, θ, κ)

end

# residual Jacobian wrt z
function rz!(rz, z, θ)
    nothing
end

# residual Jacobian wrt θ
function rθ!(rθ, z, θ)
    nothing
end

# data structure containing pre-allocated memory for interior-point solver
struct InteriorPointData{T}
    z::Vector{T} # current point
    w::Vector{T} # candidate point
    r::Vector{T} # residual
    r_norm::T    # residual norm
    r̄::Vector{T} # candidate residual
    r̄_norm::T    # candidate residual norm
    rz           # residual Jacobian wrt z
    rθ           # residual Jacobian wrt θ
    Δ::Vector{T} # search direction

    θ::Vector{T} # problem data
    idx_ineq     # indices for inequality constraints

    δz           # solution gradients
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

"""
    interior point solver
"""
function interior_point!(data,
    r_tol = 1.0e-5,
    κ_tol = 1.0e-5,
    κ_init = 1.0,
    κ_scale = 0.1,
    max_iter = 100,
    max_ls = 50,
    diff_sol = false)

    # unpack pre-allocated data
    z = data.z
    w = data.w
    θ = data.θ
    r = data.r
    r_norm = data.r_norm
    r̄ = data.r̄
    r̄_norm = data.r̄_norm
    rz = data.rz
    Δ = data.Δ

    # initialize barrier parameter
    κ = κ_init

    # compute residual, residual Jacobian
    r!(r, z, θ, κ)
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
            Δ .= rz \ r

            # initialize step length
            α = 1.0

            # candidate point
            w .= z - α * Δ

            # check inequality constraints
            iter = 0
            while inequality_check(w, idx_ineq)
                α = 0.5 * α
                iter += 1
                if iter > max_ls
                    @error "backtracking line search fail"
                    return false
                end
            end

            # reduce norm of residual
            r!(r̄, w, θ, κ)
            r̄_norm = norm(r̄, Inf)

            while r̄_norm^2.0 >= (1.0 - 0.001 * α) * r_norm^2.0
                α = 0.5 * α
                w .= z - α * Δ
                r!(r̄, w, θ, κ)
                r̄_norm = norm(r̄, Inf)

                iter += 1
                if iter > max_ls
                    @error "line search fail"
                    return false
                end
            end

            # update
            z .= w
            r .= r̄
            r_norm = r̄_norm
        end

        κ < κ_tol ? (return true) : (κ *= κ_scale)
    end

    diff_sol && differentiate_solution!(data)

    return true
end
