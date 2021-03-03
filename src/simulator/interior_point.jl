# check that inequality constraints are satisfied
function inequality_check(z, idx_ineq)
    if any(view(z, idx_ineq) .<= 0.0)
        return true
    else
        return false
    end
end

# residual
function r!(r, z, θ)
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

mutable struct ResidualData
    r::Vector # pre-allocated memory for residual
    θ::Vector # pre-allocated memory for problem data
    κ         # barrier parameter
    info::Dict   # additional info
end

function residual_data(num_var, num_data)
    ResidualData(zeros(num_var), zeros(num_data), 0.0, Dict())
end

# data structure containing pre-allocated memory for interior-point solver
struct InteriorPointData
    z # current point
    z̄ # candidate point
    r # residual
    r_norm  # residual norm
    r̄ # candidate residual
    r̄_norm   # candidate residual norm
    rz           # residual Jacobian wrt z
    rθ           # residual Jacobian wrt θ
    Δ::Vector # search direction

    idx_ineq     # indices for inequality constraints

    δz           # solution gradients

    data::ResidualData # residual data
end

function interior_point_data(num_var, num_data, idx_ineq)
    r_data = residual_data(num_var, num_data)

    InteriorPointData(
        zeros(num_var),
        zeros(num_var),
        zeros(num_var),
        0.0,
        zeros(num_var),
        0.0,
        zeros(num_var, num_var),
        zeros(num_var, num_data),
        zeros(num_var),
        idx_ineq,
        zeros(num_var, num_data),
        r_data)
end

# interior-point solver options
@with_kw mutable struct InteriorPointOptions
    r_tol = 1.0e-5
    κ_tol = 1.0e-5
    κ_init = 1.0
    κ_scale = 0.1
    max_iter = 100
    max_ls = 50
    diff_sol = false
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
            Δ .= rz \ r

            # initialize step length
            α = 1.0

            # candidate point
            z̄ .= z - α * Δ

            # check inequality constraints
            iter = 0
            while inequality_check(z̄, idx_ineq)
                α = 0.5 * α
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
                r!(r̄, w, θ)
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

        θ.κ < κ_tol ? (return true) : (θ.κ *= κ_scale)
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
