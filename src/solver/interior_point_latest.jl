abstract type LinearSolver end
abstract type AbstractIPSolver end
abstract type AbstractIPOptions end

function interior_point_options(ip_type::Symbol)
    if ip_type == :interior_point
        return Symbol("InteriorPoint113Options")
    elseif ip_type == :mehrotra
        return Symbol("MehrotraOptions")
    else
        error("Unknown ip_type.")
    end
end

mutable struct ResidualMethods
    r!
    rm!
    rz!
    rθ!
end

# residual
function r!(r, z, θ, κ)
    @warn "residual not defined"
    error()
    nothing
end

# residual mehrotra
function rm!(r, z, Δz, θ, κ)
    @warn "residual mehrotra not defined"
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

function r_update!(r, r̄)
    r .= r̄
end

# optimization spaces
abstract type Space end

# Euclidean
struct Euclidean <: Space
    n::Int
end

function candidate_point!(z̄::Vector{T}, ::Euclidean, z::Vector{T}, Δ::Vector{T}, α::T) where T
    z̄ .= z - α .* Δ
end

function update_point!(z::Vector{T}, ::Space, z̄::Vector{T}) where T
    z .= z̄
end

function mapping!(δz, s::Euclidean, δzs, z) # TODO: make allocation free
    δz .= δzs
end

# interior-point solver options
@with_kw mutable struct InteriorPoint113Options{T} <: AbstractIPOptions
    r_tol::T = 1.0e-5
    κ_tol::T = 1.0e-5
    κ_init::T = 1.0
    κ_scale::T = 0.1
    ls_scale::T = 0.5
    max_iter_inner::Int = 100
    max_iter_outer::Int = 10
    max_ls::Int = 50
    max_time::T = 60.0
    diff_sol::Bool = false
    res_norm::Real = Inf
    reg::Bool = false
    reg_pr_init = 0.0
    reg_du_init = 0.0
    solver::Symbol = :lu_solver
    verbose::Bool = false
    warn::Bool = false
end

# regularize Jacobian / Hessian
function regularize!(v_pr, v_du, reg_pr, reg_du)
    v_pr .+= reg_pr
    v_du .-= reg_du
end

mutable struct InteriorPoint113{T} <: AbstractIPSolver
    s::Space
    methods::ResidualMethods
    z::Vector{T}               # current point
    z̄::Vector{T}               # candidate point
    r#::Vector{T}               # residual
    r_norm::T                  # residual norm
    r̄#::Vector{T}               # candidate residual
    r̄_norm::T                  # candidate residual norm
    rz#::SparseMatrixCSC{T,Int} # residual Jacobian wrt z
    rθ#::SparseMatrixCSC{T,Int} # residual Jacobian wrt θ
    Δ::Vector{T}               # search direction
    idx_ineq::Vector{Int}      # indices for inequality constraints
    idx_soc::Vector{Vector{Int}} # indices for second-order cone constraints

    ix
    iy1
    iy2
    idyn
    irst
    ibil

    idx_pr::Vector{Int}        # indices for primal variables
    idx_du::Vector{Int}        # indices for dual variables
    δz::Matrix{T}              # solution gradients (this is always dense)
    δzs::Matrix{T}             # solution gradients (in optimization space; δz = δzs for Euclidean)
    θ::Vector{T}               # problem data
    κ::Vector{T}               # barrier parameter
    num_var::Int
    num_data::Int
    solver::LinearSolver
    v_pr
    v_du
    reg_pr
    reg_du
    iterations::Int
    opts::InteriorPoint113Options
end

function interior_point_latest(z, θ;
        s = Euclidean(length(z)),
        num_var = length(z),
        num_data = length(θ),
        idx_ineq = collect(1:0),
        idx_soc = Vector{Int}[],
        ix = collect(1:0), # useless
        iy1 = collect(1:0), # useless
        iy2 = collect(1:0), # useless
        idyn = collect(1:0), # useless
        irst = collect(1:0), # useless
        ibil = collect(1:0), # useless
        idx_pr = collect(1:s.n),
        idx_du = collect(1:0),
        r! = r!, rm! = rm!, rz! = rz!, rθ! = rθ!,
        r  = zeros(s.n),
        rz = spzeros(s.n, s.n),
        rθ = spzeros(s.n, num_data),
        reg_pr = [0.0], reg_du = [0.0],
        v_pr = view(rz, CartesianIndex.(idx_pr, idx_pr)),
        v_du = view(rz, CartesianIndex.(idx_du, idx_du)),
        opts::InteriorPoint113Options = InteriorPoint113Options()) where T

    # Indices
    nx = length(ix)
    ny = length(iy1)
    ix = SVector{nx, Int}(ix)
    iy1 = SVector{ny, Int}(iy1)
    iy2 = SVector{ny, Int}(iy2)
    idyn = SVector{nx, Int}(idyn)
    irst = SVector{ny, Int}(irst)
    ibil = SVector{ny, Int}(ibil)

    rz!(rz, z, θ) # compute Jacobian for pre-factorization

    InteriorPoint113(
        s,
        ResidualMethods(r!, rm!, rz!, rθ!),
        z,
        zeros(length(z)),
        r,
        0.0,
        deepcopy(r),
        0.0,
        rz,
        rθ,
        zeros(s.n),
        idx_ineq,
        idx_soc,
        ix,
        iy1,
        iy2,
        idyn,
        irst,
        ibil,
        idx_pr,
        idx_du,
        zeros(length(z), num_data),
        zeros(s.n, num_data),
        θ,
        zeros(1),
        num_var,
        num_data,
        eval(opts.solver)(rz),
        v_pr,
        v_du,
        reg_pr, reg_du,
        0,
        opts)
end

# interior point solver
function interior_point_solve!(ip::InteriorPoint113{T}) where T

    # space
    s = ip.s

    # methods
    r! = ip.methods.r!
    rm! = ip.methods.rm!
    rz! = ip.methods.rz!
    rθ! = ip.methods.rθ!

    # options
    opts = ip.opts
    r_tol = opts.r_tol
    κ_tol = opts.κ_tol
    κ_init = opts.κ_init
    κ_scale = opts.κ_scale
    ls_scale = opts.ls_scale
    max_iter_inner = opts.max_iter_inner
    max_iter_outer = opts.max_iter_outer
    max_time = opts.max_time
    max_ls = opts.max_ls
    diff_sol = opts.diff_sol
    res_norm = opts.res_norm
    reg = opts.reg
    verbose = opts.verbose
    warn = opts.warn

    # unpack pre-allocated data
    z = ip.z
    z̄ = ip.z̄
    r = ip.r
    # r_norm = ip.r_norm
    # r̄ = ip.r̄
    # r̄_norm = ip.r̄_norm
    rz = ip.rz
    Δ = ip.Δ
    idx_ineq = ip.idx_ineq
    idx_soc = ip.idx_soc
    ix = ip.ix
    iy1 = ip.iy1
    iy2 = ip.iy2
    idyn = ip.idyn
    irst = ip.irst
    ibil = ip.ibil

    θ = ip.θ
    κ = ip.κ
    v_pr = ip.v_pr
    v_du = ip.v_du
    reg_pr = ip.reg_pr
    reg_du = ip.reg_du
    solver = ip.solver
    ip.iterations = 0
    comp = false

    # initialize barrier parameter
    κ[1] = κ_init

    # initialize regularization
    reg_pr[1] = opts.reg_pr_init
    reg_du[1] = opts.reg_du_init

    # compute residual, residual Jacobian
    warn && @warn "changed"
    # r!(r, z, θ, κ[1])
    r!(r, z, θ, 0.0)

    # r_norm = norm(r, res_norm)
    warn && @warn "get rid of this"
    κ_vio = general_bilinear_violation(z, idx_ineq, idx_soc, iy1, iy2)
    r_vio = residual_violation(ip, r)


    elapsed_time = 0.0

    for i = 1:max_iter_outer
        elapsed_time >= max_time && break
        for j = 1:max_iter_inner
            elapsed_time >= max_time && break
            elapsed_time += @elapsed begin
                # check for converged residual
                warn && @warn "changed the kickout condition"
                # if r_norm < r_tol
                if max(r_vio, κ_vio) < r_tol
                    break
                end
                ip.iterations += 1

                # compute residual Jacobian
                rz!(rz, z, θ)

                # regularize (fixed, TODO: adaptive)
                reg && regularize!(v_pr, v_du, reg_pr[1], reg_du[1])

                # compute step
                warn && @warn "this shoudln't be useful"
                r!(r, z, θ, 0.0)
                linear_solve!(solver, Δ, rz, r)

                # initialize step length
                # α_ls = 1.0

                # # candidate point
                # candidate_point!(z̄, s, z, Δ, α_ls)
                #
                # # check cones
                # iter = 0
                # while cone_check(z̄, idx_ineq, idx_soc)
                #     α_ls *= ls_scale
                #     candidate_point!(z̄, s, z, Δ, α_ls)
                #
                #     iter += 1
                #     if iter > max_ls
                #         @error "backtracking line search fail"
                #         return false
                #     end
                # end
                α_ineq = ineq_step_length(z, Δ, idx_ineq; τ = 1.0)
                α_soc = soc_step_length(z, Δ, idx_soc; τ = 1.0)
                α = min(α_ineq, α_soc)
                # candidate_point!(z̄, s, z, Δ, α)
                # cautious = (α_ls <= α_ineq) && (α_ls <= α_soc)

                # ################################################################
                # ################################################################
                μ, σ = centering(z, Δ, iy1, iy2, α)
                αaff = α
                # Compute corrector residual
                ip.methods.rm!(r, z, Δ, θ, max(σ*μ, κ_tol/5)) # here we set κ = σ*μ, Δ = Δaff
                # Compute corrector search direction
                linear_solve!(solver, Δ, rz, r)
                α_ineq = ineq_step_length(z, Δ, idx_ineq; τ = 0.99)
                α_soc = soc_step_length(z, Δ, idx_soc; τ = 0.99)
                α = min(α_ineq, α_soc)
                # ################################################################


                # reduce norm of residual
                warn && @warn "changed"
                # r!(r̄, z̄, θ, κ[1])
                #####
                #####
                # r!(r̄, z, θ, κ[1])
                # candidate_point!(z̄, s, z, Δ, α)
                candidate_point!(z, s, z, Δ, α)
                #####
                # r̄_norm = norm(r̄, res_norm)

                # warn && @warn "relaxed line search"
                # # while r̄_norm >= (1.0 - 0.001 * α) * r_norm
                # while r̄_norm >= Inf*(1.0 - 0.001 * α) * r_norm
                #     @show "line search"
                #     α *= ls_scale
                #     candidate_point!(z̄, s, z, Δ, α)
                #
                #     r!(r̄, z̄, θ, κ[1])
                #     r̄_norm = norm(r̄, Inf)
                #
                #     iter += 1
                #     if iter > max_ls
                #         @error "line search fail"
                #         return false
                #     end
                # end

                # update
                # update_point!(z, s, z̄)

                # just check that the residual is low even with κ = 0
                r!(r, z, θ, 0.0)
                r_absolute = norm(r, Inf)

                # r_update!(r, r̄)

                κ_vio = general_bilinear_violation(z, idx_ineq, idx_soc, iy1, iy2)
                r_vio = residual_violation(ip, r)
                # r_norm = r̄_norm
                verbose && println(
                    "out:", i,
                    "   in:", j,
                    "   αaff:", scn(αaff),
                    "   α:", scn(α),
                    "   caut:", cautious,
                    "   μσ:", scn(μ*σ),
                    # "   σ:", scn(σ),
                    "   κ:", scn(κ[1]),
                    "   κ_vio:", scn(κ_vio, digits=0),
                    "   r_vio:", scn(r_vio, digits=0),
                    "   r_abs:", scn(r_absolute, digits=0),
                    # "   res∞:", scn(r_norm),
                    )
            end
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
            # r_norm = norm(r, res_norm)
        end
    end
end

function interior_point_solve!(ip::InteriorPoint113{T}, z::AbstractVector{T}, θ::AbstractVector{T}) where T
    ip.z .= z
    ip.θ .= θ
    interior_point_solve!(ip)
end

function differentiate_solution!(ip::InteriorPoint113)
    s = ip.s
    z = ip.z
    θ = ip.θ
    rz = ip.rz
    rθ = ip.rθ
    δz = ip.δz
    δzs = ip.δzs

    κ = ip.κ

    ip.methods.rz!(rz, z, θ) #TODO: maybe not needed
    ip.methods.rθ!(rθ, z, θ)

    linear_solve!(ip.solver, δzs, rz, rθ)
    @inbounds @views @. δzs .*= -1.0
    mapping!(δz, s, δzs, z)

    nothing
end

linear_solve!(solver::LinearSolver, x::Vector{T}, A::Array{T, 2}, b::Vector{T}) where T = linear_solve!(solver, x, sparse(A), b)





################################################################################
# New methods
################################################################################

function residual_violation(ip::InteriorPoint113{T}, r::AbstractVector{T}) where {T}
    max(norm(r[ip.idyn], Inf), norm(r[ip.irst], Inf))
end

function centering(z::AbstractVector{T}, Δaff::AbstractVector{T},
		iy1::SVector{n,Int}, iy2::SVector{n,Int}, αaff::T) where {n,T}
        # using EQ (7) from CVXOPT, we have that s ∘ z  = μ e
        # for positive orthant: s .* z  .= μ
        # for second order cone: s' * z  = μ; s0 * z1 + z0 * s1 = 0
        # μ only depends on the dot products
        # See Section 5.1.3 in CVXOPT
        # The CVXOPT linear and quadratic cone program solvers
	μaff = (z[iy1] - αaff * Δaff[iy1])' * (z[iy2] - αaff * Δaff[iy2]) / n
	μ = z[iy1]' * z[iy2] / n
	σ = clamp(μaff / μ, 0.0, 1.0)^3
	return μ, σ
end

function general_bilinear_violation(z::AbstractVector{T}, idx_ineq, idx_soc, iy1, iy2) where {T}
    nc = Int(length(idx_soc) / 2)
    nsoc = Int(length(idx_soc) / 2)
    nineq = Int(length(idx_ineq) / 2)

    # Split between primals and duals
    idx_soc_p = idx_soc[1:nsoc]
    idx_soc_d = idx_soc[nsoc+1:2nsoc]
    idx_ineq_1 = intersect(iy1, idx_ineq)
    idx_ineq_2 = intersect(iy2, idx_ineq)

    κ_vio = 0.0
    for i = 1:nc
        soc_vio_s = z[idx_soc_p[i]]' * z[idx_soc_d[i]]
        soc_vio_v = z[idx_soc_p[i][1]] * z[idx_soc_d[i][2:end]] + z[idx_soc_d[i][1]] * z[idx_soc_p[i][2:end]]
        κ_vio = max(κ_vio, abs(soc_vio_s))
        κ_vio = max(κ_vio, maximum(abs.(soc_vio_v)))
    end
    ineq_vio = maximum(abs.(z[idx_ineq_1] .* z[idx_ineq_2]))
    κ_vio = max(κ_vio, ineq_vio)
end

function soc_step_length(λ::AbstractVector{T}, Δ::AbstractVector{T};
        τ::T = 0.99, ϵ::T = 1e-14) where {T}
    # check Section 8.2 CVXOPT
    # The CVXOPT linear and quadratic cone program solvers

    # Adding to slack ϵ to make sure that we never get out of the cone
    λ0 = λ[1] #- ϵ
    λ_λ = λ0^2      - λ[2:end]' * λ[2:end]
    if λ_λ < 0.0
        error("should always be positive")
    end
    λ_λ += ϵ
    λ_Δ = λ0 * Δ[1] - λ[2:end]' * Δ[2:end] + ϵ

    ρs = λ_Δ / λ_λ
    ρv = Δ[2:end] / sqrt(λ_λ)
    ρv -= (λ_Δ / sqrt(λ_λ) + Δ[1]) / (λ0 / sqrt(λ_λ) + 1) * λ[2:end] / λ_λ
    # we make sre that the inverse always exists with ϵ,
    # if norm(ρv) - ρs) is negative (Δ is pushing towards a more positive cone)
        # the computation is ignored and we get the maximum value for α = 1.0
    # else we have α = τ / norm(ρv) - ρs)
    # we add ϵ to the denumerator to ensure strict positivity and avoid 1e-16 errors.
    α = 1.0
    if norm(ρv) - ρs > 0.0
        α = min(α, τ / (norm(ρv) - ρs))
    end
    return α
end

function soc_step_length(z::AbstractVector{T}, Δ::AbstractVector{T},
		idx_soc::AbstractVector{Vector{Int}}; τ::T=0.99) where {T}
        # We need to make this much more efficient (allocation free)
        # by choosing the right structure for idx_soc.
    α = 1.0
    for idx in idx_soc
        # we need -Δ here because we will taking the step x - α Δ
        α = min(α, soc_step_length(z[idx], -Δ[idx], τ = τ))
    end
    return α
end

function ineq_step_length(z::AbstractVector{T}, Δ::AbstractVector{T},
		idx_ineq::AbstractVector{Int}; τ::T=0.9995) where {T}
        # We need to make this much more efficient (allocation free)
        # by choosing the right structure for idx_ineq.

    n = Int(length(idx_ineq) / 2)

    ατ_p = 1.0 # primal
    ατ_d = 1.0 # dual
    for i = 1:n
        ip = idx_ineq[i]
        id = idx_ineq[i + n]
        if Δ[ip] > 0.0
            ατ_p = min(ατ_p, τ * z[ip] / Δ[ip])
        end
        if Δ[id] > 0.0
            ατ_d = min(ατ_d, τ * z[id] / Δ[id])
        end
    end
    α = min(ατ_p, ατ_d)
    return α
end
