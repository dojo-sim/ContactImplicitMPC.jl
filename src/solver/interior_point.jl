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
@with_kw mutable struct InteriorPointOptions{T} #<: AbstractIPOptions
    r_tol::T = 1.0e-5
    κ_tol::T = 1.0e-5
    ls_scale::T = 0.5
    max_iter::Int = 100
    max_ls::Int = 3
    max_time::T = 1e5
    diff_sol::Bool = false
    reg::Bool = false
    ϵ_min = 0.05 # ∈ [0.005, 0.25]
        # smaller -> faster
        # larger  -> slower, more robust
    κ_reg = 1e-3 # bilinear constraint violation level at which regularization is triggered [1e-3, 1e-4]
    γ_reg = 1e-1 # regularization scaling parameters ∈ [0, 0.1]:
        # 0   -> faster & ill-conditioned
        # 0.1 -> slower & better-conditioned
        # simulation choose γ_reg = 0.1
        # MPC choose γ_reg = 0.0
    solver::Symbol = :lu_solver
    undercut::T = 5.0 # the solver will aim at reaching κ_vio = κ_tol / undercut
        # simulation: Inf
        # MPC: 5.0
    verbose::Bool = false
    warn::Bool = false
end

mutable struct InteriorPoint{T,R,RZ,Rθ}
    s::Space
    idx::IndicesOptimization
    methods::ResidualMethods
    z::Vector{T}               # current point
    r::R                       # residual
    rz::RZ                     # residual Jacobian wrt z
    rθ::Rθ                     # residual Jacobian wrt θ
    Δ::Vector{T}               # search direction
    δz::Matrix{T}              # solution gradients (this is always dense)
    δzs::Matrix{T}             # solution gradients (in optimization space; δz = δzs for Euclidean)
    θ::Vector{T}               # problem data
    solver::LinearSolver
    reg_val::T
    iterations::Int
    opts::InteriorPointOptions{T}
    κ::Vector{T}
end

function interior_point(z, θ;
        s = Euclidean(length(z)),
        idx = IndicesOptimization(),
        r! = r!, rz! = rz!, rθ! = rθ!,
        r  = zeros(idx.nΔ),
        rz = spzeros(idx.nΔ, idx.nΔ),
        rθ = spzeros(idx.nΔ, length(θ)),
        opts::InteriorPointOptions = InteriorPointOptions()) where T

    rz!(rz, z, θ) # compute Jacobian for pre-factorization
    num_data = length(θ)

    InteriorPoint{typeof.([z[1], r, rz, rθ])...}(
        s,
        idx,
        ResidualMethods(r!, rz!, rθ!),
        z,
        r,
        rz,
        rθ,
        zeros(idx.nΔ),
        zeros(idx.nz, num_data),
        zeros(idx.nΔ, num_data),
        θ,
        eval(opts.solver)(rz),
        0.0,
        0,
        opts,
        zeros(1),
        )
end

# interior point solver
function interior_point_solve!(ip::InteriorPoint{T,R,RZ,Rθ}) where {T,R,RZ,Rθ}

    # space
    s = ip.s

    # options
    opts = ip.opts
    r_tol = opts.r_tol
    κ_tol = opts.κ_tol
    ls_scale = opts.ls_scale
    max_iter = opts.max_iter
    max_time = opts.max_time
    max_ls = opts.max_ls
    diff_sol = opts.diff_sol
    ϵ_min = opts.ϵ_min
    κ_reg = opts.κ_reg
    γ_reg = opts.γ_reg
    reg = opts.reg
    verbose = opts.verbose
    warn = opts.warn

    # unpack pre-allocated data
    z = ip.z
    r = ip.r
    rz = ip.rz
    Δ = ip.Δ
    θ = ip.θ

    # indices
    idx = ip.idx
    ortz = idx.ortz
    ortΔ = idx.ortΔ
    socz = idx.socz
    socΔ = idx.socΔ
    bil = idx.bil
    ortr = idx.ortr
    socr = idx.socr

    # initialization
    solver = ip.solver
    ip.iterations = 0
    no_progress = 0

    # compute residual, residual Jacobian
    # ip.methods.r!(r, z, θ, 0.0)
    # least_squares!(ip, z, θ, r, rz) # seems to be harmful for performance (failure and iteration count)
    z .= initial_state!(z, ortz, socz) # decrease failure rate for linearized case

    ip.methods.r!(r, z, θ, 0.0)

    κ_vio = bilinear_violation(ip, r)
    r_vio = residual_violation(ip, r)
    elapsed_time = 0.0

    for j = 1:max_iter
        elapsed_time >= max_time && break
        elapsed_time += @elapsed begin
            # check for converged residual
            if (r_vio < r_tol) && (κ_vio < κ_tol)
                break
            end
            ip.iterations += 1

            # Compute regularization level
            # κ_vio = bilinear_violation(ip, r)
            ip.reg_val = κ_vio < κ_reg ? κ_vio * γ_reg : 0.0

            # compute residual Jacobian
            rz!(ip, rz, z, θ, reg = ip.reg_val) # this is not adapted to the second order cone

            # compute step
            linear_solve!(solver, Δ, rz, r, reg = ip.reg_val)

            α_ort = ort_step_length(z, Δ, ortz, ortΔ; τ = 1.0)
            α_soc = soc_step_length(z, Δ, socz, socΔ; τ = 1.0, verbose = false)
            α = min(α_ort, α_soc)
            μ, σ = general_centering(z, Δ, ortz, ortΔ, socz, socΔ, α)
            αaff = α

            # Compute corrector residual
            ip.methods.r!(r, z, θ, max(σ * μ, κ_tol/opts.undercut)) # here we set κ = σ*μ, Δ = Δaff
            general_correction_term!(r, Δ, ortr, socr, ortΔ, socΔ)

            # Compute corrector search direction
            linear_solve!(solver, Δ, rz, r, reg = ip.reg_val, fact = false)
            τ = max(0.95, 1 - max(r_vio, κ_vio)^2)

            α_ort = ort_step_length(z, Δ, ortz, ortΔ; τ = τ)
            α_soc = soc_step_length(z, Δ, socz, socΔ; τ = min(τ, 0.99), verbose = false)
            α = min(α_ort, α_soc)

            # reduce norm of residual
            candidate_point!(z, s, z, Δ, α)
            κ_vio_cand = 0.0
            r_vio_cand = 0.0
            for i = 1:max_ls
                ip.methods.r!(r, z, θ, 0.0)
                κ_vio_cand = bilinear_violation(ip, r)
                r_vio_cand = residual_violation(ip, r)
                if (r_vio_cand <= r_vio) || (κ_vio_cand <= κ_vio)
                    break
                end
                verbose && println("linesearch $i")
                # backtracking
                candidate_point!(z, s, z, Δ, -α * ls_scale^i)
            end

            κ_vio = κ_vio_cand
            r_vio = r_vio_cand
            # nc = Int(idx.ny/4)
            # nb = Int(idx.ny/2)
            verbose && println("iter:", j,
                "  r: ", scn(norm(r, Inf)),
                "  r_vio: ", scn(r_vio),
                "  κ_vio: ", scn(κ_vio),
                # "  κ: ", scn(norm(r[idx.bil[1:nc]], Inf)),
                # "  κ: ", scn(norm(r[idx.bil[nc .+ (1:nb)]], Inf)),
                # "  κ: ", scn(norm(r[idx.bil[nc + nb .+ (1:nc)]], Inf)),
                "  Δ: ", scn(norm(Δ)),
                "  α: ", scn(norm(α)))

            # verbose && println(
            #     "in:", j,
            #     "   αaff:", scn(αaff, digits = 0),
            #     "   α:", scn(α, digits = 0),
            #     "   μσ:", scn(μ*σ, digits = 0),
            #     "   κ_vio:", scn(κ_vio, digits = 0),
            #     "   r_vio:", scn(r_vio, digits = 0),
            #     )
        end
    end
    if (r_vio < r_tol) && (κ_vio < κ_tol)
        # differentiate solution
        ########################################################################
        if false && R <: RLin
        # if R <: RLin
            # regularize solution so that k_vio == k_tol for all bilinear constraints
            ip.methods.r!(r, z, θ, 0.0)
            for i in eachindex(idx.ortz[1])
                if z[idx.ortz[1][i]] <= z[idx.ortz[2][i]]
                    z[idx.ortz[1][i]] *= κ_tol / r.rbil[i]
                else
                    z[idx.ortz[2][i]] *= κ_tol / r.rbil[i]
                end
            end
        else
            # # regularize solution so that k_vio == k_tol for all bilinear constraints
            # ip.methods.r!(r, z, θ, 0.0)
            # for i in eachindex(idx.ortz[1])
            #     if z[idx.ortz[1][i]] <= z[idx.ortz[2][i]]
            #         z[idx.ortz[1][i]] *= κ_tol / r[idx.ortr][i]
            #     else
            #         z[idx.ortz[2][i]] *= κ_tol / r[idx.ortr][i]
            #     end
            # end
            # ip.methods.r!(r, z, θ, 0.0)
            # @show scn.(r[idx.ortr])
        end
        ########################################################################
        diff_sol && differentiate_solution!(ip, reg = max(ip.reg_val, κ_tol * γ_reg))
        return true
    else
        return false
    end
end

function rz!(ip::InteriorPoint, rz::AbstractMatrix{T}, z::AbstractVector{T},
        θ::AbstractVector{T}; reg = 0.0) where {T}
    z_reg = deepcopy(z)
    ortz = ip.idx.ortz
    socz = ip.idx.socz
    for i in eachindex(ortz) # primal-dual
        z_reg[ortz[i]] = max.(z[ortz[i]], reg)
    end
    ip.methods.rz!(rz, z_reg, θ)
    return nothing
end

function rθ!(ip::InteriorPoint, rθ::AbstractMatrix{T}, z::AbstractVector{T},
        θ::AbstractVector{T}) where {T}
    ip.methods.rθ!(rθ, z, θ)
    return nothing
end

function general_correction_term!(r::AbstractVector{T}, Δ::AbstractVector{T},
    ortr::Vector{Int}, socr::Vector{Int},
    ortΔ::Vector{Vector{Int}}, socΔ::Vector{Vector{Vector{Int}}}) where {T}
    # @warn "define residual order"
    r[ortr] .+= Δ[ortΔ[1]] .* Δ[ortΔ[2]]
    r[socr] .+= vcat(
        [second_order_cone_product(
            Δ[socΔ[2][i]],
            Δ[socΔ[1][i]],
        ) for i in eachindex(socΔ[1])]...)
    return nothing
end

function least_squares!(ip::InteriorPoint, z::AbstractVector{T}, θ::AbstractVector{T},
        r::AbstractVector{T}, rz::AbstractMatrix{T}) where {T}
    # doing nothing gives the best result if z_t is correctly initialized with z_t-1 in th simulator
        dyn = ip.idx.dyn
        rst = ip.idx.rst
        A = rz[[dyn; rst], :]
        z .+= A' * ((A * A') \ r[[dyn; rst]])
    return nothing
end

function initial_state!(z::AbstractVector{T}, ortz::Vector{Vector{Int}},
        socz::Vector{Vector{Vector{Int}}}; ϵ::T=1e-20) where {T}

    # Split between primals and duals
    socz_p = socz[1]
    socz_d = socz[2]
    ortz_p = ortz[1]
    ortz_d = ortz[2]

    # ineq
    y1 = z[ortz_p]
    y2 = z[ortz_d]
    δy1 = max(-1.5 * minimum(y1), 0)
    δy2 = max(-1.5 * minimum(y2), 0)

    y1h = y1 .+ δy1
    y2h = y2 .+ δy2
    δhy1 = 0.5 * y1h'*y2h / (sum(y2h) + ϵ)
    δhy2 = 0.5 * y1h'*y2h / (sum(y1h) + ϵ)

    y10 = y1h .+ δhy1
    y20 = y2h .+ δhy2
    z[ortz_p] .= y10
    z[ortz_d] .= y20

    # soc
    for i in eachindex(socz_p)
        e = [1; zeros(length(socz_p[i]) - 1)] # identity element
        y1 = z[socz_p[i]]
        y2 = z[socz_d[i]]
        δy1 = max(-1.5 * (y1[1] - norm(y1[2:end])), 0)
        δy2 = max(-1.5 * (y2[1] - norm(y2[2:end])), 0)

        y1h = y1 + δy1 * e
        y2h = y2 + δy2 * e
        δhy1 = 0.5 * y1h'*y2h / ((y2h[1] + norm(y2h[2,end])) + ϵ)
        δhy2 = 0.5 * y1h'*y2h / ((y1h[1] + norm(y1h[2,end])) + ϵ)

        y10 = y1h + δhy1 * e
        y20 = y2h + δhy2 * e
        z[socz_p[i]] .= y10
        z[socz_d[i]] .= y20
    end
    return z
end

function interior_point_solve!(ip::InteriorPoint, z::AbstractVector{T}, θ::AbstractVector{T}) where T
    ip.z .= z
    ip.θ .= θ
    interior_point_solve!(ip)
end

function differentiate_solution!(ip::InteriorPoint; reg = 0.0)
    s = ip.s
    z = ip.z
    θ = ip.θ
    rz = ip.rz
    rθ = ip.rθ
    δz = ip.δz
    δzs = ip.δzs

    rz!(ip, rz, z, θ, reg = reg)
    rθ!(ip, rθ, z, θ)

    linear_solve!(ip.solver, δzs, rz, rθ, reg = reg)
    @inbounds @views @. δzs .*= -1.0
    mapping!(δz, s, δzs, z)
    return nothing
end


function residual_violation(ip::InteriorPoint, r::AbstractVector{T}) where {T}
    dyn = ip.idx.dyn
    rst = ip.idx.rst
    max(norm(r[dyn], Inf), norm(r[rst], Inf))
end

function general_centering(z::AbstractVector{T}, Δaff::AbstractVector{T},
        ortz::Vector{Vector{Int}}, ortΔ::Vector{Vector{Int}},
        socz::Vector{Vector{Vector{Int}}}, socΔ::Vector{Vector{Vector{Int}}}, αaff::T) where {T}
        # See Section 5.1.3 in CVXOPT
        # μ only depends on the dot products (no cone product)
        # The CVXOPT linear and quadratic cone program solvers

    n = length(ortz[1]) + sum(length.(socΔ[1]))
    # ineq
    μ = z[ortz[1]]' * z[ortz[2]]
    μaff = (z[ortz[1]] - αaff * Δaff[ortΔ[1]])' * (z[ortz[2]] - αaff * Δaff[ortΔ[2]])
    # soc
    for i in eachindex(socz[1])
        μ += z[socz[1][i]]' * z[socz[2][i]]
        μaff += (z[socz[1][i]] - αaff * Δaff[socΔ[1][i]])' * (z[socz[2][i]] - αaff * Δaff[socΔ[2][i]])
    end
    μ /= n
    μaff /= n
	σ = clamp(μaff / μ, 0.0, 1.0)^3
	return μ, σ
end

function bilinear_violation(ip::InteriorPoint, r::AbstractVector{T}) where {T}
    bil = ip.idx.bil
    return norm(r[bil], Inf)
end

function soc_value(u::AbstractVector{T}) where {T}
    u0 = u[1]
    u1 = u[2:end]
    return (u0^2 - u1' * u1)
end

function soc_step_length(λ::AbstractVector{T}, Δ::AbstractVector{T};
        τ::T = 0.99, ϵ::T = 1e-14, verbose::Bool = false) where {T}
    # check Section 8.2 CVXOPT
    # The CVXOPT linear and quadratic cone program solvers

    # Adding to slack ϵ to make sure that we never get out of the cone
    λ0 = λ[1] #- ϵ
    λ_λ = max(λ0^2 - λ[2:end]' * λ[2:end], 1e-25)
    verbose && println(
        "    vλ:", scn(soc_value(λ), digits = 0, exp_digits = 2),
        "    vλ+Δ:", scn(soc_value(λ+Δ), digits = 0, exp_digits = 2),
        "    λ_λ: ", scn(λ_λ, digits = 0, exp_digits = 2),
        "    λ:", scn.(λ, digits = 0, exp_digits = 2),
        "    Δ:", scn.(Δ, digits = 0, exp_digits = 2),
        )
    if λ_λ < 0.0
        @show λ_λ
        @warn "should always be positive"
        # error("should always be positive")
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
    verbose && println(
        "     α:", scn(α, digits = 0, exp_digits = 2))
    verbose && cone_plot(λ, Δ, show_α = true)
    return α
end

function soc_step_length(z::AbstractVector{T}, Δ::AbstractVector{T},
		socz::Vector{Vector{Vector{Int}}}, socΔ::Vector{Vector{Vector{Int}}};
        τ::T=0.99, verbose::Bool = false) where {T}
        # We need to make this much more efficient (allocation free)
    α = 1.0
    for i in eachindex(socz) # primal-dual
        for j in eachindex(socz[i]) # number of cones
            # we need -Δ here because we will taking the step x - α Δ
            α = min(α, soc_step_length(z[socz[i][j]], -Δ[socΔ[i][j]], τ = τ, verbose = verbose))
        end
    end
    return α
end

function ort_step_length(z::AbstractVector{T}, Δ::AbstractVector{T},
		ortz::Vector{Vector{Int}}, ortΔ::Vector{Vector{Int}};
        τ::T=0.9995) where {T}
        # We need to make this much more efficient (allocation free)
    α = 1.0
    for i in eachindex(ortz) # primal-dual
        for j in eachindex(ortz[i])
            k = ortz[i][j] # z
            ks = ortΔ[i][j] # Δz
            if Δ[ks] > 0.0
                α = min(α, τ * z[k] / Δ[ks])
            end
        end
    end
    return α
end
