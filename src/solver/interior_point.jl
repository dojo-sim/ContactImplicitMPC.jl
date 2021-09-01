abstract type LinearSolver end
abstract type AbstractIPSolver end
abstract type AbstractIPOptions end

function interior_point_options(ip_type::Symbol)
    if ip_type == :interior_point
        return Symbol("InteriorPointOptions")
    elseif ip_type == :interior_point_base
        return Symbol("InteriorPointBaseOptions")
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
@with_kw mutable struct InteriorPointOptions{T} <: AbstractIPOptions
    r_tol::T = 1.0e-5
    κ_tol::T = 1.0e-5
    ls_scale::T = 0.5
    max_iter_inner::Int = 100
    # max_ls::Int = 3
    max_ls::Int = 1
    max_time::T = 60.0
    diff_sol::Bool = false
    reg::Bool = false
    reg_pr_init = 0.0
    reg_du_init = 0.0
    ϵ_min = 0.05 # ∈ [0.005, 0.25]
        # smaller -> faster
        # larger  -> slower, more robust
    κ_reg = 1e-3 # bilinear constraint violation level at which regularization is triggered [1e-3, 1e-4]
    γ_reg = 1e-1 # regularization scaling parameters ∈ [0, 0.1]:
        # 0   -> faster & ill-conditioned
        # 0.1 -> slower & better-conditioned
    solver::Symbol = :lu_solver
    verbose::Bool = false
    warn::Bool = false
end

# regularize Jacobian / Hessian
function regularize!(v_pr, v_du, reg_pr, reg_du)
    v_pr .+= reg_pr
    v_du .-= reg_du
end

mutable struct InteriorPoint{T} <: AbstractIPSolver
    s::Space
    methods::ResidualMethods
    z::Vector{T}                # current point
    r#::Vector{T}               # residual
    r̄#::Vector{T}               # candidate residual   #useless
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
    num_var::Int
    num_data::Int
    solver::LinearSolver
    v_pr
    v_du
    reg_pr
    reg_du
    reg_val
    iterations::Int
    opts::InteriorPointOptions

    κ::Vector{T} # useless
end

function interior_point(z, θ;
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
        opts::InteriorPointOptions = InteriorPointOptions()) where T

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

    InteriorPoint(
        s,
        ResidualMethods(r!, rm!, rz!, rθ!),
        z,
        r,
        deepcopy(r),
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
        num_var,
        num_data,
        eval(opts.solver)(rz),
        v_pr,
        v_du,
        reg_pr, reg_du,
        0.0,
        0,
        opts,
        [opts.κ_tol],
        )
end

# interior point solver
function interior_point_solve!(ip::InteriorPoint{T}) where T

    # space
    s = ip.s
    nquat = ip.num_var - s.n

    # options
    opts = ip.opts
    r_tol = opts.r_tol
    κ_tol = opts.κ_tol
    ls_scale = opts.ls_scale
    max_iter_inner = opts.max_iter_inner
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
    idx_ineq = ip.idx_ineq
    idx_soc = ip.idx_soc
    ix = ip.ix
    iy1 = ip.iy1
    iy2 = ip.iy2
    idyn = ip.idyn
    irst = ip.irst
    ibil = ip.ibil

    θ = ip.θ
    v_pr = ip.v_pr
    v_du = ip.v_du
    reg_pr = ip.reg_pr
    reg_du = ip.reg_du
    solver = ip.solver
    ip.iterations = 0
    comp = false

    # initialize regularization
    reg_pr[1] = opts.reg_pr_init
    reg_du[1] = opts.reg_du_init


    # compute residual, residual Jacobian
    ip.methods.r!(r, z, θ, 0.0)
    # TODO need least squares init /!|
    least_squares!(ip, z, θ, r, rz) # seems to be harmful for performance (failure and iteration count)
    z .= initial_state!(z, iy1, iy2, idx_ineq, idx_soc) # decrease failure rate for linearized case
    # println("z1: ", scn(norm(z), digits = 7))

    ip.methods.r!(r, z, θ, 0.0)

    κ_vio = bilinear_violation(ip, r, nquat = nquat)
    r_vio = residual_violation(ip, r)
    elapsed_time = 0.0

    for j = 1:max_iter_inner
        elapsed_time >= max_time && break
        elapsed_time += @elapsed begin
            # check for converged residual
            if (r_vio < r_tol) && (κ_vio < κ_tol)
                break
            end
            ip.iterations += 1

            # Compute regularization level
            κ_vio = bilinear_violation(ip, r, nquat = nquat)
            ip.reg_val = κ_vio < κ_reg ? κ_vio * γ_reg : 0.0

            # compute residual Jacobian
            warn && @warn "changed"
            # ip.methods.rz!(rz, z, θ)
            rz!(ip, rz, z, θ, reg = ip.reg_val) # this is not adapted to the second order cone

            # regularize (fixed, TODO: adaptive)
            reg && regularize!(v_pr, v_du, reg_pr[1], reg_du[1])

            # compute step
            # TODO need to use reg_val here
            warn && @warn "changed"
            # linear_solve!(solver, Δ, rz, r)
            linear_solve!(solver, Δ, rz, r, reg = ip.reg_val)
            # println("Δaff1: ", scn(norm(Δ), digits = 7))

            α_ineq = ineq_step_length(z, Δ, idx_ineq; τ = 1.0, nquat = nquat)
            α_soc = soc_step_length(z, Δ, idx_soc; τ = 1.0, nquat = nquat, verbose = false)
            α = min(α_ineq, α_soc)
            # println("αaff1: ", scn(norm(α), digits = 7))

            μ, σ = general_centering(z, Δ, idx_ineq, idx_soc, iy1, iy2, α, nquat = nquat)
            # println("σ*μ: ", scn(norm(σ*μ), digits = 7))
            αaff = α

            # Compute corrector residual
            # warn && @warn "changed"
            ip.methods.r!(r, z, θ, max(σ * μ, κ_tol/50)) # here we set κ = σ*μ, Δ = Δaff
            # ip.methods.r!(r, z, θ, max(σ * μ, κ_vio/100 , κ_tol/5)) # here we set κ = σ*μ, Δ = Δaff

            # println("r: ", scn(norm(r.rdyn), digits = 7))
            # println("r: ", scn(norm(r.rrst), digits = 7))
            # println("r: ", scn(norm(r.rbil), digits = 7))
            general_correction_term!(r, Δ, ibil, idx_ineq, idx_soc, iy1, iy2, nquat = nquat)
            # println("r: ", scn(norm(r.rbil), digits = 7))

            # Compute corrector search direction
            warn && @warn "changed"
            # linear_solve!(solver, Δ, rz, r)
            linear_solve!(solver, Δ, rz, r, reg = ip.reg_val, fact = false)
            # println("Δ: ", scn(norm(Δ), digits = 7))
            τ = max(0.95, 1 - max(r_vio, κ_vio)^2)
            # println("τ: ", scn(norm(τ), digits = 7))

            α_ineq = ineq_step_length(z, Δ, idx_ineq; τ = τ, nquat = nquat)
            α_soc = soc_step_length(z, Δ, idx_soc; τ = min(τ, 0.99), nquat = nquat, verbose = false)
            α = min(α_ineq, α_soc)
            # println("α: ", scn(norm(α), digits = 7))

            # reduce norm of residual
            candidate_point!(z, s, z, Δ, α)
            κ_vio_cand = 0.0
            r_vio_cand = 0.0
            for i = 1:max_ls
                ip.methods.r!(r, z, θ, 0.0)
                κ_vio_cand = bilinear_violation(ip, r, nquat = nquat)
                r_vio_cand = residual_violation(ip, r)
                if (r_vio_cand <= r_vio) || (κ_vio_cand <= κ_vio)
                    break
                end
                verbose && println("linesearch $i")
                # backtracking
                candidate_point!(z, s, z, Δ, -α * ls_scale^i)
            end
            # println("z: ", scn(norm(z), digits = 7))
            κ_vio = κ_vio_cand
            r_vio = r_vio_cand

            verbose && println("iter:", j,
                "  r: ", scn(norm(r, Inf)),
                "  r_vio: ", scn(r_vio),
                "  κ_vio: ", scn(κ_vio),
                "  Δ: ", scn(norm(Δ)),
                # "  Δ[ix]: ", scn(norm(Δ[ix])),
                # "  Δ[iy1]: ", scn(norm(Δ[iy1])),
                # "  Δ[iy2]: ", scn(norm(Δ[iy2])),
                # "  Δaff: ", scn(norm(Δaff)),
                # "  τ: ", scn(norm(ip.τ)),
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
        diff_sol && differentiate_solution!(ip, reg = max(ip.reg_val, κ_tol * γ_reg))
        return true
    else
        return false
    end
end

function rz!(ip::AbstractIPSolver, rz::AbstractMatrix{T}, z::AbstractVector{T},
        θ::AbstractVector{T}; reg = 0.0) where {T}
    z_reg = deepcopy(z)
    z_reg[ip.iy1] = max.(z[ip.iy1], reg)
    z_reg[ip.iy2] = max.(z[ip.iy2], reg)
    ip.methods.rz!(rz, z_reg, θ)
    return nothing
end

function general_correction_term!(r::AbstractVector{T}, Δ, ibil, idx_ineq, idx_soc,
        iy1, iy2; nquat::Int = 0) where {T}
    nc = Int(length(idx_soc) / 2)
    nsoc = Int(length(idx_soc) / 2)
    nineq = Int(length(idx_ineq) / 2)

    # Split between primals and duals
    idx_soc_p = idx_soc[1:nsoc]
    idx_soc_d = idx_soc[nsoc+1:2nsoc]
    idx_ineq_1 = intersect(iy1, idx_ineq)
    idx_ineq_2 = intersect(iy2, idx_ineq)

    n_soc = length(idx_soc_p)
    n_ort = length(idx_ineq_1)
    r[ibil[1:n_ort] .- nquat] .+= Δ[idx_ineq_1 .- nquat] .* Δ[idx_ineq_2 .- nquat]
    r[ibil[n_ort+1:end] .- nquat] .+= vcat(
        [second_order_cone_product(
            Δ[idx_soc_d[i] .- nquat],
            # Δη1[(i - 1) * ne .+ (1:ne)],
            Δ[idx_soc_p[i] .- nquat],
            # [Δs2[i]; Δb1[(i-1) * (ne - 1) .+ (1:(ne - 1))]]
        ) for i = 1:nc]...)
    return nothing
end

function least_squares!(ip::AbstractIPSolver, z::AbstractVector{T}, θ::AbstractVector{T},
        r::AbstractVector{T}, rz::AbstractMatrix{T}) where {T}
    # doing nothing gives the best result if z_t is correctly initialized with z_t-1 in th simulator
        # A = rz[[ip.idyn; ip.irst], [ip.ix; ip.iy1; ip.iy2]]
        # z[[ip.ix; ip.iy1; ip.iy2]] .+= A' * ((A * A') \ r[[ip.idyn; ip.irst]])
    return nothing
end

function initial_state!(z::AbstractVector{T}, iy1::SVector{ny,Int}, iy2::SVector{ny,Int},
        idx_ineq, idx_soc; ϵ::T=1e-20) where {T,ny}

    nc = Int(length(idx_soc) / 2)
    nsoc = Int(length(idx_soc) / 2)
    nineq = Int(length(idx_ineq) / 2)

    # Split between primals and duals
    idx_soc_p = idx_soc[1:nsoc]
    idx_soc_d = idx_soc[nsoc+1:2nsoc]
    idx_ineq_1 = intersect(iy1, idx_ineq)
    idx_ineq_2 = intersect(iy2, idx_ineq)

    n_soc = length(idx_soc_p)
    n_ort = length(idx_ineq_1)

    # ineq
    y1 = z[idx_ineq_1]
    y2 = z[idx_ineq_2]
    δy1 = max(-1.5 * minimum(y1), 0)
    δy2 = max(-1.5 * minimum(y2), 0)

    y1h = y1 .+ δy1
    y2h = y2 .+ δy2
    δhy1 = 0.5 * y1h'*y2h / (sum(y2h) + ϵ)
    δhy2 = 0.5 * y1h'*y2h / (sum(y1h) + ϵ)

    y10 = y1h .+ δhy1
    y20 = y2h .+ δhy2
    z[idx_ineq_1] .= y10
    z[idx_ineq_2] .= y20

    # soc
    for i = 1:nc
        e = [1; zeros(length(idx_soc_p[i]) - 1)] # identity element
        y1 = z[idx_soc_p[i]]
        y2 = z[idx_soc_d[i]]
        δy1 = max(-1.5 * (y1[1] - norm(y1[2:end])), 0)
        δy2 = max(-1.5 * (y2[1] - norm(y2[2:end])), 0)

        y1h = y1 + δy1 * e
        y2h = y2 + δy2 * e
        δhy1 = 0.5 * y1h'*y2h / ((y2h[1] + norm(y2h[2,end])) + ϵ)
        δhy2 = 0.5 * y1h'*y2h / ((y1h[1] + norm(y1h[2,end])) + ϵ)

        y10 = y1h + δhy1 * e
        y20 = y2h + δhy2 * e
        z[idx_soc_p[i]] .= y10
        z[idx_soc_d[i]] .= y20
    end
    return z
end

function interior_point_solve!(ip::AbstractIPSolver, z::AbstractVector{T}, θ::AbstractVector{T}) where T
    ip.z .= z
    ip.θ .= θ
    interior_point_solve!(ip)
end

function differentiate_solution!(ip::AbstractIPSolver; reg = 0.0)
    s = ip.s
    z = ip.z
    θ = ip.θ
    rz = ip.rz
    rθ = ip.rθ
    δz = ip.δz
    δzs = ip.δzs

    κ = ip.κ

    rz!(ip, rz, z, θ, reg = reg)
    rθ!(ip, rθ, z, θ)

    linear_solve!(ip.solver, δzs, rz, rθ, reg = reg)
    @inbounds @views @. δzs .*= -1.0
    mapping!(δz, s, δzs, z)

    nothing
end


################################################################################
# New methods
################################################################################

function residual_violation(ip::AbstractIPSolver, r::AbstractVector{T}; nquat::Int = 0) where {T}
    max(norm(r[ip.idyn[1:end-nquat]], Inf), norm(r[ip.irst .- nquat], Inf))
end

function general_centering(z::AbstractVector{T}, Δaff::AbstractVector{T},
		idx_ineq, idx_soc, iy1::SVector{n,Int}, iy2::SVector{n,Int}, αaff::T; nquat::Int = 0) where {n,T}
        # See Section 5.1.3 in CVXOPT
        # μ only depends on the dot products (no cone product)
        # The CVXOPT linear and quadratic cone program solvers
    nc = Int(length(idx_soc) / 2)
    nsoc = Int(length(idx_soc) / 2)
    nineq = Int(length(idx_ineq) / 2)

    # Split between primals and duals
    idx_soc_p = idx_soc[1:nsoc]
    idx_soc_d = idx_soc[nsoc+1:2nsoc]
    idx_ineq_1 = intersect(iy1, idx_ineq)
    idx_ineq_2 = intersect(iy2, idx_ineq)

    # ineq
    μ = z[idx_ineq_1]' * z[idx_ineq_2]
    μaff = (z[idx_ineq_1] - αaff * Δaff[idx_ineq_1 .- nquat])' * (z[idx_ineq_2] - αaff * Δaff[idx_ineq_2 .- nquat])
    # soc
    for i = 1:nc
        μ += z[idx_soc_p[i]]' * z[idx_soc_d[i]]
        μaff += (z[idx_soc_p[i]] - αaff * Δaff[idx_soc_p[i] .- nquat])' * (z[idx_soc_d[i]] - αaff * Δaff[idx_soc_d[i] .- nquat])
    end
    μ /= n
    μaff /= n

	σ = clamp(μaff / μ, 0.0, 1.0)^3
	return μ, σ
end

function bilinear_violation(ip::AbstractIPSolver, r::AbstractVector{T}; nquat::Int = 0) where {T}
    norm(r[ip.ibil .- nquat], Inf)
end

function general_bilinear_violation(z::AbstractVector{T}, idx_ineq, idx_soc, iy1, iy2) where {T}
    # USELESS
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
        # scalar part
        soc_vio_s = z[idx_soc_p[i]]' * z[idx_soc_d[i]]
        # vector part
        soc_vio_v = z[idx_soc_p[i][1]] * z[idx_soc_d[i][2:end]] + z[idx_soc_d[i][1]] * z[idx_soc_p[i][2:end]]
        κ_vio = max(κ_vio, abs(soc_vio_s))
        κ_vio = max(κ_vio, maximum(abs.(soc_vio_v)))
    end
    ineq_vio = maximum(abs.(z[idx_ineq_1] .* z[idx_ineq_2]))
    κ_vio = max(κ_vio, ineq_vio)
end

function soc_value(u::AbstractVector)
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
		idx_soc::AbstractVector{Vector{Int}}; τ::T=0.99, nquat::Int = 0, verbose::Bool = false) where {T}
        # We need to make this much more efficient (allocation free)
        # by choosing the right structure for idx_soc.
    α = 1.0
    for (i,idx) in enumerate(idx_soc)
        # we need -Δ here because we will taking the step x - α Δ
        α = min(α, soc_step_length(z[idx], -Δ[idx .- nquat], τ = τ, verbose = verbose))
        # α = min(α, soc_step_length(z[idx], -Δ[idx], τ = τ, verbose = verbose))
    end
    return α
end

function ineq_step_length(z::AbstractVector{T}, Δ::AbstractVector{T},
		idx_ineq::AbstractVector{Int}; τ::T=0.9995, nquat::Int=0) where {T}
        # We need to make this much more efficient (allocation free)
        # by choosing the right structure for idx_ineq.

    n = Int(length(idx_ineq) / 2)

    ατ_p = 1.0 # primal
    ατ_d = 1.0 # dual
    for i = 1:n
        ip = idx_ineq[i]
        id = idx_ineq[i + n]
        if Δ[ip .- nquat] > 0.0
            ατ_p = min(ατ_p, τ * z[ip] / Δ[ip .- nquat])
        end
        if Δ[id .- nquat] > 0.0
            ατ_d = min(ατ_d, τ * z[id] / Δ[id .- nquat])
        end
    end
    α = min(ατ_p, ατ_d)
    return α
end
