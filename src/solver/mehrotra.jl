# interior-point solver options
@with_kw mutable struct MehrotraOptions{T} <: AbstractIPOptions
    r_tol::T = 1.0e-5
    κ_tol::T = 1.0e-5
    max_iter::Int = 100
    max_time::T = 1e5
    diff_sol::Bool = false
    res_norm::Real = Inf
    reg::Bool = false
    ϵ_min = 0.05 # ∈ [0.005, 0.25]
        # smaller -> faster
        # larger  -> slower, more robust
    κ_reg = 1e-3 # bilinear constraint violation level at which regularization is triggered [1e-3, 1e-4]
    γ_reg = 1e-1 # regularization scaling parameters ∈ [0, 0.1]:
        # 0   -> faster & ill-conditioned
        # 0.1 -> slower & better-conditioned
    undercut::T = 5.0 # the solver will aim at reaching κ_vio = κ_tol / undercut
        # simulation choose undercut = Inf
        # MPC choose undercut = 5.0
    solver::Symbol = :lu_solver
    verbose::Bool = false
end

mutable struct Mehrotra{T,nx,ny,R,RZ,Rθ} <: AbstractIPSolver
    s::Space
    oss::OptimizationSpace13
    methods::ResidualMethods
    z::Vector{T}                 # current point
    Δaff::Vector{T}              # affine search direction
    Δ::Vector{T}                 # corrector search direction
    r::R                            # residual
    rz::RZ                           # residual Jacobian wrt z
    rθ::Rθ                           # residual Jacobian wrt θ

    iz                         # indices of z = (iw1, iort, isoc)
    iΔz                        # indices of Δz = (iw1, iort, isoc)
    ir                         # indices of the residual

    idx_ineq::Vector{Int}        # indices for inequality constraints
    idx_ort::Vector{Vector{Int}} # indices for inequality constraints split between primal and dual
    idx_orts::Vector{Vector{Int}} # indices for inequality constraints split between primal and dual
    idx_soc # indices for second-order cone constraints
    idx_socs # indices for second-order cone constraints
    δz::Matrix{T}                # solution gradients (this is always dense)
    δzs::Matrix{T}               # solution gradients (in optimization space; δz = δzs for Euclidean)
    θ::Vector{T}                 # problem data
    κ::Vector{T}                 # barrier parameter
    num_var::Int
    num_data::Int
    solver::LinearSolver
    ix::SVector{nx,Int}
    iy1::SVector{ny,Int}
    iy2::SVector{ny,Int}
    idyn::SVector{nx,Int}
    irst::SVector{ny,Int}
    ibil::SVector{ny,Int}
    reg_val::T
    τ::T
    σ::T
    μ::T
    αaff::T
    α::T
    iterations::Int
    opts::MehrotraOptions
end

function mehrotra(z::AbstractVector{T}, θ::AbstractVector{T};
        s = Euclidean(length(z)),
        oss = OptimizationSpace13(),
        num_var = length(z),
        num_data = length(θ),
        iz = nothing,
        iΔz = nothing,
        ir = nothing,
        idx_ineq = collect(1:0),
        idx_ort = [collect(1:0), collect(1:0)],
        idx_orts = [collect(1:0), collect(1:0)],
        idx_soc = Vector{Int}[],
        idx_socs = Vector{Int}[],
        ix = collect(1:0),
        iy1 = collect(1:0),
        iy2 = collect(1:0),
        idyn = collect(1:0),
        irst = collect(1:0),
        ibil = collect(1:0),
        r! = r!, rz! = rz!, rθ! = rθ!,
        r  = zeros(s.n),
        rz = spzeros(s.n, s.n),
        rθ = spzeros(s.n, num_data),
        opts::MehrotraOptions = MehrotraOptions()) where T

    rz!(rz, z, θ) # compute Jacobian for pre-factorization

    # Search direction
    Δaff = zeros(s.n)
    Δ = zeros(s.n)

    # Indices
    nx = length(ix)
    ny = length(iy1)
    ix = SVector{nx, Int}(ix)
    iy1 = SVector{ny, Int}(iy1)
    iy2 = SVector{ny, Int}(iy2)
    idyn = SVector{nx, Int}(idyn)
    irst = SVector{ny, Int}(irst)
    ibil = SVector{ny, Int}(ibil)
    ny == 0 && @warn "ny == 0, we will get NaNs during the Mehrotra solve."

    Ts = typeof.((r, rz, rθ))
    Mehrotra{T,nx,ny,Ts...}(
        s,
        oss,
        ResidualMethods(r!, rz!, rθ!),
        z,
        Δaff,
        Δ,
        r,
        rz,
        rθ,
        iz,
        iΔz,
        ir,
        idx_ineq,
        idx_ort,
        idx_orts,
        idx_soc,
        idx_socs,
        zeros(length(z), num_data),
        zeros(s.n, num_data),
        θ,
        zeros(1),
        num_var,
        num_data,
        eval(opts.solver)(rz),
        ix,
        iy1,
        iy2,
        idyn,
        irst,
        ibil,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0,
        opts)
end

function interior_point_solve!(ip::Mehrotra{T,nx,ny,R,RZ,Rθ}) where {T,nx,ny,R,RZ,Rθ}

    # space
    s = ip.s

    # options
    opts = ip.opts
    r_tol = opts.r_tol
    κ_tol = opts.κ_tol
    max_iter = opts.max_iter
    max_time = opts.max_time
    diff_sol = opts.diff_sol
    res_norm = opts.res_norm
    reg = opts.reg
    ϵ_min = opts.ϵ_min
    κ_reg = opts.κ_reg
    γ_reg = opts.γ_reg
    verbose = opts.verbose

    # unpack pre-allocated data
    z = ip.z
    Δaff = ip.Δaff
    Δ = ip.Δ
    r = ip.r
    rz = ip.rz
    idx_ineq = ip.idx_ineq
    idx_ort = ip.idx_ort
    idx_soc = ip.idx_soc
    θ = ip.θ
    κ = ip.κ
    ix = ip.ix
    iy1 = ip.iy1
    iy2 = ip.iy2
    solver = ip.solver

    iw1 = ip.iz[1]
    iort = ip.iz[2]
    isoc = ip.iz[3]
    iw1s = ip.iΔz[1]
    iorts = ip.iΔz[2]
    isocs = ip.iΔz[3]
    ibil_ort = ip.ir[5]
    ibil_soc = ip.ir[6]

    # Initialization
    ip.iterations = 0
    ip.reg_val = 0.0
    comp = false

    # compute residual, residual Jacobian
    ip.methods.r!(r, z, θ, 0.0) # here we set κ = 0, Δ = 0

    least_squares!(ip, z, θ, r, rz) # this one uses indices from global scope in nonlinear mode
    z .= initial_state!(z, ix, iy1, iy2)

    ip.methods.r!(r, z, θ, 0.0) # here we set κ = 0, Δ = 0

    r_vio = residual_violation(ip, r)
    κ_vio = bilinear_violation(ip, r)
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
            κ_vio = bilinear_violation(ip, r)
            ip.reg_val = κ_vio < κ_reg ? κ_vio * γ_reg : 0.0

            # compute residual Jacobian
            rz!(ip, rz, z, θ, reg = ip.reg_val)

            # compute affine search direction
            linear_solve!(solver, Δaff, rz, r, reg = ip.reg_val)
            αaff = step_length(z, Δaff, iy1, iy2; τ = 1.0)
            centering!(ip, z, Δaff, iy1, iy2, αaff)

            # Compute corrector residual
            ip.methods.r!(r, z, θ, max(ip.σ*ip.μ, κ_tol/opts.undercut)) # here we set κ = σ*μ, Δ = Δaff
            general_correction_term!(r, Δaff, ibil_ort, ibil_soc, iorts, isocs)

            # Compute corrector search direction
            linear_solve!(solver, Δ, rz, r, reg = ip.reg_val, fact=false)

            progress!(ip, max(r_vio, κ_vio), ϵ_min=ϵ_min)
            α = step_length(z, Δ, iy1, iy2; τ = ip.τ)

            # candidate point
            candidate_point!(z, s, z, Δ, α)

            # update
            ip.methods.r!(r, z, θ, 0.0) # we set κ= 0.0 to measure the bilinear constraint violation
            verbose && println("iter:", j,
                "  r: ", scn(norm(r, res_norm)),
                "  r_vio: ", scn(r_vio),
                "  κ_vio: ", scn(κ_vio),
                "  Δ: ", scn(norm(Δ)),
                "  Δaff: ", scn(norm(Δaff)),
                "  τ: ", scn(norm(ip.τ)),
                "  α: ", scn(norm(α)))

            r_vio = residual_violation(ip, r)
            κ_vio = bilinear_violation(ip, r)
        end
    end
    verbose && println("iter : ", ip.iterations,
                     "  r_vio: ", scn(r_vio),
                     "  κ_vio: ", scn(κ_vio),
                     )

    if (r_vio > r_tol) || (κ_vio > κ_tol)
        @show (r_vio > r_tol)
        @show (κ_vio > κ_tol)
        @error "Mehrotra solver failed to reduce residual below r_tol."
        return false
    end
    # @show R
    # @show scn.(κ)
    # @show scn(opts.κ_tol)
    # differentiate solution
    # @warn "wrong reg_val"
    # diff_sol && differentiate_solution!(ip, reg = ip.reg_val)
    # We regularize the Jacobian rz in the implicit function theorem to regularize the gradients.
    # We choose reg = κ_tol * γ as this is small (typically 1e-5) and always strictly positive,
    # avoiding NaN gradients when the bilinear constraints are satisfied perfectly rbil = 0.
    # I think that κ_tol * γ_reg > ip.reg_val is ALWAYS true.
    diff_sol && differentiate_solution!(ip, reg = max(ip.reg_val, κ_tol * γ_reg))
    return true
end

function initial_state!(z::AbstractVector{T}, ix::SVector{nx,Int},
        iy1::SVector{ny,Int}, iy2::SVector{ny,Int}; ϵ::T=1e-20) where {T,nx,ny}

    xt = z[ix]
    y1t = z[iy1]
    y2t = z[iy2]

    δy1 = max(-1.5 * minimum(y1t), 0)
    δy2 = max(-1.5 * minimum(y2t), 0)

    y1h = y1t .+ δy1
    y2h = y2t .+ δy2

    δhy1 = 0.5 * y1h'*y2h / (sum(y2h) + ϵ)
    δhy2 = 0.5 * y1h'*y2h / (sum(y1h) + ϵ)

    x0 = xt
    y10 = y1h .+ δhy1
    y20 = y2h .+ δhy2
    z[ix] .= x0
    z[iy1] .= y10
    z[iy2] .= y20
    return z
end

function progress!(ip::Mehrotra{T}, violation::T; ϵ_min::T = 0.05) where {T}
    ϵ = min(ϵ_min, violation^2)
    ip.τ = 1 - ϵ
end

function step_length(z::AbstractVector{T}, Δ::AbstractVector{T},
		iy1::SVector{n,Int}, iy2::SVector{n,Int}; τ::T=0.9995) where {n,T}
    ατ_p = 1.0
    ατ_d = 1.0
    for i in eachindex(iy1)
        if Δ[iy1[i]] > 0.0
            ατ_p = min(ατ_p, τ * z[iy1[i]] / Δ[iy1[i]])
        end
        if Δ[iy2[i]] > 0.0
            ατ_d = min(ατ_d, τ * z[iy2[i]] / Δ[iy2[i]])
        end
    end
    α = min(ατ_p, ατ_d)
    return α
end

function centering!(ip::Mehrotra{T}, z::AbstractVector{T}, Δaff::AbstractVector{T},
		iy1::SVector{n,Int}, iy2::SVector{n,Int}, αaff::T) where {n,T}
	μaff = (z[iy1] - αaff * Δaff[iy1])' * (z[iy2] - αaff * Δaff[iy2]) / n
	ip.μ = z[iy1]' * z[iy2] / n
	ip.σ = (μaff / ip.μ)^3
	return nothing
end
