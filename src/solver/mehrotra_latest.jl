# interior-point solver options
@with_kw mutable struct Mehrotra113Options{T} <: AbstractIPOptions
    r_tol::T = 1.0e-5
    κ_tol::T = 1.0e-5
    max_iter_inner::Int = 100 #TODO rename
    max_time::T = 60.0
    diff_sol::Bool = false
    reg::Bool = false
    ϵ_min::T = 0.05 # ∈ [0.005, 0.25]
        # smaller -> faster
        # larger  -> slower, more robust
    κ_reg::T = 1e-3 # bilinear constraint violation level at which regularization is triggered [1e-3, 1e-4]
    γ_reg::T = 1e-0 # regularization scaling parameters ∈ [0, 0.1]:
        # 0   -> faster & ill-conditioned
        # 0.1 -> slower & better-conditioned
    ls_relaxation::Real = Inf # In the line search (ls); we accept a step if:
    # better_off = (ls_relaxation * r_vio >= r_vio_cand) || (ls_relaxation * κ_vio >= κ_vio_cand)
    # 1.0  ->  the nominal value for ls_relaxation is 1.0: we need to make progress at each step
    # Inf  ->  the extreme value for ls_relaxation is Inf: we always take the full step
    # even if this degrades all constraint satisfaction
    τ_ls::T = 2.0 # line search decay factor
    max_iter_ls::Int = 3 # maximum number of backtracking steps
    solver::Symbol = :lu_solver
    verbose::Bool = false


    res_norm::Real = Inf #useless
    κ_init::T = 1.0   # useless
    reg_pr_init = 0.0 #useless
    reg_du_init = 0.0 #useless
end

mutable struct Mehrotra113{T,nx,ny,R,RZ,Rθ} <: AbstractIPSolver
    s::Space
    methods::ResidualMethods
    z::Vector{T}                 # current point
    Δaff::Vector{T}              # affine search direction
    Δ::Vector{T}                 # corrector search direction
    r::R                            # residual
    r̄::R #useless                            # candidate residual
    rz::RZ                           # residual Jacobian wrt z
    rθ::Rθ                           # residual Jacobian wrt θ
    idx_ineq::Vector{Int}        # indices for inequality constraints
    idx_soc::Vector{Vector{Int}} # indices for second-order cone constraints
    idx_pr::Vector{Int}          # indices for primal variables
    idx_du::Vector{Int}          # indices for dual variables
    δz::Matrix{T}                # solution gradients (this is always dense)
    δzs::Matrix{T}               # solution gradients (in optimization space; δz = δzs for Euclidean)
    θ::Vector{T}                 # problem data
    κ::Vector{T}                 # barrier parameter
    num_var::Int
    num_data::Int
    solver::LinearSolver
    v_pr # view
    v_du # view
    z_y1 # view into z corresponding to the first set of variables in the bilinear constraints y1 .* y2 = 0 (/!\this is z space)
    z_y2 # view into z corresponding to the second set of variables in the bilinear constraints y1 .* y2 = 0 (/!\this is z space)
    Δaff_y1 # view into Δaff corresponding to the first set of variables in the bilinear constraints y1 .* y2 = 0 (/!\this is Δ space)
    Δaff_y2 # view into Δaff corresponding to the second set of variables in the bilinear constraints y1 .* y2 = 0 (/!\this is Δ space)
    Δ_y1 # view into Δ corresponding to the first set of variables in the bilinear constraints y1 .* y2 = 0 (/!\this is Δ space)
    Δ_y2 # view into Δ corresponding to the second set of variables in the bilinear constraints y1 .* y2 = 0 (/!\this is Δ space)
    ix::SVector{nx,Int}
    iy1::SVector{ny,Int}
    iy2::SVector{ny,Int}
    idyn::SVector{nx,Int}
    irst::SVector{ny,Int}
    ibil::SVector{ny,Int}
    reg_pr #useless
    reg_du #useless
    reg_val::T
    τ::T
    σ::T
    μ::T
    αaff::T
    α::T
    iterations::Int
    opts::Mehrotra113Options
end

function mehrotra_latest(z::AbstractVector{T}, θ::AbstractVector{T};
        s = Euclidean(length(z)),
        num_var = length(z),
        num_data = length(θ),
        idx_ineq = collect(1:0),
        idx_soc = Vector{Int}[],
        idx_pr = collect(1:s.n),
        idx_du = collect(1:0),
        ix = collect(1:0),
        iy1 = collect(1:0),
        iy2 = collect(1:0),
        idyn = collect(1:0),
        irst = collect(1:0),
        ibil = collect(1:0),
        r! = r!, rm! = rm!, rz! = rz!, rθ! = rθ!,
        r  = zeros(s.n),
        rz = spzeros(s.n, s.n),
        rθ = spzeros(s.n, num_data),
        reg_pr = [0.0], reg_du = [0.0],
        v_pr = view(rz, CartesianIndex.(idx_pr, idx_pr)),
        v_du = view(rz, CartesianIndex.(idx_du, idx_du)),
        opts::Mehrotra113Options = Mehrotra113Options()) where T

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
    ny == 0 && @warn "ny == 0, we will get NaNs during the Mehrotra113 solve."

    # Views
    z_y1 = view(z, iy1)
    z_y2 = view(z, iy2)
    Δaff_y1 = view(Δaff, iy1) # TODO this should be in Δ space
    Δaff_y2 = view(Δaff, iy2) # TODO this should be in Δ space
    Δ_y1 = view(Δ, iy1) # TODO this should be in Δ space
    Δ_y2 = view(Δ, iy2) # TODO this should be in Δ space

    Ts = typeof.((r, rz, rθ))
    Mehrotra113{T,nx,ny,Ts...}(
        s,
        ResidualMethods(r!, rm!, rz!, rθ!),
        z,
        Δaff,
        Δ,
        r,
        deepcopy(r), #useless
        rz,
        rθ,
        idx_ineq,
        idx_soc,
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
        z_y1,
        z_y2,
        Δaff_y1,
        Δaff_y2,
        Δ_y1,
        Δ_y2,
        ix,
        iy1,
        iy2,
        idyn,
        irst,
        ibil,
        reg_pr, reg_du,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0,
        opts,
        )
end

# interior point solver

function interior_point_solve!(ip::Mehrotra113{T}, z::AbstractVector{T}, θ::AbstractVector{T}) where T
    ip.z .= z
    ip.θ .= θ
    interior_point_solve!(ip)
end

function interior_point_solve!(ip::Mehrotra113{T,nx,ny,R,RZ,Rθ}) where {T,nx,ny,R,RZ,Rθ}

    # space
    s = ip.s

    # options
    opts = ip.opts
    r_tol = opts.r_tol
    κ_tol = opts.κ_tol
    max_iter_inner = opts.max_iter_inner
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
    idx_soc = ip.idx_soc
    θ = ip.θ
    κ = ip.κ
    v_pr = ip.v_pr
    v_du = ip.v_du
    z_y1 = ip.z_y1
    z_y2 = ip.z_y2
    Δaff_y1 = ip.Δaff_y1
    Δaff_y2 = ip.Δaff_y2
    Δ_y1 = ip.Δ_y1
    Δ_y2 = ip.Δ_y2
    ix = ip.ix
    iy1 = ip.iy1
    iy2 = ip.iy2
    reg_pr = ip.reg_pr
    reg_du = ip.reg_du
    solver = ip.solver

    # Initialization
    ip.iterations = 0
    ip.reg_val = 0.0
    comp = false

    # initialize regularization
    reg_pr[1] = opts.reg_pr_init
    reg_du[1] = opts.reg_du_init

    # compute residual, residual Jacobian
    ip.methods.rm!(r, z, 0.0 .* Δaff, θ, 0.0) # here we set κ = 0, Δ = 0
    comp && println("**** rl:", scn(norm(r, res_norm), digits=4))

    least_squares!(ip, z, θ, r, rz) # this one uses indices from global scope in nonlinear mode
    # zcopy = deepcopy(z)
    # z .= initial_state!(z, ix, iy1, iy2; comp = comp)
    # @warn "wrong initial state for SOC"
    z .= general_initial_state!(z, ix, iy1, iy2, idx_ineq, idx_soc)
    test_positivity(z, iy1, iy2, idx_ineq, idx_soc)
    # @show scn(norm(z - zcopy))

    ip.methods.rm!(r, z, 0.0 .* Δaff, θ, 0.0) # here we set κ = 0, Δ = 0
    comp && println("**** rinit:", scn(norm(r, res_norm), digits=4))

    r_vio = residual_violation(ip, r)
    # κ_vio = bilinear_violation(ip, r)
    @warn "changed bilinear violation"
    κ_vio = general_bilinear_violation(ip, z, idx_ineq, idx_soc, iy1, iy2)
    elapsed_time = 0.0

    for j = 1:max_iter_inner
        elapsed_time >= max_time && break
        elapsed_time += @elapsed begin
            # check for converged residual
            if (r_vio < r_tol) && (κ_vio < κ_tol)
                break
            end
            ip.iterations += 1
            comp && println("************************** ITERATION :", ip.iterations)

            # Compute regularization level
            # κ_vio = bilinear_violation(ip, r)
            # @warn "changed bilinear violation"
            κ_vio = general_bilinear_violation(ip, z, idx_ineq, idx_soc, iy1, iy2)
            ip.reg_val = κ_vio < κ_reg ? κ_vio * γ_reg : 0.0

            # compute residual Jacobian
            rz!(ip, rz, z, θ, reg = ip.reg_val)
            # display_structure(rz, ip.ibil, iy1, iy2, idx_ineq, idx_soc)

            # regularize (fixed, TODO: adaptive)
            reg && regularize!(v_pr, v_du, reg_pr[1], reg_du[1])

            # compute affine search direction
            @warn "recompute r!"
            ip.methods.r!(r, z, θ, κ_vio/10) # we set κ= 0.0 to measure the bilinear constraint violation
            linear_solve!(solver, Δaff, rz, r, reg = ip.reg_val)
            # @show scn(norm(rz * Δaff - r, Inf))

            # αaff = step_length(z, Δaff, iy1, iy2; τ = 1.0)
            αaff = general_step_length(z, Δaff, iy1, iy2, idx_ineq, idx_soc;
                    τ_ineq = 1.0, τ_soc = 1.0)
            # @show αaff_ - αaff
            ##########################################################################################
            ##########################################################################################
            @warn "now we only do the first step"

            r_vio = residual_violation(ip, r)
            # κ_vio = bilinear_violation(ip, r)
            @warn "changed bilinear violation"
            κ_vio = general_bilinear_violation(ip, z, idx_ineq, idx_soc, iy1, iy2)
            candidate_point!(z, s, z, Δaff, αaff)
            @show min_value(z, iy1, iy2, idx_ineq, idx_soc)
            test_positivity(z, iy1, iy2, idx_ineq, idx_soc)
            α = αaff
            ##########################################################################################
            #
            # centering!(ip, z, Δaff, iy1, iy2, αaff)
            #
            # # Compute corrector residual
            # # @warn "changed"
            # # @show "###### objective"
            # # @show ip.σ*ip.μ
            # # @show κ_tol/5
            # ip.methods.rm!(r, z, Δaff, θ, max(ip.σ*ip.μ, κ_tol/5)) # here we set κ = σ*μ, Δ = Δaff
            # # ip.methods.rm!(r, z, Δaff, θ, max(ip.σ*ip.μ, 0.0)) # here we set κ = σ*μ, Δ = Δaff
            # # println("μ: ", scn(ip.μ, digits=6))
            # # println("σ: ", scn(ip.σ, digits=6))
            #
            # # Compute corrector search direction
            # linear_solve!(solver, Δ, rz, r, reg = ip.reg_val, fact=false)
            # fraction_to_boundary!(ip, max(r_vio, κ_vio), ϵ_min=ϵ_min)
            # # α = step_length(z, Δ, iy1, iy2; τ = ip.τ)
            # # α = ineq_step_length(z, Δ, idx_ineq; τ = ip.τ) # too slow
            # # @warn "wrong fraction to boundary, should be capped at 0.99"
            # α = general_step_length(z, Δ, iy1, iy2, idx_ineq, idx_soc;
            #         τ_ineq = ip.τ, τ_soc = min(ip.τ, 0.99))
            #         # τ_ineq = ip.τ, τ_soc = min(ip.τ, 0.50))
            # # zcopy = deepcopy(z)
            # # test_value(z, iy1, iy2, idx_ineq, idx_soc)
            # # test_value(zcopy - α * Δ, iy1, iy2, idx_ineq, idx_soc)
            # # test_positivity(zcopy - α * Δ, iy1, iy2, idx_ineq, idx_soc)
            #
            # # @show α - α_
            # comp && println("**** Δ1:", scn(norm(α*Δ[ix]), digits=4))
            # comp && println("**** Δ2:", scn(norm(α*Δ[iy1]), digits=4))
            # comp && println("**** Δ3:", scn(norm(α*Δ[iy2]), digits=4))

            verbose && println("iter:", j,
                "  r: ", scn(norm(r, res_norm)),
                "  r_vio: ", scn(r_vio),
                "  κ_vio: ", scn(κ_vio),
                "  Δ: ", scn(norm(Δ)),
                # "  Δ[ix]: ", scn(norm(Δ[ix])),
                # "  Δ[iy1]: ", scn(norm(Δ[iy1])),
                # "  Δ[iy2]: ", scn(norm(Δ[iy2])),
                "  Δaff: ", scn(norm(Δaff)),
                "  τ: ", scn(norm(ip.τ)),
                "  α: ", scn(norm(α)))

            # # Line search
            # # candidate point with full step
            # # @show "@@@@@@@@@@@@@@@@@@@"
            # # @show α
            # candidate_point!(z, s, z, Δ, α)
            # # soc_correction(z, idx_soc, κ_tol)
            # test_positivity(z, iy1, iy2, idx_ineq, idx_soc)
            # r_vio_cand = 0.0
            # κ_vio_cand = 0.0
            # for ls = 1:ip.opts.max_iter_ls
            #     # check that we are making progress at least on one of the residual or bilinear constraints
            #     ip.methods.r!(r, z, θ, 0.0) # we set κ= 0.0 to measure the bilinear constraint violation
            #     r_vio_cand = residual_violation(ip, r)
            #     # κ_vio_cand = bilinear_violation(ip, r)
            #     # @warn "changed bilinear violation"
            #     κ_vio_cand = general_bilinear_violation(ip, z, idx_ineq, idx_soc, iy1, iy2)
            #
            #     better_off = (ip.opts.ls_relaxation * r_vio >= r_vio_cand) ||
            #         (ip.opts.ls_relaxation * κ_vio >= κ_vio_cand)
            #     # if we are better off after the step, we take it
            #     better_off && break
            #     # otherwise we backtrack
            #     candidate_point!(z, s, z, Δ, -α / ip.opts.τ_ls^ls)
            #     @show ls
            # end
            # test_positivity(z, iy1, iy2, idx_ineq, idx_soc)
            # # We set the violations to their updated values
            @warn "changed r!"
            ip.methods.r!(r, z, θ, 0.0) # we set κ= 0.0 to measure the bilinear constraint violation
            ip.methods.r!(r, z, θ, κ_vio/10) # we set κ= 0.0 to measure the bilinear constraint violation

            # r_vio = r_vio_cand
            # κ_vio = κ_vio_cand
        end
    end
    verbose && println("iter : ", ip.iterations,
                     "  r_vio: ", scn(r_vio),
                     "  κ_vio: ", scn(κ_vio),
                     )

    if (r_vio > r_tol) || (κ_vio > κ_tol)
        @error "Mehrotra113 solver failed to reduce residual below r_tol."
        return false
    end
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

# Harmful
# clipping!(z, iy1, iy2, max(ip.σ*ip.μ, κ_tol/5), κΣ = 1e10)
# function clipping!(z::Vector{T}, iy1, iy2, κ; κΣ = 1e10)
#     # z[iy2] .= max.( min.(z[iy2], κΣ * κ ./ z[iy1]) , κ ./ (κΣ * z[iy1]) )
#     return nothing
# end

function least_squares!(ip::Mehrotra113{T}, z::AbstractVector{T}, θ::AbstractVector{T},
        r::AbstractVector{T}, rz::AbstractMatrix{T}) where {T}
    # doing nothing gives the best result if z_t is correctly initialized with z_t-1 in th simulator
        # A = rz[[ip.idyn; ip.irst], [ip.ix; ip.iy1; ip.iy2]]
        # z[[ip.ix; ip.iy1; ip.iy2]] .+= A' * ((A * A') \ r[[ip.idyn; ip.irst]])
    return nothing
end

function residual_violation(ip::Mehrotra113{T}, r::AbstractVector{T}) where {T}
    max(norm(r[ip.idyn], Inf), norm(r[ip.irst], Inf))
end

function bilinear_violation(ip::Mehrotra113{T}, r::AbstractVector{T}) where {T}
    norm(r[ip.ibil], Inf)
end



function initial_state!(z::AbstractVector{T}, ix::SVector{nx,Int},
        iy1::SVector{ny,Int}, iy2::SVector{ny,Int}; comp::Bool=true, ϵ::T=1e-20) where {T,nx,ny}

    xt = z[ix]
    y1t = z[iy1]
    y2t = z[iy2]
    # comp && println("**** xt :", scn(norm(xt), digits=4))
    # comp && println("**** y1t :", scn(norm(y1t), digits=4))
    # comp && println("**** y2t :", scn(norm(y2t), digits=4))

    δy1 = max(-1.5 * minimum(y1t), 0)
    δy2 = max(-1.5 * minimum(y2t), 0)
    # comp && println("**** δw2:", scn(norm(δy1), digits=4))
    # comp && println("**** δw3:", scn(norm(δy2), digits=4))

    y1h = y1t .+ δy1
    y2h = y2t .+ δy2
    # comp && println("**** w2h:", scn.(y1h[1:3], digits=4))
    # comp && println("**** w3h:", scn.(y2h[1:3], digits=4))

    δhy1 = 0.5 * y1h'*y2h / (sum(y2h) + ϵ)
    δhy2 = 0.5 * y1h'*y2h / (sum(y1h) + ϵ)
    # comp && println("**** sum(y1h):", scn(sum(y1h), digits=4))
    # comp && println("**** sum(y2h):", scn(sum(y2h), digits=4))
    # comp && println("****  dot:", scn(norm(0.5 * y1h'*y2h), digits=4))
    # comp && println("**** δhw2:", scn(norm(δhy1), digits=4))
    # comp && println("**** δhw3:", scn(norm(δhy2), digits=4))

    x0 = xt
    y10 = y1h .+ δhy1
    y20 = y2h .+ δhy2
    z[ix] .= x0
    z[iy1] .= y10
    z[iy2] .= y20
    return z
end

function fraction_to_boundary!(ip::Mehrotra113{T}, violation::T; ϵ_min::T = 0.05) where {T}
    # Name comes from IPOPT paper Eq. 8
    # On the Implementation of an Interior-Point Filter Line-Search
    # Algorithm for Large-Scale Nonlinear Programming
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

function centering!(ip::Mehrotra113{T}, z::AbstractVector{T}, Δaff::AbstractVector{T},
		iy1::SVector{n,Int}, iy2::SVector{n,Int}, αaff::T) where {n,T}
        # using EQ (7) from CVXOPT, we have that s ∘ z  = μ e
        # for positive orthant: s .* z  .= μ
        # for second order cone: s' * z  = μ; s0 * z1 + z0 * s1 = 0
        # μ only depends on the dot products
        # See Section 5.1.3 in CVXOPT
        # The CVXOPT linear and quadratic cone program solvers
	μaff = (z[iy1] - αaff * Δaff[iy1])' * (z[iy2] - αaff * Δaff[iy2]) / n
    # @show "#########centering "
    # @show μaff
	ip.μ = z[iy1]' * z[iy2] / n
	ip.σ = clamp(μaff / ip.μ, 0.0, 1.0)^3
	return nothing
end

function rz!(ip::Mehrotra113{T}, rz::AbstractMatrix{T}, z::AbstractVector{T},
        θ::AbstractVector{T}; reg = 0.0) where {T}
    z_reg = deepcopy(z)
    z_reg[ip.iy1] = max.(z[ip.iy1], reg)
    z_reg[ip.iy2] = max.(z[ip.iy2], reg)
    ip.methods.rz!(rz, z_reg, θ)
    return nothing
end

function differentiate_solution!(ip::Mehrotra113; reg = 0.0)
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
# Linearized Solver
################################################################################

function least_squares!(ip::Mehrotra113{T}, z::Vector{T}, θ::AbstractVector{T},
		r::RLin{T}, rz::RZLin{T}) where {T}
	# @warn "wrong"
	δθ = θ - r.θ0
	δrdyn = r.rdyn0 - r.rθdyn * δθ
	δrrst = r.rrst0 - r.rθrst * δθ

	δw1 = rz.A1 * δrdyn + rz.A2 * δrrst
	δw2 = rz.A3 * δrdyn + rz.A4 * δrrst
	δw3 = rz.A5 * δrdyn + rz.A6 * δrrst

	@. @inbounds z[r.ix]  .= r.x0  .+ δw1
	@. @inbounds z[r.iy1] .= r.y10 .+ δw2
	@. @inbounds z[r.iy2] .= r.y20 .+ δw3
	return nothing
end

function residual_violation(ip::Mehrotra113{T}, r::RLin{T}) where {T}
    max(norm(r.rdyn, Inf), norm(r.rrst, Inf))
end

function bilinear_violation(ip::Mehrotra113{T}, r::RLin{T}) where {T}
    norm(r.rbil, Inf)
end

function rz!(ip::Mehrotra113{T}, rz::RZLin{T}, z::AbstractVector{T},
		θ::AbstractVector{T}; reg::T = 0.0) where {T}
	rz!(rz, z, θ, reg = reg)
	return nothing
end

function rz!(rz::RZLin{T}, z::AbstractVector{T},
		θ::AbstractVector{T}; reg::T = 0.0) where {T}
	rz!(rz, z, reg = reg)
	return nothing
end

function rθ!(ip::Mehrotra113{T}, rθ::RθLin{T}, z::AbstractVector{T},
		θ::AbstractVector{T}) where {T}
	return nothing
end

function rθ!(rθ::RθLin{T}, z::AbstractVector{T},
		θ::AbstractVector{T}) where {T}
	return nothing
end

################################################################################
# SOC code
################################################################################

function soc_correction(z, idx_soc, tol)
    nc = Int(length(idx_soc) / 2)
    for i = 1:2nc
        ne = length(idx_soc[i])
        v = soc_value(z[idx_soc[i]])
        # if v <= 0.0
        if v < tol/10
            @show "@@@@@@@@@@@@ soc_correction"
            @show v
            δ = tol/10 - v
            @show sqrt(δ)
            z[idx_soc[i]] .+= sqrt(δ)*[1.0; zeros(ne - 1)]
        end
    end
    return nothing
end

function display_structure(rz, ibil, iy1, iy2, idx_ineq, idx_soc)
    nc = Int(length(idx_soc) / 2)
    nsoc = Int(length(idx_soc) / 2)
    nineq = Int(length(idx_ineq) / 2)

    # Split between primals and duals
    idx_soc_p = idx_soc[1:nsoc]
    idx_soc_d = idx_soc[nsoc+1:2nsoc]
    idx_ineq_1 = intersect(iy1, idx_ineq)
    idx_ineq_2 = intersect(iy2, idx_ineq)

    plt = plot()
    plot!(Gray.(1e100*abs.(rz[Vector(ibil), Vector(iy1)])))
    plot!(Gray.(1e100*abs.(rz[Vector(ibil), Vector(iy2)])))
    display(plt)

    return nothing
end

function general_bilinear_violation(ip::Mehrotra113{T}, z::AbstractVector{T}, idx_ineq, idx_soc, iy1, iy2) where {T}
    # norm(r[ip.ibil], Inf)
    # @warn "changed κ_vio"
    nc = Int(length(idx_soc) / 2)
    nsoc = Int(length(idx_soc) / 2)
    nineq = Int(length(idx_ineq) / 2)

    # Split between primals and duals
    idx_soc_p = idx_soc[1:nsoc]
    idx_soc_d = idx_soc[nsoc+1:2nsoc]
    idx_ineq_1 = intersect(iy1, idx_ineq)
    idx_ineq_2 = intersect(iy2, idx_ineq)
    # @show "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    κ_vio = 0.0
    for i = 1:nc
        soc_vio = z[idx_soc_p[i]]' * z[idx_soc_d[i]] / length(idx_soc_d[i])
        # @show soc_vio
        κ_vio = max(κ_vio, soc_vio)
    end
    ineq_vio = maximum(abs.(z[idx_ineq_1] .* z[idx_ineq_2]))
    # @show ineq_vio
    κ_vio = max(κ_vio, ineq_vio)
end


function soc_value(u::AbstractVector)
    u0 = u[1]
    u1 = u[2:end]
    return (u0^2 - u1' * u1)
end

function test_positivity(z::AbstractVector{T}, iy1, iy2, idx_ineq::AbstractVector{Int},
        idx_soc::AbstractVector{Vector{Int}}, ϵ = 1e-13) where {T}

    nc = Int(length(idx_soc) / 2)
    nsoc = Int(length(idx_soc) / 2)
    nineq = Int(length(idx_ineq) / 2)

    # Split between primals and duals
    idx_soc_p = idx_soc[1:nsoc]
    idx_soc_d = idx_soc[nsoc+1:2nsoc]
    idx_ineq_1 = intersect(iy1, idx_ineq)
    idx_ineq_2 = intersect(iy2, idx_ineq)

    for i = 1:nc
        # @show soc_value(z[idx_soc_p[i]])
        # @show soc_value(z[idx_soc_d[i]])
        @assert soc_value(z[idx_soc_p[i]]) >= -ϵ #0.0
        @assert soc_value(z[idx_soc_d[i]]) >= -ϵ #0.0
    end
    @assert all(z[idx_ineq_1] .>= 0.0)
    @assert all(z[idx_ineq_2] .>= 0.0)

    return nothing
end

function test_value(z::AbstractVector{T}, iy1, iy2, idx_ineq::AbstractVector{Int},
        idx_soc::AbstractVector{Vector{Int}}) where {T}

    nc = Int(length(idx_soc) / 2)
    nsoc = Int(length(idx_soc) / 2)
    nineq = Int(length(idx_ineq) / 2)

    # Split between primals and duals
    idx_soc_p = idx_soc[1:nsoc]
    idx_soc_d = idx_soc[nsoc+1:2nsoc]
    idx_ineq_1 = intersect(iy1, idx_ineq)
    idx_ineq_2 = intersect(iy2, idx_ineq)

    for i = 1:nc
        @show scn(soc_value(z[idx_soc_p[i]]))
        @show scn(soc_value(z[idx_soc_d[i]]))
    end
    @show scn.(z[idx_ineq_1])
    @show scn.(z[idx_ineq_2])

    return nothing
end


function min_value(z::AbstractVector{T}, iy1, iy2, idx_ineq::AbstractVector{Int},
        idx_soc::AbstractVector{Vector{Int}}) where {T}

    nc = Int(length(idx_soc) / 2)
    nsoc = Int(length(idx_soc) / 2)
    nineq = Int(length(idx_ineq) / 2)

    # Split between primals and duals
    idx_soc_p = idx_soc[1:nsoc]
    idx_soc_d = idx_soc[nsoc+1:2nsoc]
    idx_ineq_1 = intersect(iy1, idx_ineq)
    idx_ineq_2 = intersect(iy2, idx_ineq)

    val = Inf
    for i = 1:nc
        val = min(val, soc_value(z[idx_soc_p[i]]))
        val = min(val, soc_value(z[idx_soc_d[i]]))
    end
    val = min(val, minimum(z[idx_ineq_1]))
    val = min(val, minimum(z[idx_ineq_2]))
    return val
end

function general_initial_state!(z::AbstractVector{T}, ix::SVector{nx,Int},
        iy1::SVector{ny,Int}, iy2::SVector{ny,Int}, idx_ineq, idx_soc;
        ϵ::T=1e-20) where {T,nx,ny}

    nc = Int(length(idx_soc) / 2)
    nsoc = Int(length(idx_soc) / 2)
    nineq = Int(length(idx_ineq) / 2)

    # Split between primals and duals
    idx_soc_p = idx_soc[1:nsoc]
    idx_soc_d = idx_soc[nsoc+1:2nsoc]
    idx_ineq_1 = intersect(iy1, idx_ineq)
    idx_ineq_2 = intersect(iy2, idx_ineq)
    @assert length(idx_ineq_1) == nineq
    @assert length(idx_ineq_2) == nineq

    # Equalities remain unchanged
    # z[ix] .= z[ix]

    # Inequalities
    y1t = z[idx_ineq_1]
    y2t = z[idx_ineq_2]

    δy1 = max(-1.5 * minimum(y1t), 0)
    δy2 = max(-1.5 * minimum(y2t), 0)
    y1h = y1t .+ δy1
    y2h = y2t .+ δy2

    δhy1 = 0.5 * y1h'*y2h / (sum(y2h) + ϵ)
    δhy2 = 0.5 * y1h'*y2h / (sum(y1h) + ϵ)
    y10 = y1h .+ δhy1
    y20 = y2h .+ δhy2

    z[idx_ineq_1] .= y10
    z[idx_ineq_2] .= y20

    # SOC
    for i = 1:nc
        n = length(idx_soc_p[i])
        @assert n == length(idx_soc_d[i])
        e = [1.0; zeros(n - 1)]
        # extract primal and dual cones
        p = z[idx_soc_p[i]]
        d = z[idx_soc_d[i]]

        # compute the minimal α such that
        # p + α e ≧ 0 with e the identity vector in the SOC
        αp = max(0.0, soc_step_length(p, e))
        αd = max(0.0, soc_step_length(d, e))
        println("αp : ", scn(αp))
        println("αd : ", scn(αd))

        # Ensure strict positivity
        δp = 1.5 * αp
        δd = 1.5 * αd
        ph = p + δp * e
        dh = d + δd * e
        println("δp : ", scn(δp))
        println("δd : ", scn(δd))
        println("ph : ", scn(soc_value(ph)))
        println("dh : ", scn(soc_value(dh)))

        # Equalizer
        δhp = 0.5 * p' * d  / (e' * d + ϵ)
        δhd = 0.5 * p' * d  / (e' * p + ϵ)
        p0 = ph + δhp * e
        d0 = dh + δhd * e
        println("δhp : ", scn(δhp))
        println("δhd : ", scn(δhd))
        println("p0 : ", scn(soc_value(p0)))
        println("d0 : ", scn(soc_value(d0)))

        z[idx_soc_p[i]] .= p0
        z[idx_soc_d[i]] .= d0
        println("z[p] : ", scn(soc_value(z[idx_soc_p[i]])))
        println("z[d] : ", scn(soc_value(z[idx_soc_d[i]])))
    end

    for i = 1:nc
        @assert soc_value(z[idx_soc_p[i]]) >= 0.0
        @assert soc_value(z[idx_soc_d[i]]) >= 0.0
    end
    @assert all(z[idx_ineq_1] .>= 0.0)
    @assert all(z[idx_ineq_2] .>= 0.0)

    return z
end

function general_step_length(z::AbstractVector{T}, Δ::AbstractVector{T}, iy1, iy2,
		idx_ineq::AbstractVector{Int}, idx_soc::AbstractVector{Vector{Int}};
        τ_ineq::T=0.9995, τ_soc::T=0.99) where {n,T}
    # We need to make this much more efficient (allocation free)
    # by choosing the right structure for idx_ine and idx_soc.

    α_ineq = ineq_step_length(z, Δ, idx_ineq; τ = τ_ineq)
    @warn "changed α_soc"
    # α_soc = soc_step_length(z, Δ, idx_soc; τ = τ_soc)
    α = α_ineq
    # @show α_ineq
    while min_value(z - α * Δ, iy1, iy2, idx_ineq, idx_soc) < 0.0 && (α > 1e-17)
        α *= 0.5
    end

    # @show "##### general_step_length"
    # @show α_ineq
    # @show α_soc
    # α = min(α_ineq, α_soc)
    return α
end










# function soc_step_size_slow(λ::AbstractVector{T}, Δ::AbstractVector{T};
#         τ::T = 0.99) where {T}
#     m = length(λ)
#     J = Diagonal([1; - ones(m - 1)])
#     λbar = λ / sqrt(λ' * J * λ + 1e-15)
#     ρ = 1 / sqrt(λ' * J * λ) * [λbar' * J * Δ; Δ[2:end] - (λbar' * J * Δ + Δ[1]) / (λbar[1] + 1) * λbar[2:end]]
#     α = τ * min(1.0, 1.0 / (norm(ρ[2:end]) - ρ[1]))
#     return α
# end
#


# function soc_value(u::AbstractVector)
#     u0 = u[1]
#     u1 = u[2:end]
#     return (u0^2 - u1' * u1)
# end
#
Random.seed!(100)
m = 5
λ = [3; 1e-1rand(m - 1)]
Δ = 1e2rand(m)
α1 = soc_step_length(λ, Δ, τ = 1.0, ϵ = 0.0)
val1 = soc_value(λ + α1 * Δ)
α2 = soc_step_length(λ, Δ, τ = 1.0, ϵ = 1e-10)
val2 = soc_value(λ + α2 * Δ)
α3 = soc_step_length(λ, Δ, τ = 0.99)
val3 = soc_value(λ + α3 * Δ)

@test norm(val1) < 1e-10
@test norm(val2) < 1e-8

################################################################################
# Mehrotra113 Structure
################################################################################

# STRUCT
#
# Mehrotra113
# Mehrotra113Options
#
#
# METHODS
#
# mehrotra
# interior_point_solve!
# least_squares!
# residual_violation
# bilinear_violation
# initial_state!
# fraction_to_boundary!
# step_length!
# centering!
# interior_point_solve!
# rz!
# differentiate_solution!
