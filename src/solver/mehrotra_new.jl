# interior-point solver options
@with_kw mutable struct MehrotraOptions{T} <: AbstractIPOptions
    r_tol::T = 1.0e-5
    κ_tol::T = 1.0e-5
    κ_init::T = 1.0                   # useless
    # κ_scale::T = 0.1                  # useless
    # ls_scale::T = 0.5                 # useless
    max_iter_inner::Int = 100
    # max_iter_outer::Int = 1           # useless
    # max_ls::Int = 50                  # useless
    max_time::T = 60.0
    diff_sol::Bool = false
    res_norm::Real = Inf
    reg::Bool = false
    reg_pr_init = 0.0
    reg_du_init = 0.0
    ϵ_min = 0.05 # ∈ [0.005, 0.25]
        # smaller -> faster
        # larger -> slower, more robust
    κ_reg = 1e-3 # bilinear constraint violation level at which regularization is triggered [1e-3, 1e-4]
    γ_reg = 1e-1 # regularization scaling parameters ∈ [0, 0.1]:
        # 0   -> faster & ill-conditioned
        # 0.1 -> slower & better-conditioned
    solver::Symbol = :lu_solver
    verbose::Bool = false
end

mutable struct Mehrotra{T,nx,ny,R,RZ,Rθ} <: AbstractIPSolver
    s::Space
    methods::ResidualMethods
    z::Vector{T}                 # current point
    Δaff::Vector{T}              # affine search direction
    Δ::Vector{T}                 # corrector search direction
    r::R                          # residual
    rbil                         # corrector residual
    r_merit::T                   # residual norm
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
    # nbil::Int #useless
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
    ibil
    reg_pr
    reg_du
    iterations::Int
    opts::MehrotraOptions
end

function mehrotra(z::AbstractVector{T}, θ::AbstractVector{T};
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
        ibil = collect(1:0),
        r! = r!, rm! = rm!, rz! = rz!, rθ! = rθ!,
        r  = zeros(s.n),
        rz = spzeros(s.n, s.n),
        rθ = spzeros(s.n, num_data),
        reg_pr = [0.0], reg_du = [0.0],
        v_pr = view(rz, CartesianIndex.(idx_pr, idx_pr)),
        v_du = view(rz, CartesianIndex.(idx_du, idx_du)),
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
    ibil = SVector{ny, Int}(ibil)
    ny == 0 && @warn "ny == 0, we will get NaNs during the Mehrotra solve."

    # Views
    z_y1 = view(z, iy1)
    z_y2 = view(z, iy2)
    Δaff_y1 = view(Δaff, iy1) # TODO this should be in Δ space
    Δaff_y2 = view(Δaff, iy2) # TODO this should be in Δ space
    Δ_y1 = view(Δ, iy1) # TODO this should be in Δ space
    Δ_y2 = view(Δ, iy2) # TODO this should be in Δ space
    rbil = bilinear_res(r, ibil)

    Ts = typeof.((r, rz, rθ))
    Mehrotra{T,nx,ny,Ts...}(
        s,
        ResidualMethods(r!, rm!, rz!, rθ!),
        z,
        # zeros(length(z)),
        Δaff,
        Δ,
        r,
        rbil,
        0.0,
        deepcopy(r), #useless
        # 0.0,
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
        ibil,
        reg_pr, reg_du,
        0,
        opts)
end

function bilinear_res(r::AbstractVector, ibil)
    view(r, ibil)
end

# interior point solver
function interior_point_solve!(ip::Mehrotra{T,nx,ny,R,RZ,Rθ}) where {T,nx,ny,R,RZ,Rθ}

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
    rbil = ip.rbil
    r_merit = ip.r_merit
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
    comp = false
    reg_val = 0.0

    # initialize regularization
    reg_pr[1] = opts.reg_pr_init
    reg_du[1] = opts.reg_du_init

    # δθ = θ - r.θ0
    # # comp = true
    # comp && println("**** δθ:", scn(norm(δθ), digits=4))
    # # comp = false
    # comp && println("****  θ[μ,h]:", scn.(θ[end-1:end], digits=4))
    # comp && println("****  θ:", scn(norm(θ), digits=4))
    # comp && println("****  z:", scn(norm(z), digits=4))

    # compute residual, residual Jacobian
    ip.methods.rm!(r, z, 0.0 .* Δaff, θ, 0.0) # here we set κ = 0, Δ = 0
    comp && println("**** rl:", scn(norm(r, res_norm), digits=4))

    least_squares!(ip.z, ip.θ, ip.r, ip.rz)
    z .= initial_state!(z, ix, iy1, iy2, comp = comp)

    ip.methods.rm!(r, z, 0.0 .* Δaff, θ, 0.0) # here we set κ = 0, Δ = 0
    comp && println("**** rinit:", scn(norm(r, res_norm), digits=4))

    r_merit = norm(r, res_norm)
    elapsed_time = 0.0

    for j = 1:max_iter_inner
        elapsed_time >= max_time && break
        elapsed_time += @elapsed begin
            # check for converged residual
            if r_merit < r_tol
                break
            end
            ip.iterations += 1
            comp && println("************************** ITERATION :", ip.iterations)

            # Compute regularization level
            κ_vio = maximum(r.rbil)
            κ_vio < κ_reg ? reg_val = κ_vio * γ_reg : reg_val = 0.0

            # compute residual Jacobian
            rz!(rz, z, θ, reg = reg_val)

            # regularize (fixed, TODO: adaptive)
            reg && regularize!(v_pr, v_du, reg_pr[1], reg_du[1])

            # compute affine search direction
            linear_solve!(solver, Δaff, rz, r, reg = reg_val)

            αaff = step_length(z, Δaff, iy1, iy2, τ=1.0)
            μaff = (z_y1 - αaff * Δaff[iy1])' * (z_y2 - αaff * Δaff[iy2]) / ny

            μ = z_y1'*z_y2 / length(z_y1)
            σ = (μaff / μ)^3

            # Compute corrector residual
            rm!(r, z, Δaff, θ, σ*μ) # here we set κ = σ*μ, Δ = Δaff
            # @show norm(Δaff)

            # Compute corrector search direction
            linear_solve!(solver, Δ, rz, r, reg = reg_val)

            τ = progress(r_merit, ϵ_min=ϵ_min)
            α = step_length(z, Δ, iy1, iy2, τ=τ)

            comp && println("**** Δ1:", scn(norm(α*Δ[ix]), digits=4))
            comp && println("**** Δ2:", scn(norm(α*Δ[iy1]), digits=4))
            comp && println("**** Δ3:", scn(norm(α*Δ[iy2]), digits=4))

            # verbose && println("iter:", j,
            #     "  r: ", scn(norm(r, res_norm)),
            #     "  Δ: ", scn(norm(Δ)),
            #     # "  Δ[ix]: ", scn(norm(Δ[ix])),
            #     # "  Δ[iy1]: ", scn(norm(Δ[iy1])),
            #     # "  Δ[iy2]: ", scn(norm(Δ[iy2])),
            #     "  Δaff: ", scn(norm(Δaff)),
            #     "  τ: ", scn(norm(τ)),
            #     "  α: ", scn(norm(α)))

            # candidate point
            candidate_point!(z, s, z, Δ, α)
            # update
            r!(r, z, θ, 0.0) # we set κ= 0.0 to measure the bilinear constraint violation
            r_merit = norm(r, res_norm)
        end
    end
    # verbose && println("iter : ", ip.iterations)
    # verbose && println("r final: ", scn(r_merit))

    if r_merit > r_tol
        @error "Mehrotra solver failed to reduce residual below r_tol."
        return false
    end

    # differentiate solution
    diff_sol && differentiate_solution!(ip, reg = reg_val)
    return true
end

# TODO maybe we will need to implement this merit function to use κ_tol > b and r_tol > a
# function merit(rlin::AbstractVector, rbil::AbstractVector, t::Real)
# 	a = norm(rlin, t)
# 	b = norm(rbil, t)
# 	return a, b
# end

function initial_state!(z, ix, iy1, iy2; comp::Bool=true)
    xt = z[ix]
    y1t = z[iy1]
    y2t = z[iy2]

    δy1 = max(-1.5 * minimum(y1t), 0)
    δy2 = max(-1.5 * minimum(y2t), 0)
    comp && println("**** δw2:", scn(norm(δy1), digits=4))
    comp && println("**** δw3:", scn(norm(δy2), digits=4))

    y1h = y1t .+ δy1
    y2h = y2t .+ δy2
    # comp && println("**** w2h:", scn.(y1h[1:3], digits=4))
    # comp && println("**** w3h:", scn.(y2h[1:3], digits=4))

    δhy1 = 0.5 * y1h'*y2h / sum(y2h)
    δhy2 = 0.5 * y1h'*y2h / sum(y1h)
    comp && println("****  dot:", scn(norm(0.5 * y1h'*y2h), digits=4))
    comp && println("**** δhw2:", scn(norm(δhy1), digits=4))
    comp && println("**** δhw3:", scn(norm(δhy2), digits=4))

    x0 = xt
    y10 = y1h .+ δhy1
    y20 = y2h .+ δhy2
    z[ix] .= x0
    z[iy1] .= y10
    z[iy2] .= y20
    return z
end

function progress(merit; ϵ_min=0.05)
    ϵ = min(ϵ_min, merit^2)
    τ = 1 - ϵ
    return τ
end

# function step_length(w2::S, w3::S, Δw2::S, Δw3::S; τ::Real=0.9995) where {S}
#     ατ_p = 1.0
#     ατ_d = 1.0
#     for i in eachindex(w2)
#         if Δw2[i] > 0.0
#             ατ_p = min(ατ_p, τ * w2[i] / Δw2[i])
#         end
#         if Δw3[i] > 0.0
#             ατ_d = min(ατ_d, τ * w3[i] / Δw3[i])
#         end
#     end
#     α = min(ατ_p, ατ_d)
#     return α
# end

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

function least_squares!(ip::Mehrotra{T,nx,ny,R,RZ,Rθ}) where {T,nx,ny,R,RZ,Rθ}
	least_squares!(ip.z, ip.θ, ip.r, ip.rz)
	return nothing
end

function least_squares!(z::Vector{T}, θ::AbstractVector{T}, r::RLin{T}, rz::RZLin{T}) where {T}
	δθ = r.θ0 - θ
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

function interior_point_solve!(ip::Mehrotra{T}, z::AbstractVector{T}, θ::AbstractVector{T}) where T
    ip.z .= z
    ip.θ .= θ
    interior_point_solve!(ip)
end

function differentiate_solution!(ip::Mehrotra; reg = 0.0)
    s = ip.s
    z = ip.z
    θ = ip.θ
    rz = ip.rz
    rθ = ip.rθ
    δz = ip.δz
    δzs = ip.δzs

    κ = ip.κ

    ip.methods.rz!(rz, z, θ, reg = reg)
    ip.methods.rθ!(rθ, z, θ)

    linear_solve!(ip.solver, δzs, rz, rθ, reg = reg)
    @inbounds @views @. δzs .*= -1.0
    mapping!(δz, s, δzs, z)

    nothing
end
