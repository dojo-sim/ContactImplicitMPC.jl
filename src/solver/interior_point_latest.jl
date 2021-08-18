abstract type LinearSolver end
abstract type AbstractIPSolver end
abstract type AbstractIPOptions end

function interior_point_options(ip_type::Symbol)
    if ip_type == :interior_point
        return Symbol("InteriorPoint116Options")
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
@with_kw mutable struct InteriorPoint116Options{T} <: AbstractIPOptions
    r_tol::T = 1.0e-5
    κ_tol::T = 1.0e-5
    ls_scale::T = 0.5
    max_iter_inner::Int = 100
    max_ls::Int = 3
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

mutable struct InteriorPoint116{T} <: AbstractIPSolver
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
    opts::InteriorPoint116Options

    κ::Vector{T} # useless
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
        opts::InteriorPoint116Options = InteriorPoint116Options()) where T

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

    InteriorPoint116(
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
function interior_point_solve!(ip::InteriorPoint116{T}) where T

    # space
    s = ip.s

    # methods
    # r! = ip.methods.r!
    # rm! = ip.methods.rm!
    # rz! = ip.methods.rz!
    # rθ! = ip.methods.rθ!

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

    # TODO need least squares init /!|

    # compute residual, residual Jacobian
    ip.methods.r!(r, z, θ, 0.0)
    κ_vio = general_bilinear_violation(z, idx_ineq, idx_soc, iy1, iy2)
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
            κ_vio = general_bilinear_violation(z, idx_ineq, idx_soc, iy1, iy2)
            ip.reg_val = κ_vio < κ_reg ? κ_vio * γ_reg : 0.0

            # compute residual Jacobian
            # TODO need to use reg_val here
            # @warn "changed"
            ip.methods.rz!(rz, z, θ)
            # TODO rz!(ip, rz, z, θ, reg = ip.reg_val) # this is not adapted to the second order cone

            # regularize (fixed, TODO: adaptive)
            reg && regularize!(v_pr, v_du, reg_pr[1], reg_du[1])

            # compute step
            # TODO need to use reg_val here
            linear_solve!(solver, Δ, rz, r)

            α_ineq = ineq_step_length(z, Δ, idx_ineq; τ = 1.0)
            α_soc = soc_step_length(z, Δ, idx_soc; τ = 1.0)
            α = min(α_ineq, α_soc)

            μ, σ = centering(z, Δ, iy1, iy2, α)
            αaff = α

            # Compute corrector residual
            ip.methods.r!(r, z, θ, max(σ*μ, κ_tol/5)) # here we set κ = σ*μ, Δ = Δaff
            general_correction_term!(r, Δ, ibil, idx_ineq, idx_soc, iy1, iy2)

            # Compute corrector search direction
            # TODO need to use reg_val here
            linear_solve!(solver, Δ, rz, r)
            τ = max(0.95, 1 - max(r_vio, κ_vio)^2)
            α_ineq = ineq_step_length(z, Δ, idx_ineq; τ = τ)
            α_soc = soc_step_length(z, Δ, idx_soc; τ = min(τ, 0.99))
            α = min(α_ineq, α_soc)

            # reduce norm of residual
            candidate_point!(z, s, z, Δ, α)
            κ_vio_cand = 0.0
            r_vio_cand = 0.0
            for i = 1:max_ls
                ip.methods.r!(r, z, θ, 0.0)
                κ_vio_cand = general_bilinear_violation(z, idx_ineq, idx_soc, iy1, iy2)
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

            verbose && println(
                "in:", j,
                "   αaff:", scn(αaff, digits = 0),
                "   α:", scn(α, digits = 0),
                "   μσ:", scn(μ*σ, digits = 0),
                "   κ_vio:", scn(κ_vio, digits = 0),
                "   r_vio:", scn(r_vio, digits = 0),
                )
        end
    end
    if (r_vio < r_tol) && (κ_vio < κ_tol)
        # differentiate solution
        diff_sol && differentiate_solution!(ip)
        return true
    else
        return false
    end
end

function rz!(ip::InteriorPoint116{T}, rz::AbstractMatrix{T}, z::AbstractVector{T},
        θ::AbstractVector{T}; reg = 0.0) where {T}
    z_reg = deepcopy(z)
    z_reg[ip.iy1] = max.(z[ip.iy1], reg)
    z_reg[ip.iy2] = max.(z[ip.iy2], reg)
    ip.methods.rz!(rz, z_reg, θ)
    return nothing
end

function general_correction_term!(r::AbstractVector{T}, Δ, ibil, idx_ineq, idx_soc, iy1, iy2) where {T}
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

    r[ibil[1:n_ort]] .+= Δ[idx_ineq_1] .* Δ[idx_ineq_2]
    r[ibil[n_ort+1:end]] .+= vcat(
        [second_order_cone_product(
            Δ[idx_soc_d[i]],
            # Δη1[(i - 1) * ne .+ (1:ne)],
            Δ[idx_soc_p[i]],
            # [Δs2[i]; Δb1[(i-1) * (ne - 1) .+ (1:(ne - 1))]]
        ) for i = 1:nc]...)
    return nothing
end

function residual_mehrotra(model::ContactModel, env::Environment{<:World, LinearizedCone}, z, Δz, θ, κ)
	ix, iy1, iy2 = linearization_var_index(model, env)
	idyn, irst, ibil, ialt = linearization_term_index(model, env)
	rm = residual(model, env, z, θ, κ)
	rm[ibil] .+= Δz[iy1] .* Δz[iy2]
	return rm
end

function residual_mehrotra(model::ContactModel, env::Environment{<:World, NonlinearCone}, z, Δz, θ, κ)
	nc = model.dim.c
	nb = nc * friction_dim(env)
	ne = dim(env)

	idyn, irst, ibil, ialt = linearization_term_index(model, env)
	Δq2, Δγ1, Δb1, Δη1, Δs1, Δs2 = unpack_z(model, env, Δz)
	rm = residual(model, env, z, θ, κ)

	# Positive orthant
	rm[ibil[1:nc]] .+= Δγ1 .* Δs1
	# SOC
	rm[ibil[nc .+ (1:nc + nb)]] .+= vcat(
		[second_order_cone_product(
			Δη1[(i - 1) * ne .+ (1:ne)],
			[Δs2[i]; Δb1[(i-1) * (ne - 1) .+ (1:(ne - 1))]
		]) for i = 1:nc]...)
	return rm
end



function interior_point_solve!(ip::InteriorPoint116{T}, z::AbstractVector{T}, θ::AbstractVector{T}) where T
    ip.z .= z
    ip.θ .= θ
    interior_point_solve!(ip)
end

function differentiate_solution!(ip::InteriorPoint116)
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

function residual_violation(ip::InteriorPoint116{T}, r::AbstractVector{T}) where {T}
    max(norm(r[ip.idyn], Inf), norm(r[ip.irst], Inf))
end

function centering(z::AbstractVector{T}, Δaff::AbstractVector{T},
		iy1::SVector{n,Int}, iy2::SVector{n,Int}, αaff::T) where {n,T}
        # See Section 5.1.3 in CVXOPT
        # μ only depends on the dot products (no cone product)
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
    λ_λ = λ0^2 - λ[2:end]' * λ[2:end]
    verbose && println(
        "    vλ:", scn(soc_value(λ), digits = 0, exp_digits = 2),
        "    vλ+Δ:", scn(soc_value(λ+Δ), digits = 0, exp_digits = 2),
        "    λ_λ: ", scn(λ_λ, digits = 0, exp_digits = 2),
        "    λ:", scn.(λ, digits = 0, exp_digits = 2),
        "    Δ:", scn.(Δ, digits = 0, exp_digits = 2),
        )
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
    verbose && println(
        "     α:", scn(α, digits = 0, exp_digits = 2))
    verbose && cone_plot(λ, Δ, show_α = true)
    return α
end

function soc_step_length(z::AbstractVector{T}, Δ::AbstractVector{T},
		idx_soc::AbstractVector{Vector{Int}}; τ::T=0.99, verbose::Bool = false) where {T}
        # We need to make this much more efficient (allocation free)
        # by choosing the right structure for idx_soc.
    α = 1.0
    for (i,idx) in enumerate(idx_soc)
        # we need -Δ here because we will taking the step x - α Δ
        α = min(α, soc_step_length(z[idx], -Δ[idx], τ = τ, verbose = false))
        # α = min(α, soc_step_length(z[idx], -Δ[idx], τ = τ, verbose = verbose))
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

function cone_plot(λ::AbstractVector{T}, Δ::AbstractVector{T};
		ϵ = 1e-20, N::Int = 100, verbose::Bool = false, show_α::Bool = false) where {T}
	plt = plot(aspect_ratio = 1.0, fmt = :svg)

	# unit circle
	θ = Vector(0:2π/N:(1.0+1/N)*2π)
	circ = [cos.(θ) , sin.(θ)]
	plot!(plt, circ[1], circ[2], label = false)

	# Keypoints
	λp = project(λ)
	λΔp = project(λ + Δ)
	verbose && @show soc_value(λ )
	verbose && @show soc_value(λ + Δ)
	scatter!(plt, λΔp[1:1],  λΔp[2:2],  markersize = 4.0, color = :orange, label = "λ+Δ")
	scatter!(plt, λp[1:1],   λp[2:2],   markersize = 8.0, color = :red,    label = "λ")

	if show_α
		α = soc_step_length(λ, Δ, τ = 1.0, ϵ = ϵ)
		λαΔp = project(λ + α * Δ)
		verbose && @show α
		verbose && @show soc_value(λ + α * Δ)
		scatter!(plt, λαΔp[1:1], λαΔp[2:2], markersize = 2.0, color = :black,   label = "λ+αΔ")
	end

	# Line
	A = Vector(0:1/N:1.0)
	L = [[project(λ + a * Δ)[1] for a in A], [project(λ + a * Δ)[2] for a in A]]
	plot!(plt, L[1], L[2], label = false)

	display(plt)
	return nothing
end

function project(x::AbstractVector{T}; ϵ::T = 1e-15) where {T}
	if length(x) == 3
		project_3D(x, ϵ = ϵ)
	elseif length(x) == 2
		project_2D(x, ϵ = ϵ)
	end
end

function project_3D(x::AbstractVector{T}; ϵ::T = 1e-15) where {T}
	@assert length(x) == 3
	xs = x[1]
	xv = x[2:end]
	r = norm(xv) / (xs + ϵ)
	v = xv ./ norm(xv)
	p = r * v
	return p
end

function project_2D(x::AbstractVector{T}; ϵ::T = 1e-15) where {T}
	@assert length(x) == 2
	xs = x[1]
	xv = x[2:end]
	r = norm(xv) / (xs + ϵ)
	v = xv ./ norm(xv)
	p = r * v
	return [p[1], 0.0]
end

#
# λ = [1 + 1e-30, 1/sqrt(2), 1/sqrt(2)]
# Δ = [2.0,-1,2]
# cone_plot(λ, Δ, verbose = true)
