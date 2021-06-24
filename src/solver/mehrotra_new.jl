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
    ϵ_min = 0.05
    solver::Symbol = :lu_solver
    verbose::Bool = false
end

mutable struct Mehrotra{T} <: AbstractIPSolver
    s::Space
    methods::ResidualMethods
    z::Vector{T}                 # current point
    # z̄::Vector{T}                 # candidate point
    Δaff::Vector{T}              # affine search direction
    Δ::Vector{T}                 # corrector search direction
    r                            # residual
    rm                           # corrector residual
    rbil                         # corrector residual
    r_merit::T                   # residual norm
    r̄ #useless                            # candidate residual
    # r̄_merit::T                   # candidate residual norm
    rz                           # residual Jacobian wrt z
    rθ                           # residual Jacobian wrt θ
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
    nbil::Int
    solver::LinearSolver
    v_pr # view
    v_du # view
    z_y1 # view into z corresponding to the first set of variables in the bilinear constraints y1 .* y2 = 0 (/!\this is z space)
    z_y2 # view into z corresponding to the second set of variables in the bilinear constraints y1 .* y2 = 0 (/!\this is z space)
    Δaff_y1 # view into Δaff corresponding to the first set of variables in the bilinear constraints y1 .* y2 = 0 (/!\this is Δ space)
    Δaff_y2 # view into Δaff corresponding to the second set of variables in the bilinear constraints y1 .* y2 = 0 (/!\this is Δ space)
    Δ_y1 # view into Δ corresponding to the first set of variables in the bilinear constraints y1 .* y2 = 0 (/!\this is Δ space)
    Δ_y2 # view into Δ corresponding to the second set of variables in the bilinear constraints y1 .* y2 = 0 (/!\this is Δ space)
    ix
    iy1
    iy2
    ibil
    reg_pr
    reg_du
    iterations::Int
    opts::MehrotraOptions
end

function mehrotra(z, θ;
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
        rm = deepcopy(r),
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
    nbil = length(iy1)
    iy1 = SVector{nbil, Int}(iy1)
    iy2 = SVector{nbil, Int}(iy2)
    ibil = SVector{nbil, Int}(ibil)
    nbil == 0 && @warn "nbil == 0, we will get NaNs during the Mehrotra solve."

    # Views
    z_y1 = view(z, iy1)
    z_y2 = view(z, iy2)
    Δaff_y1 = view(Δaff, iy1) # TODO this should be in Δ space
    Δaff_y2 = view(Δaff, iy2) # TODO this should be in Δ space
    Δ_y1 = view(Δ, iy1) # TODO this should be in Δ space
    Δ_y2 = view(Δ, iy2) # TODO this should be in Δ space
    rbil = bilinear_res(rm, ibil)

    Mehrotra(
        s,
        ResidualMethods(r!, rm!, rz!, rθ!),
        z,
        # zeros(length(z)),
        Δaff,
        Δ,
        r,
        rm, # rm
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
        nbil,
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
function interior_point_solve!(ip::Mehrotra{T}) where T

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
    verbose = opts.verbose

    # unpack pre-allocated data
    z = ip.z
    # z̄ = ip.z̄
    Δaff = ip.Δaff
    Δ = ip.Δ
    r = ip.r
    rm = ip.rm
    rbil = ip.rbil
    nbil = ip.nbil
    r_merit = ip.r_merit
    # r̄ = ip.r̄
    # r̄_merit = ip.r̄_merit
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
    ip.iterations = 0
    comp = false

    # initialize regularization
    reg_pr[1] = opts.reg_pr_init
    reg_du[1] = opts.reg_du_init

    δθ = θ - r.θ0
    comp = true
    comp && println("**** δθ:", scn(norm(δθ), digits=4))
    comp && println("****  θ[μ,h]:", scn.(θ[end-1:end], digits=4))
    comp && println("****  θ:", scn(norm(θ), digits=4))
    comp && println("****  z:", scn(norm(z), digits=4))
    comp = false
    # compute residual, residual Jacobian
    # @warn "bad init"

    # @warn "bad κ"
    ip.methods.rm!(r, z, 0.0 .* Δaff, θ, 0.0) # here we set κ = 0, Δ = 0
    # ip.methods.rm!(r, z, 0.0 .* Δaff, θ, 1e-8) # here we set κ = 0, Δ = 0
    comp && println("**** rl:", scn(norm(r, res_norm), digits=4))

	nx = length(r.ix)
	ny = length(r.iy1)
	A = [r.Dx r.Dy1 zeros(nx,ny);
		 r.Rx r.Ry1 Diagonal(r.Ry2);
         # zeros(ny,nx) Diagonal(r.y2) Diagonal(r.y1);
		 ]
    comp && println("****  A:", scn(norm(A), digits=4))
	x = [r.rdyn0; r.rrst0] - [r.rθdyn * δθ; r.rθrst * δθ]
    comp && println("****  x:", scn(norm(x), digits=4))
	wt = A' * ((A * A') \ x)
    comp && println("**** wt:", scn(norm(wt), digits=4))
    comp && println("**** Awt-x:", scn(norm(A*wt .- x), digits=4))
    # z[[ix; iy1; iy2]] .= z[[ix; iy1; iy2]] + wt
	z[[ix; iy1; iy2]] .= [r.x0; r.y10; r.y20] + wt
    comp && println("**** z+wt:", scn(norm(z), digits=4))

    z .= initial_state!(z, ix, iy1, iy2, comp = comp)
    comp && println("**** w1:", scn(norm(z[ix]), digits=4))
    comp && println("**** w2:", scn(norm(z[iy1]), digits=4))
    comp && println("**** w3:", scn(norm(z[iy2]), digits=4))

    ip.methods.rm!(r, z, 0.0 .* Δaff, θ, 0.0) # here we set κ = 0, Δ = 0
    comp && println("**** rinit:", scn(norm(r, res_norm), digits=4))

    # plt = plot()
    # # plot!(plt, abs.(r))
    # plot!(plt, abs.([Vector(r.rdyn); Vector(r.rrst); Vector(r.rbil); ]))
    # display(plt)

    r_merit = norm(r, res_norm)
    # @show r_merit
    elapsed_time = 0.0

    for j = 1:max_iter_inner
        elapsed_time >= max_time && break
        elapsed_time += @elapsed begin

            # @show j
            # @show norm(z)
            # @show norm(θ)
            # @show minimum(abs.(z))
            # plt = plot()
            # plot!(z)
            # display(plt)

            # check for converged residual
            if r_merit < r_tol
                break
            end
            ip.iterations += 1
            # compute residual Jacobian
            rz!(rz, z, θ)
            # @show norm(rz.Dx)
            # @show norm(rz.Dy1)
            # @show norm(rz.Rx)
            # @show norm(rz.Ry1)
            # @show norm(rz.Ry2)
            # @show norm(rz.y1)
            # @show norm(rz.y2)
            # @show minimum(rz.y1)
            # @show minimum(rz.y2)
            # @show minimum(abs.(rz.y2))
            # @show minimum(abs.(rz.y1))

            # regularize (fixed, TODO: adaptive)
            reg && regularize!(v_pr, v_du, reg_pr[1], reg_du[1])

            # compute affine search direction
            linear_solve!(solver, Δaff, rz, r)
            A_ = [rz.Dx rz.Dy1 zeros(length(ix), length(iy1));
                  rz.Rx rz.Ry1 Diagonal(rz.Ry2)]
            M_ = [rz.Dx rz.Dy1 zeros(length(ix), length(iy1));
                  rz.Rx rz.Ry1 Diagonal(rz.Ry2);
                  zeros(length(iy2), length(ix)) Diagonal(rz.y2) Diagonal(rz.y1)]
            # @warn "overwriting the linear_solve!"
            Δaff[[ix; iy1; iy2]] = M_ \ [r.rdyn; r.rrst; r.rbil]

            comp && println("**** A:", scn(norm(A_), digits=4))
            comp && println("**** M:", scn(norm(M_), digits=4))
            comp && println("**** cond(Dx):", scn(cond(rz.Dx), digits=4))
            comp && println("**** My1:", scn(norm(rz.y1), digits=4))
            comp && println("**** My2:", scn(norm(rz.y2), digits=4))
            comp && println("**** My1:", scn.(rz.y1, digits=1))
            comp && println("**** My2:", scn.(rz.y2, digits=1))
            comp && println("**** rdyn:", scn(norm(r.rdyn), digits=4))
            comp && println("**** rrst:", scn(norm(r.rrst), digits=4))
            comp && println("**** rbil:", scn(norm(r.rbil), digits=4))

            comp && println("**** Δaff1:", scn(norm(Δaff[ix]), digits=4))
            comp && println("**** Δaff2:", scn(norm(Δaff[iy1]), digits=4))
            comp && println("**** Δaff3:", scn(norm(Δaff[iy2]), digits=4))
            comp && println("**** Δaff2:", scn.(Δaff[iy1], digits=1))
            comp && println("**** Δaff3:", scn.(Δaff[iy2], digits=1))

            err = M_ * [Δaff[ix]; Δaff[iy1]; Δaff[iy2]] - [r.rdyn; r.rrst; r.rbil]
            err2 = M_ * [Δaff[ix]; Δaff[iy1]; Δaff[iy2]] - [r.rdyn; 0.0*r.rrst; 0.0*r.rbil]
            comp && println("**** err:", scn(norm(err), digits=4))
            comp && println("**** err2:", scn(norm(err2), digits=4))
            comp && println("**** errx:", scn(norm(err[ix]), digits=4))
            comp && println("**** erry1:", scn(norm(err[iy1]), digits=4))
            comp && println("**** erry2:", scn(norm(err[iy2]), digits=4))
            comp && println("**** err:", scn.(err, digits=1))

            # @show norm(Δaff)
            # plt = plot()
            # plot!(Δaff)
            # display(plt)
            αaff = step_length(z_y1, z_y2, Δaff_y1, Δaff_y2, τ=1.0)
            μaff = (z_y1 - αaff * Δaff[iy1])' * (z_y2 - αaff * Δaff[iy2]) / nbil
            # @show nbil
            # @show norm(αaff)
            # @show norm(μaff)

            μ = z_y1'*z_y2 / length(z_y1)
            σ = (μaff / μ)^3

            # Compute corrector residual
            # rm!(rm, z, Δaff, θ, σ*μ) # here we set κ = σ*μ, Δ = Δaff
            #@@@
            # rm!(rm, z, Δaff, θ, max(min(1e-8, r_tol/10), σ*μ)) # here we set κ = σ*μ, Δ = Δaff
            rm!(rm, z, Δaff, θ, σ*μ) # here we set κ = σ*μ, Δ = Δaff
            # @show norm(Δaff)


            # Compute corrector search direction
            linear_solve!(solver, Δ, rz, rm)

            # @warn "overwriting the linear_solve!"
            Δ[[ix; iy1; iy2]] = M_ \ [rm.rdyn; rm.rrst; rm.rbil]

            τ = progress(r_merit, ϵ_min=ϵ_min)
            α = step_length(z_y1, z_y2, Δ_y1, Δ_y2, τ=τ)
            comp && println("**** Δ1:", scn(norm(α*Δ[ix]), digits=4))
            comp && println("**** Δ2:", scn(norm(α*Δ[iy1]), digits=4))
            comp && println("**** Δ3:", scn(norm(α*Δ[iy2]), digits=4))
            #
            # @show norm(rm.rdyn)
            # @show norm(rm.rrst)
            # @show norm(rm.rbil)
            # @show norm(Δ)
            # @show τ
            # @show α
            # @show norm(r, res_norm)
            # @show norm(Δ)
            # @show norm(Δaff)
            # @show norm(τ)
            # @show norm(α)
            verbose && println("iter:", j,
                "  r: ", scn(norm(r, res_norm)),
                "  Δ: ", scn(norm(Δ)),
                # "  Δ[ix]: ", scn(norm(Δ[ix])),
                # "  Δ[iy1]: ", scn(norm(Δ[iy1])),
                # "  Δ[iy2]: ", scn(norm(Δ[iy2])),
                "  Δaff: ", scn(norm(Δaff)),
                "  τ: ", scn(norm(τ)),
                "  α: ", scn(norm(α)))

            # candidate point
            candidate_point!(z, s, z, Δ, α)

            # update
            # @warn "bad approx"
            r!(r, z, θ, 0.0) # we set κ= 0.0 to measure the bilinear constraint violation
            # r!(r, z, θ, 1e-10) # we set κ= 0.0 to measure the bilinear constraint violation # this κ=1e-10 remove instabilities leading to divergence of the method.

            # plt = plot()
            # # plot!(plt, abs.(r))
            # plot!(plt, abs.([Vector(r.rdyn); Vector(r.rrst); Vector(r.rbil); ]))
            # display(plt)

            r_merit = norm(r, res_norm)
            # @show r_merit
            # verbose && println("iter: ", j, "   res∞: ", scn(r_merit))
        end
    end
    verbose && println("r final: ", scn(r_merit))

    if r_merit > r_tol
        @error "Mehrotra solver failed to reduce residual below r_tol."
        return false
    end

    # differentiate solution
    diff_sol && differentiate_solution!(ip)
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
    comp && println("**** w2h:", scn.(y1h[1:3], digits=4))
    comp && println("**** w3h:", scn.(y2h[1:3], digits=4))

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

function step_length(w2::S, w3::S, Δw2::S, Δw3::S; τ::Real=0.9995) where {S}
    ατ_p = 1.0
    ατ_d = 1.0
    for i in eachindex(w2)
        if Δw2[i] > 0.0
            ατ_p = min(ατ_p, τ * w2[i] / Δw2[i])
        end
        if Δw3[i] > 0.0
            ατ_d = min(ατ_d, τ * w3[i] / Δw3[i])
        end
    end
    α = min(ατ_p, ατ_d)
    return α
end

function interior_point_solve!(ip::Mehrotra{T}, z::AbstractVector{T}, θ::AbstractVector{T}) where T
    ip.z .= z
    ip.θ .= θ
    interior_point_solve!(ip)
end

function differentiate_solution!(ip::Mehrotra)
    s = ip.s
    z = ip.z
    θ = ip.θ
    rz = ip.rz
    rθ = ip.rθ
    δz = ip.δz
    δzs = ip.δzs

    κ = ip.κ

    ip.methods.rz!(rz, z, θ)
    ip.methods.rθ!(rθ, z, θ)

    linear_solve!(ip.solver, δzs, rz, rθ)
    @inbounds @views @. δzs .*= -1.0
    mapping!(δz, s, δzs, z)

    nothing
end
