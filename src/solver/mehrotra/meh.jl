abstract type LinearSolver end

# residual
function r!(r, z, θ, κ)
    @warn "residual not defined"
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


# function r_update!(r, r̄)
#     r .= r̄
# end

# # optimization spaces
# abstract type Space end
#
# # Euclidean
# struct Euclidean <: Space
#     n::Int
# end
#
# function candidate_point!(z̄::Vector{T}, ::Euclidean, z::Vector{T}, Δ::Vector{T}, α::T) where T
#     z̄ .= z - α .* Δ
# end
#
# function update_point!(z::Vector{T}, ::Space, z̄::Vector{T}) where T
#     z .= z̄
# end
#
# function mapping!(δz, s::Euclidean, δzs, z) # TODO: make allocation free
#     δz .= δzs
# end


# interior-point solver options
@with_kw mutable struct Mehrotra27Options{T}
    r_tol::T = 1.0e-5
    κ_tol::T = 1.0e-5 # useless
    max_iter_inner::Int = 100
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

# mutable struct ResidualMethods
#     r!
#     rm!
#     rz!
#     rθ!
# end

# regularize Jacobian / Hessian
function regularize!(v_pr, v_du, reg_pr, reg_du)
    v_pr .+= reg_pr
    v_du .-= reg_du
end

mutable struct Mehrotra27{T}
    s::Space
    methods::ResidualMethods
    z::Vector{T}                 # current point
    Δaff::Vector{T}              # affine search direction
    Δ::Vector{T}                 # corrector search direction
    r                            # residual
    rm                           # corrector residual
    rbil                         # corrector residual
    r_merit::T                   # residual norm
    rz                           # residual Jacobian wrt z
    rθ                           # residual Jacobian wrt θ
    idx_ineq::Vector{Int}        # indices for inequality constraints
    idx_soc::Vector{Vector{Int}} # indices for second-order cone constraints
    idx_pr::Vector{Int}          # indices for primal variables
    idx_du::Vector{Int}          # indices for dual variables
    δz::Matrix{T}                # solution gradients (this is always dense)
    δzs::Matrix{T}               # solution gradients (in optimization space; δz = δzs for Euclidean)
    θ::Vector{T}                 # problem data
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
    iy1
    iy2
    reg_pr
    reg_du
    opts::Mehrotra27Options
end

function mehrotra(z, θ;
        s = Euclidean(length(z)),
        num_var = length(z),
        num_data = length(θ),
        idx_ineq = collect(1:0),
        idx_soc = Vector{Int}[],
        idx_pr = collect(1:s.n),
        idx_du = collect(1:0),
        iy1 = collect(1:0),
        iy2 = collect(1:0),
        ibil = collect(1:0),
        r! = r!, rm! = rm!, rz! = rz!, rθ! = rθ!,
        r  = zeros(s.n),
        rm = zeros(s.n),
        rz = spzeros(s.n, s.n),
        rθ = spzeros(s.n, num_data),
        reg_pr = [0.0], reg_du = [0.0],
        v_pr = view(rz, CartesianIndex.(idx_pr, idx_pr)),
        v_du = view(rz, CartesianIndex.(idx_du, idx_du)),
        opts = Mehrotra27Options()) where T

    rz!(rz, z, θ) # compute Jacobian for pre-factorization

    # Search direction
    Δaff = zeros(s.n)
    Δ = zeros(s.n)

    # Indices
    nbil = length(iy1)
    iy1 = SVector{nbil, Int}(iy1)
    iy2 = SVector{nbil, Int}(iy2)

    # Views
    z_y1 = view(z, iy1)
    z_y2 = view(z, iy2)
    Δaff_y1 = view(Δaff, iy1) # TODO this should be in Δ space
    Δaff_y2 = view(Δaff, iy2) # TODO this should be in Δ space
    Δ_y1 = view(Δ, iy1) # TODO this should be in Δ space
    Δ_y2 = view(Δ, iy2) # TODO this should be in Δ space
    rbil = bilinear_res(rm, ibil)

    Mehrotra27(
        s,
        ResidualMethods(r!, rm!, rz!, rθ!),
        z,
        Δaff,
        Δ,
        r,
        rm, # rm
        rbil,
        0.0,
        rz,
        rθ,
        idx_ineq,
        idx_soc,
        idx_pr,
        idx_du,
        zeros(length(z), num_data),
        zeros(s.n, num_data),
        θ,
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
        iy1,
        iy2,
        reg_pr,
        reg_du,
        opts)
end

function bilinear_res(r::AbstractVector, ibil)
    view(r, ibil)
end

function bilinear_res(r::RLin, ibil)
    r.rbil
end

# interior point solver
function interior_point_solve!(ip::Mehrotra27{T}) where T

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
    Δaff = ip.Δaff
    Δ = ip.Δ
    r = ip.r
    rm = ip.rm
    rbil = ip.rbil
    r_merit = ip.r_merit
    rz = ip.rz
    idx_ineq = ip.idx_ineq
    idx_soc = ip.idx_soc
    θ = ip.θ
    v_pr = ip.v_pr
    v_du = ip.v_du
    z_y1 = ip.z_y1
    z_y2 = ip.z_y2
    Δaff_y1 = ip.Δaff_y1
    Δaff_y2 = ip.Δaff_y2
    Δ_y1 = ip.Δ_y1
    Δ_y2 = ip.Δ_y2
    reg_pr = ip.reg_pr
    reg_du = ip.reg_du
    nbil = ip.nbil
    solver = ip.solver

    # initialize regularization
    reg_pr[1] = opts.reg_pr_init
    reg_du[1] = opts.reg_du_init

    # compute residual, residual Jacobian
    r!(r, z, θ, 0.0) # here we set κ = 0, Δ = 0
    r_merit = norm(r, res_norm)

    elapsed_time = 0.0

    for j = 1:max_iter_inner
        elapsed_time >= max_time && break
        elapsed_time += @elapsed begin
            # check for converged residual
            if r_merit < r_tol
                break
            end

            # compute residual Jacobian
            rz!(rz, z, θ)

            # regularize (fixed, TODO: adaptive)
            reg && regularize!(v_pr, v_du, reg_pr[1], reg_du[1])

            # compute step
            linear_solve!(solver, Δaff, rz, r)
            αaff = step_length(z_y1, z_y2, Δaff_y1, Δaff_y2, τ=1.0)
            μaff = (z_y1 - αaff * Δaff_y1)' * (z_y2 - αaff * Δaff_y2) / length(z_y1) # TODO too slow

            μ = z_y1'*z_y2 / nbil
            σ = (μaff / μ)^3

            r_update!(rm, r)
            rbil += Δaff_y1 .* Δaff_y2 .- σ*μ # TODO too slow
            linear_solve!(solver, Δ, rz, rm) # TODO this should not recompute the factorization of rz
            τ = progress(r_merit, ϵ_min = ϵ_min) # TODO too slow
            α = step_length(z_y1, z_y2, Δ_y1, Δ_y2; τ = τ)

            # candidate point
            candidate_point!(z, s, z, Δ, α)

            # update
            r!(r, z, θ, 0.0) # we set κ= 0.0 to measure the bilinear constraint violation
            r_merit = norm(r, res_norm)
            verbose && println("iter: ", j, "   res∞: ", scn(r_merit))
        end
    end

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

function progress(merit::T; ϵ_min::T=0.05) where {T<:Real}
    ϵ = min(ϵ_min, merit^2)
    τ = 1.0 - ϵ
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

function interior_point_solve!(ip::Mehrotra27{T}, z::AbstractVector{T}, θ::AbstractVector{T}) where T
    ip.z .= z
    ip.θ .= θ
    interior_point_solve!(ip)
end

function differentiate_solution!(ip::Mehrotra27)
    s = ip.s
    z = ip.z
    θ = ip.θ
    rz = ip.rz
    rθ = ip.rθ
    δz = ip.δz
    δzs = ip.δzs

    ip.methods.rz!(rz, z, θ) #TODO: maybe not needed
    ip.methods.rθ!(rθ, z, θ)

    linear_solve!(ip.solver, δzs, rz, rθ)
    @inbounds @views @. δzs .*= -1.0
    mapping!(δz, s, δzs, z)

    nothing
end

linear_solve!(solver::LinearSolver, x::Vector{T}, A::Array{T, 2}, b::Vector{T}) where T = linear_solve!(solver, x, sparse(A), b)
