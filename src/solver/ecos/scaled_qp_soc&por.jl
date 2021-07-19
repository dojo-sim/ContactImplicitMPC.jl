################################################################################
# Quadratic Program
################################################################################
using Random
using LinearAlgebra

"""
https://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf

min   1/2 * x' * P * x + a' * x
s.t.  G * x + s = h
      A * x = b
      s ≧ 0

Primal variables:
s ∈ [m]
x ∈ [n]
y ∈ [p]
z ∈ [m]

Problem data:
a ∈ [n]
h ∈ [m]
b ∈ [p]
P ∈ [n,n]
G ∈ [m,n]
A ∈ [p,n]

"""

################################################################################
# Structures and Constructors
################################################################################

abstract type AbstractConeVector36{T}
end

mutable struct PORVector36{T} <: AbstractConeVector36{T}
    v::AbstractVector{T}
end

mutable struct SOCVector36{T} <: AbstractConeVector36{T}
    v::AbstractVector{T}
end


abstract type AbstractCone36{T}
end

dim(c::AbstractCone36) = c.m

mutable struct PositiveOrthantCone36{T} <: AbstractCone36{T}
    m::Int # dimension of the cone
    s::PORVector36{T}
    z::PORVector36{T}
    λ::PORVector36{T}
    W::AbstractMatrix{T}
    Wi::AbstractMatrix{T}
end

function positive_orthant_cone(m::Int)
    s = PORVector36([1.0; 1e-2 * rand(m - 1) / sqrt(m)])
    z = PORVector36([1.0; 1e-2 * rand(m - 1) / sqrt(m)])
    λ = PORVector36([1.0; 1e-2 * rand(m - 1) / sqrt(m)])
    W = zeros(m,m)
    Wi = zeros(m,m)
    @assert value(s) > 0.0
    @assert value(z) > 0.0
    @assert value(λ) > 0.0
    return PositiveOrthantCone36(m, s, z, λ, W, Wi)
end

mutable struct SecondOrderCone36{T} <: AbstractCone36{T}
    m::Int # dimension of the cone
    s::SOCVector36{T}
    z::SOCVector36{T}
    λ::SOCVector36{T}
    W::AbstractMatrix{T}
    Wi::AbstractMatrix{T}
end

function second_order_cone(m::Int)
    s = SOCVector36([1.0; 1e-2 * rand(m - 1) / sqrt(m)])
    z = SOCVector36([1.0; 1e-2 * rand(m - 1) / sqrt(m)])
    λ = SOCVector36([1.0; 1e-2 * rand(m - 1) / sqrt(m)])
    W = zeros(m,m)
    Wi = zeros(m,m)
    @assert value(s) > 0.0
    @assert value(z) > 0.0
    @assert value(λ) > 0.0
    return SecondOrderCone36(m, s, z, λ, W, Wi)
end


mutable struct CartesianCone36
    nc::Int # number of individual cones
    m::Int # dimension of the cartesian product of all the cones
    ms::AbstractVector{Int} # dimension of each cone
    ind::AbstractVector{AbstractVector{Int}} # indices for each cone
    cones::AbstractVector{AbstractCone36}
end

function CartesianCone36(cones::AbstractVector{C}) where {C <: AbstractCone36}
    nc = length(cones)
    ms = dim.(cones)
    m = sum(ms)
    s = cumsum(ms)
    ind = [Vector(1 .+ s[i]-ms[i]:s[i]) for i=1:nc]

    return CartesianCone36(nc, m, ms, ind, cones)
end

function CartesianCone36(ms::AbstractVector{Int}, types::AbstractVector{Symbol})
    nc = length(ms)

    cones = Vector{AbstractCone36}(undef, nc)
    for i = 1:nc
        cones[i] = eval(types[i])(ms[i])
    end
    return CartesianCone36(cones)
end


mutable struct QuadraticProgram36{T}
    nc::Int
    n::Int
    p::Int
    m::Int
    ms::AbstractVector{Int}
    types::AbstractVector{Symbol}
    ix::AbstractVector{Int}
    iy::AbstractVector{Int}
    iz::AbstractVector{Int}
    a::AbstractVector{T}
    h::AbstractVector{T}
    b::AbstractVector{T}
    P::AbstractMatrix{T}
    G::AbstractMatrix{T}
    A::AbstractMatrix{T}
    M::AbstractMatrix{T}
end

function quadratic_program(n::Int, p::Int, ms::AbstractVector{Int},
        types::AbstractVector{Symbol}; seed::Int = 100)
    Random.seed!(seed)
    nc = length(ms)
    m = sum(ms)
    ix = Vector(1:n)
    iy = Vector(n .+ (1:p))
    iz = Vector(n+p .+ (1:m))

    a = rand(n)
    h = rand(m)
    b = rand(p)
    P_ = rand(n, n)
    P = P_ * P_'
    G = rand(m, n)
    A = rand(p, n)
    M = [P   A'         G'        ;
         A   zeros(p,p) zeros(p,m);
         G   zeros(m,p) zeros(m,m)]
    @assert rank(A) == p
    QuadraticProgram36(nc, n, p, m, ms, types, ix, iy, iz, a, h, b, P, G, A, M)
end


# interior-point solver options
@with_kw mutable struct QuadraticProgramOptions36{T}
    iter::Int = 10
    tol::T = 1.0e-8
    # γ captures the second order term for the correction step: γ = 1.0 fully captured, γ = 0.0 not captured
    γ::T = 1.0
    # Step size selection α ≈ τ * α_max
    τ_por::T = 1.0 - 1e-8 # recommend 0.99 however it is slow
    τ_soc::T = 1.0 - 1e-2 # recommend 0.99 however it is slow
end


mutable struct Solution36{T}
    n::Int
    p::Int
    m::Int
    x::AbstractVector{T}
    y::AbstractVector{T}
    cc::CartesianCone36
end

function solution(n::Int, p::Int, ms::AbstractVector{Int}, types::AbstractVector{Symbol})
    x = ones(n)
    y = ones(p)
    cc = CartesianCone36(ms, types)
    Solution36(n, p, m, x, y, cc)
end


mutable struct Residual36{T}
    n::Int
    p::Int
    m::Int
    rx::AbstractVector{T}
    ry::AbstractVector{T}
    rcc::CartesianCone36
    # only uses rcc.cones[i].z === rz
    # only uses rcc.cones[i].s === rs
end

function residual(n::Int, p::Int, ms::AbstractVector{Int}, types::AbstractVector{Symbol})
    rx = zeros(n)
    ry = zeros(p)
    rcc = CartesianCone36(ms, types)
    Residual36(n, p, m, rx, ry, rcc)
end


mutable struct Direction36{T}
    n::Int
    p::Int
    m::Int
    Δx::AbstractVector{T}
    Δy::AbstractVector{T}
    Δcc::CartesianCone36
    # only uses Δcc.cones[i].z === Δz
    # only uses Δcc.cones[i].s === Δs
end

function direction(n::Int, p::Int, ms::AbstractVector{Int}, types::AbstractVector{Symbol})
    Δx = zeros(n)
    Δy = zeros(p)
    Δcc = CartesianCone36(ms, types)
    Direction36(n, p, m, Δx, Δy, Δcc)
end

################################################################################
# Cone Methods
################################################################################

function por_value(u::AbstractVector)
    return minimum(u)
end

function soc_value(u::AbstractVector)
    u0 = u[1]
    u1 = u[2:end]
    return (u0^2 - u1' * u1)
end

function value(u::PORVector36)
    # if value >= 0 then u ≥ 0
    # if value >  0 then u > 0
    return por_value(u.v)
end

function value(u::SOCVector36)
    # if value >= 0 then u ≥ 0
    # if value >  0 then u > 0
    return soc_value(u.v)
end

function cone_prod(u::PORVector36, w::PORVector36)
    return u.v .* w.v
end

function cone_prod(u::SOCVector36, w::SOCVector36)
    return [u.v' * w.v; u.v[1] * w.v[2:end] + w.v[1] * u.v[2:end]]
end

function cone_div(u::PORVector36, w::PORVector36)
    return w.v ./ u.v
end

function cone_div(u::SOCVector36, w::SOCVector36)
    # check CVXOPT after equation 36b for inverse-free formula.
    m = length(u.v)
    A = [u.v[1]     u.v[2:end]';
         u.v[2:end] u.v[1]*I(m - 1)]
    return A \ w.v
end

function identity_vector(c::PositiveOrthantCone36{T}) where {T}
    e = ones(T, c.m)
    return e
end

function identity_vector(c::SecondOrderCone36{T}) where {T}
    e = [ones(T, 1) ; zeros(T, c.m - 1)]
    return e
end

function identity_vector(cs::CartesianCone36)
    m = cs.m
    e = vcat(identity_vector.(cs.cones)...)
    return e
end

function scaling!(cc::CartesianCone36)
    for i in eachindex(cc.cones)
        scaling!(cc.cones[i])
    end
    return nothing
end

function scaling!(c::PositiveOrthantCone36{T}) where {T}
    m = c.m
    W = Diagonal(sqrt.(c.s.v) ./ sqrt.(c.z.v))
    Wi = Diagonal(sqrt.(c.z.v) ./ sqrt.(c.s.v))
    λ = sqrt.(c.s.v .* c.z.v)

    # @show norm(λ - Wi * s.s, Inf)
    # @show norm(λ - W * s.z, Inf)
    @assert norm(λ - Wi * c.s.v, Inf) < 1e-8
    @assert norm(λ - W * c.z.v, Inf) < 1e-8

    c.λ.v = λ
    c.W = W
    c.Wi = Wi
    return nothing
end

function scaling!(c::SecondOrderCone36{T}) where {T}
    m = c.m
    e = identity_vector(c)
    J = Diagonal([1; - ones(m - 1)])

    z̄ = 1 / sqrt(c.z.v' * J * c.z.v + 1e-15) * c.z.v
    s̄ = 1 / sqrt(c.s.v' * J * c.s.v + 1e-15) * c.s.v

    γ = sqrt((1 + z̄' * s̄) / 2)
    w̄ = (1 / 2γ) * (s̄ + J * z̄)
    v = 1 / sqrt(2 * (w̄[1] + 1)) * (w̄ + e)

    W̄ = 2 * v * v' - J
    W̄i = 2 * J * v * v' * J - J

    κ = exp(1/4 * log((c.s.v' * J * c.s.v) / (c.z.v' * J * c.z.v)))
    W = κ * W̄
    Wi = 1/κ * W̄i

    λbar = [γ; 1/2 * (z̄[2:end] + s̄[2:end] + (z̄[1] - s̄[1]) / (w̄[1] + 1) * w̄[2:end])]
    λ = exp(1/4 * log((c.s.v' * J * c.s.v) * (c.z.v' * J * c.z.v))) * λbar

    @show norm(λ - Wi * c.s.v, Inf)
    @show norm(λ - W  * c.z.v, Inf)
    @show norm(z̄, Inf)
    @show norm(s̄, Inf)
    @show norm(W, Inf)
    @show norm(W', Inf)
    @show norm(W - W', Inf)
    # @show norm(inv(W) - Wi, Inf)

    @assert norm(W - W', Inf) < 1e-10
    # @assert norm(inv(W) - Wi, Inf) < 1e-8
    # @assert norm(λ - W  * c.z.v, Inf) < 1e-8
    # @assert norm(λ - Wi * c.s.v, Inf) < 1e-8

    c.λ.v = λ
    c.W = W
    c.Wi = Wi
    return nothing
end

function analytical_step_size(cc::CartesianCone36, Δcc::CartesianCone36;
        τ_por::T = 1-1e-8, τ_soc::T = 0.99) where {T}
    α = 1.0
    for i in eachindex(cc.ms)
        α = min(α, analytical_step_size(cc.cones[i], Δcc.cones[i], τ_por = τ_por, τ_soc = τ_soc))
    end
    return α
end

function analytical_step_size(c::PositiveOrthantCone36{T}, Δc::PositiveOrthantCone36{T};
        τ_por::T = 1-1e-8, τ_soc::T = 0.99) where {T}
    αs = por_analytical_step_size(c.s.v, Δc.s.v, τ = τ_por)
    αz = por_analytical_step_size(c.z.v, Δc.z.v, τ = τ_por)
    return min(αs, αz)
end

function analytical_step_size(c::SecondOrderCone36{T}, Δc::SecondOrderCone36{T};
        τ_por::T = 1-1e-8, τ_soc::T = 0.99) where {T}
    αs = soc_analytical_step_size(c.s.v, Δc.s.v, τ = τ_soc)
    αz = soc_analytical_step_size(c.z.v, Δc.z.v, τ = τ_soc)
    return min(αs, αz)
end

function por_analytical_step_size(λ::AbstractVector{T}, Δ::AbstractVector{T};
        τ::T = 1-1e-8) where {T}

    m = length(λ)
    α = 1.0
    for i = 1:m
        if Δ[i] < 0
            @warn "tweak"
            # α = min(α, -(1 - 1e-10) * λ[i] / Δ[i])
            α = min(α, - τ * λ[i] / Δ[i])
        end
    end
    @show α
    @show por_value(λ + α * Δ)
    return α
end

function soc_analytical_step_size(λ::AbstractVector{T}, Δ::AbstractVector{T};
        τ::T = 0.99) where {T}

    m = length(λ)
    J = Diagonal([1; - ones(m - 1)])
    λbar = λ / sqrt(λ' * J * λ + 1e-15)
    @warn "tweak"
    ρ = 1 / sqrt(λ' * J * λ) * [λbar' * J * Δ; Δ[2:end] - (λbar' * J * Δ + Δ[1]) / (λbar[1] + 1) * λbar[2:end]]
    α = τ * min(1.0, 1.0 / (norm(ρ[2:end]) - ρ[1]))
    # @show value(λ + α * Δ)
    return α
end

################################################################################
# Solver Methods
################################################################################

function residual!(r::Residual36, s::Solution36, qp::QuadraticProgram36)
    nc = qp.nc
    z_ = vcat([s.cc.cones[i].z.v for i = 1:nc]...)
    s_ = vcat([s.cc.cones[i].s.v for i = 1:nc]...)
    r.rx =       qp.P * s.x + qp.A' * s.y + qp.G' *  z_ + qp.a
    r.ry =       qp.A * s.x                             - qp.b
    rz   =  s_ + qp.G * s.x                             - qp.h

    ind = r.rcc.ind
    for i in eachindex(r.rcc.ms)
        r.rcc.cones[i].z.v = rz[ind[i]]
    end
    return nothing
end

function violation(r::Residual36, s::Solution36, qp::QuadraticProgram36)
    nc = qp.nc
    residual!(r, s, qp)
    v1 = norm(r.rx, Inf)
    v2 = norm(r.ry, Inf)
    v3 = norm(vcat([r.rcc.cones[i].z.v for i=1:nc]...), Inf)
    v4 = - min([value(s.cc.cones[i].s.v) for i=1:nc]..., 0.0)
    v5 = - min([value(s.cc.cones[i].z.v) for i=1:nc]..., 0.0)

    # W, Wi, λ = scaling(s) #changed
    v6 = norm(vcat([s.cc.cones[i].s.v' * s.cc.cones[i].z.v for i = 1:nc]...), Inf) #GOOD
    # v6 = norm(abs.(cone_prod(Wi' * s.s, W * s.z)), Inf) #GOOD
    # v7 = norm(abs.(cone_prod(s.s, s.z)), Inf) #BAD
    # @show scn(v7)

    v = [v1, v2, v3, v4, v5, v6]
    vmax, imax = findmax(v)
    return norm(v, Inf)
end

function affine_direction!(Δ::Direction36, r::Residual36, s::Solution36, qp::QuadraticProgram36)
    # W, Wi, λ = scaling(s) #changed
    nc = qp.nc

    dx = - r.rx
    dy = - r.ry
    # dz = - r.rz
    # ds = - cone_prod(λ, λ)
    @warn "changed"
    dz = - vcat([r.rcc.cones[i].z.v for i = 1:nc]...)
    @warn "changed"
    ds = - vcat([cone_prod(s.cc.cones[i].λ.v, s.cc.cones[i].λ.v) for i = 1:nc]...)

    ind = s.cc.ind
    @warn "changed"
    dz̃ = dz - vcat([s.cc.cones[i].W' * cone_div(s.cc.cones[i].λ.v, ds[ind[i]]) for i = 1:nc]...)
    # d_ = [dx; dy; dz - W' * cone_div(λ, ds)]
    d_ = [dx; dy; dz̃]

    Mreg = deepcopy(qp.M)
    for i = 1:nc
        Mreg[qp.iz[ind[i]], qp.iz[ind[i]]] += - s.cc.cones[i].W' * s.cc.cones[i].W
    end
    @show norm(Mreg)
    @show norm(Mreg[qp.iz, qp.iz])
    Δ_ = Mreg \ d_

    Δ.Δx = Δ_[qp.ix]
    Δ.Δy = Δ_[qp.iy]
    for i = 1:nc
        Δ.Δcc.cones[i].z.v = Δ_[qp.iz[ind[i]]]
        Δ.Δcc.cones[i].s.v = s.cc.cones[i].W' * (cone_div(s.cc.cones[i].λ.v, ds[ind[i]]) - s.cc.cones[i].W * Δ.Δcc.cones[i].z.v)
    end
    return nothing
end

function combined_direction!(Δ::Direction36, r::Residual36, s::Solution36,
        qp::QuadraticProgram36, σ::T; η::T = σ, γ::T = 1.0) where {T}
    nc = qp.nc

    # μ = (λ' * λ) / m
    μ = 0.0
    for i = 1:nc
        μ += (s.cc.cones[i].λ.v' * s.cc.cones[i].λ.v) / m
    end

    e = identity_vector(s.cc)

    dx = - (1 - η) * r.rx
    dy = - (1 - η) * r.ry
    # dz = - (1 - η) * r.rz
    # ds = - cone_prod(λ, λ) - γ * cone_prod(Wi' * Δ.Δs, W * Δ.Δz) + σ * μ .* e
    @warn "changed"
    dz = - (1 - η) * vcat([r.rcc.cones[i].z.v for i = 1:nc]...)
    @warn "changed"
    ds = - vcat([cone_prod(s.cc.cones[i].λ.v, s.cc.cones[i].λ.v) for i = 1:nc]...)
    ds -= γ * vcat([cone_prod(s.cc.cones[i].Wi' * Δ.Δcc.cones[i].s.v, s.cc.cones[i].W' * Δ.Δcc.cones[i].z.v) for i = 1:nc]...)
    @warn "we are treating all cones collectively, we could handle them individually here"
    ds += σ * μ * e

    ind = s.cc.ind
    @warn "changed"
    dz̃ = dz - vcat([s.cc.cones[i].W' * cone_div(s.cc.cones[i].λ.v, ds[ind[i]]) for i = 1:nc]...)
    # d_ = [dx; dy; dz - W' * cone_div(λ, ds)]
    d_ = [dx; dy; dz̃]

    Mreg = deepcopy(qp.M)
    # Mreg[qp.iz, qp.iz] += - W' * W
    for i = 1:nc
        Mreg[qp.iz[ind[i]], qp.iz[ind[i]]] += - s.cc.cones[i].W' * s.cc.cones[i].W
    end

    Δ_ = Mreg \ d_

    Δ.Δx = Δ_[qp.ix]
    Δ.Δy = Δ_[qp.iy]
    # Δ.Δz = Δ_[qp.iz]
    # Δ.Δs = W' * (cone_div(λ, ds) - W * Δ.Δz)
    for i = 1:nc
        Δ.Δcc.cones[i].z.v = Δ_[qp.iz[ind[i]]]
        Δ.Δcc.cones[i].s.v = s.cc.cones[i].W' * (cone_div(s.cc.cones[i].λ.v, ds[ind[i]]) - s.cc.cones[i].W * Δ.Δcc.cones[i].z.v)
    end
    return nothing
end

function centering(s::Solution36, Δ::Direction36)
    # n = s.n
    # p = s.p
    # m = s.m
    # W, Wi, λ = scaling(s) #changed

    α = analytical_step_size(s.cc, Δ.Δcc, τ_por = 1.0, τ_soc = 1.0)
    @show α
    @warn "issue here"
    # αs = analytical_step_size(s.s, Δ.Δs)
    # αz = analytical_step_size(s.z, Δ.Δz)
    # α = min(αs, αz)

    num = 0.0
    denum = 0.0
    for i = 1:s.cc.nc
        num += (s.cc.cones[i].Wi' * Δ.Δcc.cones[i].s.v)' * (s.cc.cones[i].W * Δ.Δcc.cones[i].z.v)
        denum += s.cc.cones[i].λ.v' * s.cc.cones[i].λ.v
    end

    ρ = 1 - α + α^2 * num / denum #good
    # ρ = 1 - α + α^2 * ((Wi' * Δ.Δs)' * (W * Δ.Δz)) / (λ' * λ) #good
    # ρ_ = ((s.s + α * Δ.Δs)' * (s.z + α * Δ.Δz)) / (s.s' * s.z) #good
    # @assert norm(ρ_ - ρ) < 1e-10
    σ = max(0, min(1, ρ))^3
    return α, ρ, σ
end

function update!(s::Solution36{T}, Δ::Direction36{T}; τ_por::T = 1-1e-8, τ_soc::T = 0.99) where {T}
    nc = s.cc.nc

    α = analytical_step_size(s.cc, Δ.Δcc, τ_por = τ_por, τ_soc = τ_soc)
    @show α
    s.x += α * Δ.Δx
    s.y += α * Δ.Δy
    for i = 1:nc
        s.cc.cones[i].z.v += α * Δ.Δcc.cones[i].z.v
        s.cc.cones[i].s.v += α * Δ.Δcc.cones[i].s.v
    end
    scaling!(s.cc)
    return nothing
end

function qp_solve!(qp::QuadraticProgram36, opts::QuadraticProgramOptions36)
    n = qp.n
    p = qp.p
    m = qp.m

    s = solution(n, p, qp.ms, qp.types)
    r = residual(n, p, qp.ms, qp.types)
    Δ = direction(n, p, qp.ms, qp.types)

    for k = 1:opts.iter
        residual!(r, s, qp)
        vio = violation(r, s, qp)
        println("iter: ", k, "   vio: ", scn(vio))
        (vio < opts.tol) && break

        affine_direction!(Δ, r, s, qp)
        α, ρ, σ = centering(s, Δ)
        combined_direction!(Δ, r, s, qp, σ, η = 0.0, γ = opts.γ)
        update!(s, Δ; τ_por = opts.τ_por, τ_soc = opts.τ_soc)
    end
    return nothing
end

################################################################################
# Problem Setup
################################################################################

n = 10
p = 3
ms = [2,3]
types = [:positive_orthant_cone, :second_order_cone]

qp = quadratic_program(n, p, ms, types)
s = solution(n, p, ms, types)
r = residual(n, p, ms, types)
Δ = direction(n, p, ms, types)

residual!(r, s, qp)
violation(r, s, qp)
scaling!(s.cc)
affine_direction!(Δ, r, s, qp)
α, ρ, σ = centering(s, Δ)
combined_direction!(Δ, r, s, qp, σ)
update!(s, Δ)

n = 10
p = 3
ms = [2,3]
types = [:positive_orthant_cone, :second_order_cone]
opts = QuadraticProgramOptions36(iter = 100, τ_por = 0.99, tol = 1e-8)
qp = quadratic_program(n, p, ms, types)
qp_solve!(qp, opts)



a = 10
a = 10
a = 10
a = 10
a = 10


################################################################################
# Test
################################################################################

m0 = 5
m1 = 10
c0 = positive_orthant_cone(m0)
c1 = second_order_cone(m1)

ms = [m0, m1]
types = [:positive_orthant_cone, :second_order_cone]
cs = CartesianCone36(ms, types)
@test typeof(cs.cones[1]) <: PositiveOrthantCone36
@test typeof(cs.cones[2]) <: SecondOrderCone36

e0 = ones(m0)
e1 = [1.0; zeros(m1 - 1)]
@test identity_vector(c0) == e0
@test identity_vector(c1) == e1
@test identity_vector(cs) == [e0; e1]

@test cs.ind == [[1,2,3,4,5,], [6,7,8,9,10,11,12,13,14,15]]


@testset "QP - SOC" begin
    m = 4

    function test_division(m::Int, seed::Int = 100)
        Random.seed!(100)
        u = [1; 1e-2 * rand(m - 1)]
        v = [1; 1e-2 * rand(m - 1)]
        # @show norm(cone_prod(u, cone_div(u, v)) - v)
        @test norm(cone_prod(u, cone_div(u, v)) - v) < 1e-10
        return nothing
    end

    test_division(m)


    function empirical_step_size(λ::AbstractVector, Δ::AbstractVector; β=0.99)
        α = 1.0
        val = value(λ + α * Δ)
        while val <= 1e-10
            α *= β
            val = value(λ + α * Δ)
        end
        @assert strictly_positive(λ + α * Δ)
        return α
    end

    Random.seed!(10)
    λ = [3.0; rand(m - 1)]
    Δ = 100 * rand(m)
    α0 = analytical_step_size(λ, Δ)
    @test α0 < 1.0
    @test abs(value(λ + α0 * Δ)) < 1e-4

    α1 = empirical_step_size(λ, Δ, β = 1-1e-6)
    @test abs(value(λ + α1 * Δ)) < 1e-4
    @test norm(α0 - α1) < 1e-4

    α2 = empirical_step_size(λ, Δ, β = 1/2)
    @test value(λ + α2 * Δ) > 0
end
