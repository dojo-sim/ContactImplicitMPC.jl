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

abstract type AbstractConeVector34{T}
end

mutable struct PORVector34{T} <: AbstractConeVector34{T}
    v::AbstractVector{T}
end

mutable struct SOCVector34{T} <: AbstractConeVector34{T}
    v::AbstractVector{T}
end


abstract type AbstractCone34{T}
end

dim(c::AbstractCone34) = c.m

mutable struct PositiveOrthantCone34{T} <: AbstractCone34{T}
    m::Int # dimension of the cone
    s::PORVector34{T}
    z::PORVector34{T}
    λ::PORVector34{T}
    W::AbstractMatrix{T}
    Wi::AbstractMatrix{T}
end

function positive_orthant_cone(m::Int)
    s = PORVector34([1.0; 1e-2 * rand(m - 1) / sqrt(m)])
    z = PORVector34([1.0; 1e-2 * rand(m - 1) / sqrt(m)])
    λ = PORVector34([1.0; 1e-2 * rand(m - 1) / sqrt(m)])
    W = zeros(m,m)
    Wi = zeros(m,m)
    @assert value(s) > 0.0
    @assert value(z) > 0.0
    @assert value(λ) > 0.0
    return PositiveOrthantCone34(m, s, z, λ, W, Wi)
end

mutable struct SecondOrderCone34{T} <: AbstractCone34{T}
    m::Int # dimension of the cone
    s::SOCVector34{T}
    z::SOCVector34{T}
    λ::SOCVector34{T}
    W::AbstractMatrix{T}
    Wi::AbstractMatrix{T}
end

function second_order_cone(m::Int)
    s = SOCVector34([1.0; 1e-2 * rand(m - 1) / sqrt(m)])
    z = SOCVector34([1.0; 1e-2 * rand(m - 1) / sqrt(m)])
    λ = SOCVector34([1.0; 1e-2 * rand(m - 1) / sqrt(m)])
    W = zeros(m,m)
    Wi = zeros(m,m)
    @assert value(s) > 0.0
    @assert value(z) > 0.0
    @assert value(λ) > 0.0
    return SecondOrderCone34(m, s, z, λ, W, Wi)
end


mutable struct CartesianCone34
    nc::Int # number of individual cones
    m::Int # dimension of the cartesian product of all the cones
    ms::AbstractVector{Int} # dimension of each cone
    ind::AbstractVector{AbstractVector{Int}} # indices for each cone
    cones::AbstractVector{AbstractCone34}
end

function CartesianCone34(cones::AbstractVector{C}) where {C <: AbstractCone34}
    nc = length(cones)
    ms = dim.(cones)
    m = sum(ms)
    s = cumsum(ms)
    ind = [Vector(s[i] .+ (-ms[i] + 1 :0)) for i=1:nc]
    ind = [Vector(1 .+ s[i]-ms[i]:s[i]) for i=1:nc]

    return CartesianCone34(nc, m, ms, ind, cones)
end

function CartesianCone34(ms::AbstractVector{Int}, types::AbstractVector{Symbol})
    nc = length(ms)

    cones = Vector{AbstractCone34}(undef, nc)
    for i = 1:nc
        cones[i] = eval(types[i])(ms[i])
    end
    return CartesianCone34(cones)
end


mutable struct QuadraticProgram34{T}
    n::Int
    p::Int
    m::Int
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

function quadratic_program(n::Int, p::Int, m::Int; seed::Int = 100)
    Random.seed!(seed)
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
    QuadraticProgram34(n, p, m, ix, iy, iz, a, h, b, P, G, A, M)
end


# interior-point solver options
@with_kw mutable struct QuadraticProgramOptions34{T}
    iter::Int = 10
    tol::T = 1.0e-8
    # γ captures the second order term for the correction step: γ = 1.0 fully captured, γ = 0.0 not captured
    γ::T = 1.0
    # Step size selection α ≈ τ * α_max
    τ_por::T = 1.0 - 1e-8 # recommend 0.99 however it is slow
    τ_soc::T = 1.0 - 1e-2 # recommend 0.99 however it is slow
end


mutable struct Solution34{T}
    n::Int
    p::Int
    m::Int
    x::AbstractVector{T}
    y::AbstractVector{T}
    cc::CartesianCone34
end

function solution(n::Int, p::Int, ms::AbstractVector{Int}, types::AbstractVector{Symbol})
    x = ones(n)
    y = ones(p)
    cc = CartesianCone34(ms, types)
    Solution34(n, p, m, x, y, cc)
end


mutable struct Residual34{T}
    n::Int
    p::Int
    m::Int
    rx::AbstractVector{T}
    ry::AbstractVector{T}
    rcc::CartesianCone34
    # only uses rcc.cones[i].z === rz
    # only uses rcc.cones[i].s === rs
end

function residual(n::Int, p::Int, ms::AbstractVector{Int}, types::AbstractVector{Symbol})
    rx = zeros(n)
    ry = zeros(p)
    rcc = CartesianCone34(ms, types)
    Residual34(n, p, m, rx, ry, rcc)
end


mutable struct Direction34{T}
    n::Int
    p::Int
    m::Int
    Δx::AbstractVector{T}
    Δy::AbstractVector{T}
    Δcc::CartesianCone34
    # only uses Δcc.cones[i].z === Δz
    # only uses Δcc.cones[i].s === Δs
end

function direction(n::Int, p::Int, ms::AbstractVector{Int}, types::AbstractVector{Symbol})
    Δx = zeros(n)
    Δy = zeros(p)
    Δcc = CartesianCone34(ms, types)
    Direction34(n, p, m, Δx, Δy, Δcc)
end

################################################################################
# Cone Methods
################################################################################

function value(u::PORVector34)
    # if value >= 0 then u ≥ 0
    # if value >  0 then u > 0
    return minimum(u.v)
end

function value(u::SOCVector34)
    # if value >= 0 then u ≥ 0
    # if value >  0 then u > 0
    u0 = u.v[1]
    u1 = u.v[2:end]
    return (u0^2 - u1' * u1)
end

function cone_prod(u::PORVector34, w::PORVector34)
    return u.v .* w.v
end

function cone_prod(u::SOCVector34, w::SOCVector34)
    return [u.v' * w.v; u.v[1] * w.v[2:end] + w.v[1] * u.v[2:end]]
end

function cone_div(u::PORVector34, w::PORVector34)
    return w.v ./ u.v
end

function cone_div(u::SOCVector34, w::SOCVector34)
    # check CVXOPT after equation 34b for inverse-free formula.
    m = length(u.v)
    A = [u.v[1]     u.v[2:end]';
         u.v[2:end] u.v[1]*I(m - 1)]
    return A \ w.v
end

function identity_vector(c::PositiveOrthantCone34{T}) where {T}
    e = ones(T, c.m)
    return e
end

function identity_vector(c::SecondOrderCone34{T}) where {T}
    e = [ones(T, 1) ; zeros(T, c.m - 1)]
    return e
end

function identity_vector(cs::CartesianCone34)
    m = cs.m
    e = vcat(identity_vector.(cs.cones)...)
    return e
end

function scaling!(cc::CartesianCone34)
    for i in eachindex(cc.cones)
        scaling!(cc.cones[i])
    end
    return nothing
end

function scaling!(c::PositiveOrthantCone34{T}) where {T}
    m = c.m
    W = Diagonal(sqrt.(c.s) ./ sqrt.(c.z))
    Wi = Diagonal(sqrt.(c.z) ./ sqrt.(c.s))
    λ = sqrt.(c.s .* c.z)

    # @show norm(λ - Wi * s.s, Inf)
    # @show norm(λ - W * s.z, Inf)
    @assert norm(λ - Wi * c.s, Inf) < 1e-8
    @assert norm(λ - W * c.z, Inf) < 1e-8

    c.λ = λ
    c.W = W
    c.Wi = Wi
    return nothing
end

function scaling!(c::SecondOrderCone34{T}) where {T}
    m = c.m
    e = identity_vector(c)
    J = Diagonal([1; - ones(m - 1)])

    z̄ = 1 / sqrt(c.z' * J * c.z) * c.z
    s̄ = 1 / sqrt(c.s' * J * c.s) * c.s

    γ = sqrt((1 + z̄' * s̄) / 2)
    w̄ = (1 / 2γ) * (s̄ + J * z̄)
    v = 1 / sqrt(2 * (w̄[1] + 1)) * (w̄ + e)

    W̄ = 2 * v * v' - J
    W̄i = 2 * J * v * v' * J - J

    κ = exp(1/4 * log((s.s' * J * s.s) / (s.z' * J * s.z)))
    W = κ * W̄
    Wi = 1/κ * W̄i

    λbar = [γ; 1/2 * (z̄[2:end] + s̄[2:end] + (z̄[1] - s̄[1]) / (w̄[1] + 1) * w̄[2:end])]
    λ = exp(1/4 * log((c.s' * J * c.s) * (c.z' * J * c.z))) * λbar

    # @show norm(λ - Wi * c.s, Inf)
    # @show norm(λ - W  * c.z, Inf)
    # @show norm(W - W', Inf)
    # @show norm(inv(W) - Wi, Inf)

    @assert norm(W - W', Inf) < 1e-10
    # @assert norm(inv(W) - Wi, Inf) < 1e-8
    @assert norm(λ - W  * c.z, Inf) < 1e-8
    @assert norm(λ - Wi * c.s, Inf) < 1e-8

    c.λ = λ
    c.W = W
    c.Wi = Wi
    return nothing
end

function analytical_step_size(cc::CartesianCone34, Δcc::CartesianCone34)
    α = 1.0
    for i in eachindex(cc.mi)
        α = min(α, analytical_step_size(cc.cones[i], Δcc.cones[i]))
    end
    return α
end

function analytical_step_size(c::PositiveOrthantCone34{T}, Δc::PositiveOrthantCone34{T}) where {T}
    αs = por_analytical_step_size(c.λ, Δc.s)
    αz = por_analytical_step_size(c.λ, Δc.z)
    return min(αs, αz)
end

function analytical_step_size(c::SecondOrderCone34{T}, Δc::SecondOrderCone34{T}) where {T}
    αs = soc_analytical_step_size(c.λ, Δc.s)
    αz = soc_analytical_step_size(c.λ, Δc.z)
    return min(αs, αz)
end

function por_analytical_step_size(λ::AbstractVector, Δ::AbstractVector)
    m = length(λ)

    α = 1.0
    for i = 1:m
        if Δ[i] < 0
            @warn "tweak"
            # α = min(α, -(1 - 1e-10) * λ[i] / Δ[i])
            α = min(α, - λ[i] / Δ[i])
        end
    end
    # @show value(λ + α * Δ)
    return α
end

function soc_analytical_step_size(λ::AbstractVector, Δ::AbstractVector)
    m = length(λ)
    J = Diagonal([1; - ones(m - 1)])
    λbar = λ / sqrt(λ' * J * λ)
    ρ = 1 / sqrt(λ' * J * λ) * [λbar' * J * Δ; Δ[2:end] - (λbar' * J * Δ + Δ[1]) / (λbar[1] + 1) * λbar[2:end]]
    α = min(1.0, 1.0 / (norm(ρ[2:end]) - ρ[1]))
    # @show value(λ + α * Δ)
    return α
end

################################################################################
# Solver Methods
################################################################################

function residual!(r::Residual34, s::Solution34, qp::QuadraticProgram34)
    r.rx =       qp.P * s.x + qp.A' * s.y + qp.G' * s.z + qp.a
    r.ry =       qp.A * s.x                             - qp.b
    rz   = s.s + qp.G * s.x                             - qp.h

    ind = r.rcc.ind
    for i in eachindex(r.rcc.mi)
        r.rcc.cones[i].z = rz[ind[i]]
    end
    return nothing
end

function violation(r::Residual34, s::Solution34, qp::QuadraticProgram34)
    nc = r.rcc.nc
    residual!(r, s, qp)
    v1 = norm(r.rx, Inf)
    v2 = norm(r.ry, Inf)
    v3 = norm(vcat([r.rcc.cones[i].z for i=1:nc]...), Inf)
    v4 = - min([value(s.cc.cones[i].s) for i=1:nc]..., 0.0)
    v5 = - min([value(s.cc.cones[i].z) for i=1:nc]..., 0.0)

    # W, Wi, λ = scaling(s) #changed
    v6 = norm(vcat([s.cc.cones[i].s' * s.cc.cones[i].z]...), Inf) #GOOD
    # v6 = norm(abs.(cone_prod(Wi' * s.s, W * s.z)), Inf) #GOOD
    # v7 = norm(abs.(cone_prod(s.s, s.z)), Inf) #BAD
    # @show scn(v7)

    v = [v1, v2, v3, v4, v5, v6]
    vmax, imax = findmax(v)
    return norm(v, Inf)
end

function affine_direction!(Δ::Direction34, r::Residual34, s::Solution34, qp::QuadraticProgram34)
    n = qp.n
    p = qp.p
    m = qp.m

    W, Wi, λ = scaling(s) #changed
    dx = - r.rx
    dy = - r.ry
    dz = - r.rz
    ds = - cone_prod(λ, λ)

    d_ = [dx; dy; dz - W' * cone_div(λ, ds)]

    Mreg = deepcopy(qp.M)
    Mreg[qp.iz, qp.iz] += - W' * W

    Δ_ = Mreg \ d_

    Δ.Δx = Δ_[qp.ix]
    Δ.Δy = Δ_[qp.iy]
    Δ.Δz = Δ_[qp.iz]
    Δ.Δs = W' * (cone_div(λ, ds) - W * Δ.Δz)
    return nothing
end

function combined_direction!(Δ::Direction34, r::Residual34, s::Solution34,
        qp::QuadraticProgram34, σ::T; η::T = σ, γ::T = 1.0) where {T}
    n = qp.n
    p = qp.p
    m = qp.m

    W, Wi, λ = scaling(s) #changed
    μ = (λ' * λ) / m
    e = identity_vector(qp)

    dx = - (1 - η) * r.rx
    dy = - (1 - η) * r.ry
    dz = - (1 - η) * r.rz
    ds = - cone_prod(λ, λ) - γ * cone_prod(Wi' * Δ.Δs, W * Δ.Δz) + σ * μ .* e

    d_ = [dx; dy; dz - W' * cone_div(λ, ds)]

    Mreg = deepcopy(qp.M)
    Mreg[qp.iz, qp.iz] += - W' * W

    Δ_ = Mreg \ d_

    Δ.Δx = Δ_[qp.ix]
    Δ.Δy = Δ_[qp.iy]
    Δ.Δz = Δ_[qp.iz]
    Δ.Δs = W' * (cone_div(λ, ds) - W * Δ.Δz)
    return nothing
end

function centering(s::Solution34, Δ::Direction34)
    n = s.n
    p = s.p
    m = s.m
    W, Wi, λ = scaling(s) #changed

    αs = analytical_step_size(s.s, Δ.Δs)
    αz = analytical_step_size(s.z, Δ.Δz)
    α = min(αs, αz)

    ρ = 1 - α + α^2 * ((Wi' * Δ.Δs)' * (W * Δ.Δz)) / (λ' * λ) #good
    # ρ_ = ((s.s + α * Δ.Δs)' * (s.z + α * Δ.Δz)) / (s.s' * s.z) #good
    # @show norm(ρ_ - ρ)
    # @assert norm(ρ_ - ρ) < 1e-10
    σ = max(0, min(1, ρ))^3
    return α, ρ, σ
end

function update!(s::Solution34{T}, Δ::Direction34{T}; τ::T = 0.99) where {T}
    W, Wi, λ = scaling(s)

    αs = analytical_step_size(s.s, Δ.Δs / τ)
    αz = analytical_step_size(s.z, Δ.Δz / τ)
    α = min(αs, αz)

    s.s += α * Δ.Δs
    s.x += α * Δ.Δx
    s.y += α * Δ.Δy
    s.z += α * Δ.Δz
    return nothing
end

function qp_solve!(qp::QuadraticProgram34, opts::QuadraticProgramOptions34)
    n = qp.n
    p = qp.p
    m = qp.m

    s = solution(n, p, m)
    r = residual(n, p, m)
    Δ = direction(n, p, m)

    for k = 1:opts.iter
        residual!(r, s, qp)
        vio = violation(r, s, qp)
        println("iter: ", k, "   vio: ", scn(vio))
        (vio < opts.tol) && break

        affine_direction!(Δ, r, s, qp)
        α, ρ, σ = centering(s, Δ)
        combined_direction!(Δ, r, s, qp, σ, η = 0.0, γ = opts.γ)
        update!(s, Δ; τ = opts.τ_por)
    end
    return nothing
end

################################################################################
# Problem Setup
################################################################################

n = 10
m = 5
p = 3

qp = quadratic_program(n, p, m)
s = solution(n, p, m)
r = residual(n, p, m)
Δ = direction(n, p, m)

residual!(r, s, qp)
violation(r, s, qp)
scaling(s)
affine_direction!(Δ, r, s, qp)
α, ρ, σ = centering(s, Δ)
combined_direction!(Δ, r, s, qp, σ)
update!(s, Δ)

n = 103
m = 53
p = 34
opts = QuadraticProgramOptions34(iter = 100, τ_por = 0.99, tol = 1e-8)
qp = quadratic_program(n, p, m, seed = 200)
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
cs = CartesianCone34(ms, types)
@test typeof(cs.cones[1]) <: PositiveOrthantCone34
@test typeof(cs.cones[2]) <: SecondOrderCone34

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
