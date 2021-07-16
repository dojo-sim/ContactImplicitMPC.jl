################################################################################
# Linear Program
################################################################################
using Random
using LinearAlgebra

"""
https://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf

min   0.5 * x' * P * x + c' * x
s.t.  G * x + s = h
      A * x = b
      s ≧ 0

Primal variables:
s ∈ [m]
x ∈ [n]
y ∈ [p]
z ∈ [m]

Problem data:
c ∈ [n]
h ∈ [m]
b ∈ [p]
P ∈ [n,n]
G ∈ [m,n]
A ∈ [p,n]

"""

################################################################################
# Methods
################################################################################

mutable struct QuadraticProgram15{T}
    n::Int
    p::Int
    m::Int
    ix::AbstractVector{Int}
    iy::AbstractVector{Int}
    iz::AbstractVector{Int}
    c::AbstractVector{T}
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

    c = rand(n)
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
    QuadraticProgram15(n, p, m, ix, iy, iz, c, h, b, P, G, A, M)
end

# interior-point solver options
@with_kw mutable struct QuadraticProgramOptions15{T}
    tol::T = 1.0e-8
    γ::T = 1.0
    τ::T = 1.0 - 1e-8 # recommend 0.99 this is slow
end

mutable struct Solution15{T}
    n::Int
    p::Int
    m::Int
    s::AbstractVector{T}
    x::AbstractVector{T}
    y::AbstractVector{T}
    z::AbstractVector{T}
end

function solution(n::Int, p::Int, m::Int;)
    s0 = ones(m)
    x0 = ones(n)
    y0 = ones(p)
    z0 = ones(m)
    @assert all(s0 .> 0.0)
    @assert all(z0 .> 0.0)
    Solution15(n, p, m, s0, x0, y0, z0)
end

mutable struct Residual15{T}
    n::Int
    p::Int
    m::Int
    rx::AbstractVector{T}
    ry::AbstractVector{T}
    rz::AbstractVector{T}
end

function residual(n::Int, p::Int, m::Int;)
    rx0 = zeros(n)
    ry0 = zeros(p)
    rz0 = zeros(m)
    Residual15(n, p, m, rx0, ry0, rz0)
end

function residual!(r::Residual15, s::Solution15, qp::QuadraticProgram15)
    r.rx =       qp.P * s.x + qp.A' * s.y + qp.G' * s.z + qp.c
    r.ry =       qp.A * s.x                             - qp.b
    r.rz = s.s + qp.G * s.x                             - qp.h

    r_ = [zeros(n); zeros(p); s.s] + qp.M * [s.x; s.y; s.z] + [qp.c; - qp.b; - qp.h]
    @assert norm(r_[qp.ix] - r.rx, Inf) < 1e-10
    @assert norm(r_[qp.iy] - r.ry, Inf) < 1e-10
    @assert norm(r_[qp.iz] - r.rz, Inf) < 1e-10
    return nothing
end

function violation(r::Residual15, s::Solution15, qp::QuadraticProgram15)
    residual!(r, s, qp)
    v1 = norm(r.rx, Inf)
    v2 = norm(r.ry, Inf)
    v3 = norm(r.rz, Inf)
    v4 = norm(- min.(s.s, 0.0), Inf)
    v5 = norm(- min.(s.z, 0.0), Inf)
    v6 = norm(abs.(s.s .* s.z), Inf)
    v = [v1, v2, v3, v4, v5, v6]
    vmax, imax = findmax(v)
    return norm(v, Inf)
end

mutable struct Direction15{T}
    n::Int
    p::Int
    m::Int
    Δx::AbstractVector{T}
    Δy::AbstractVector{T}
    Δz::AbstractVector{T}
    Δs::AbstractVector{T}
end

function direction(n::Int, p::Int, m::Int;)
    Δx0 = zeros(n)
    Δy0 = zeros(p)
    Δz0 = zeros(m)
    Δs0 = zeros(m)
    Direction15(n, p, m, Δx0, Δy0, Δz0, Δx0)
end

function affine_direction!(Δ::Direction15, r::Residual15, s::Solution15, qp::QuadraticProgram15)
    n = qp.n
    p = qp.p
    m = qp.m

    W, λ = scaling(s)
    dx = - r.rx
    dy = - r.ry
    dz = - r.rz
    ds = - s.s .* s.z
    dsS = - cone_prod(λ, λ)
    # @show norm(dsS - ds)
    @assert norm(dsS - ds) < 1e-10

    d_ = [dx; dy; dz + s.s]
    d_S = [dx; dy; dz - W' * cone_div(λ, ds)]
    # @show norm(d_ - d_S)
    @assert norm(d_ - d_S) < 1e-10

    Mreg = deepcopy(qp.M)
    Mreg[qp.iz, qp.iz] += - W' * W
    # Mreg[qp.iz, qp.iz] += - Diagonal(s.s ./ s.z) # works too
    # @show norm((W' * W) - Diagonal(s.s ./ s.z), Inf)
    # @assert norm((W' * W) - Diagonal(s.s ./ s.z), Inf) < 1e-10

    Δ_ = Mreg \ d_

    Δ.Δx = Δ_[qp.ix]
    Δ.Δy = Δ_[qp.iy]
    Δ.Δz = Δ_[qp.iz]

    Δs = W' * (- λ - W * Δ.Δz)
    ΔsS = W' * (cone_div(λ, ds) - W * Δ.Δz)
    Δ.Δs = Δs

    # @show norm(ΔsS - Δs, Inf)
    @assert norm(ΔsS - Δs, Inf) < 1e-10
    return nothing
end

function combined_direction!(Δ::Direction15, r::Residual15, s::Solution15,
        qp::QuadraticProgram15, σ::T; η::T = σ, γ::T = 1.0) where {T}
    n = qp.n
    p = qp.p
    m = qp.m

    W, λ = scaling(s)
    μS = (λ' * λ) / m
    μ = (s.s' * s.z) / m
    # @show norm(μ - μS, Inf)
    @assert norm(μ - μS, Inf) < 1e-10

    e = identity_vector(qp)

    dx = - (1 - η) * r.rx
    dy = - (1 - η) * r.ry
    dz = - (1 - η) * r.rz
    ds = - s.s .* s.z - γ * cone_prod(inv(W)' * Δ.Δs, W * Δ.Δz) + σ * μ .* e
    dsS = - cone_prod(λ, λ) - γ * cone_prod(inv(W)' * Δ.Δs, W * Δ.Δz) + σ * μ .* e
    # @show norm(dsS - ds)
    @assert norm(dsS - ds) < 1e-10

    d_ = [dx; dy; dz - ds ./ s.z]
    d_S = [dx; dy; dz - W' * cone_div(λ, ds)]
    # @show norm(d_ - d_S)
    @assert norm(d_ - d_S) < 1e-10

    Mreg = deepcopy(qp.M)
    Mreg[qp.iz, qp.iz] += - W' * W

    Δ_ = Mreg \ d_

    Δ.Δx = Δ_[qp.ix]
    Δ.Δy = Δ_[qp.iy]
    Δ.Δz = Δ_[qp.iz]
    Δ.Δs = W' * (cone_div(λ, ds) - W * Δ.Δz)
    return nothing
end

function scaling(s::Solution15)
    W = Diagonal(sqrt.(s.s) ./ sqrt.(s.z))
    λ = sqrt.(s.s .* s.z)
    # @show norm((W' * W) - Diagonal(s.s ./ s.z), Inf)
    # @assert norm((W' * W) - Diagonal(s.s ./ s.z), Inf) < 1e-10

    # @show norm(λ - inv(W) * s.s, Inf)
    # @show norm(λ - W * s.z, Inf)
    @assert norm(λ - inv(W) * s.s, Inf) < 1e-10
    @assert norm(λ - W * s.z, Inf) < 1e-10
    return W, λ
end

function cone_prod(u, v)
    # @warn "not implemented"
    return u .* v
end

function cone_div(u, v)
    # @warn "not implemented"
    return v ./ u
end

function step_size_centering(s::Solution15, Δ::Direction15)
    n = s.n
    p = s.p
    m = s.m
    W, λ = scaling(s)

    α = 1.0
    for i = 1:m
        if Δ.Δs[i] < 0
            α = min(α, -(1 - 1e-10) * s.s[i] / Δ.Δs[i])
        end
        if Δ.Δz[i] < 0
            α = min(α, -(1 - 1e-10) * s.z[i] / Δ.Δz[i])
        end
    end
    @assert minimum(s.s + α * Δ.Δs) >= -10.0
    @assert minimum(s.z + α * Δ.Δz) >= -10.0

    ρ_ = ((s.s + α * Δ.Δs)' * (s.z + α * Δ.Δz)) / (s.s' * s.z)
    ρ = 1 - α + α^2 * ((inv(W)' * Δ.Δs)' * (W * Δ.Δz)) / (λ' * λ)
    @assert norm(ρ_ - ρ) < 1e-10
    σ = max(0, min(1, ρ))^3
    return α, ρ, σ
end

function identity_vector(qp::QuadraticProgram15)
    m = qp.m
    e = ones(m)
    return e
end

function update!(s::Solution15{T}, Δ::Direction15{T}; τ::T = 0.99) where {T}
    n = s.n
    p = s.p
    m = s.m

    W, λ = scaling(s)
    α = 1.0
    for i = 1:m
        if (inv(W)' * Δ.Δs)[i] < 0
            α = min(α, - τ * λ[i] / (inv(W)' * Δ.Δs)[i])
        end
        if (W * Δ.Δz)[i] < 0
            α = min(α, - τ * λ[i] / (W * Δ.Δz)[i])
        end
    end

    α_ = 1.0
    for i = 1:m
        if Δ.Δs[i] < 0
            α_ = min(α_, -τ * s.s[i] / Δ.Δs[i])
        end
        if Δ.Δz[i] < 0
            α_ = min(α_, -τ * s.z[i] / Δ.Δz[i])
        end
    end
    # @show norm(α - α_, Inf)
    @assert norm(α - α_, Inf) < 1e-10
    s.s += α * Δ.Δs
    s.x += α * Δ.Δx
    s.y += α * Δ.Δy
    s.z += α * Δ.Δz
    return nothing
end

function qp_solve!(qp::QuadraticProgram15, opts::QuadraticProgramOptions15)
    n = qp.n
    p = qp.p
    m = qp.m

    s = solution(n, p, m)
    r = residual(n, p, m)
    Δ = direction(n, p, m)

    for k = 1:10
        residual!(r, s, qp)
        vio = violation(r, s, qp)
        println("iter: ", k, "   vio: ", scn(vio))
        (vio < opts.tol) && break

        affine_direction!(Δ, r, s, qp)
        α, ρ, σ = step_size_centering(s, Δ)
        combined_direction!(Δ, r, s, qp, σ, η = 0.0, γ = opts.γ)
        update!(s, Δ; τ = opts.τ)

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
opts = QuadraticProgramOptions15()

s = solution(n, p, m)
r = residual(n, p, m)
Δ = direction(n, p, m)

residual!(r, s, qp)
violation(r, s, qp)
scaling(s)
affine_direction!(Δ, r, s, qp)
α, ρ, σ = step_size_centering(s, Δ)
combined_direction!(Δ, r, s, qp, σ)
update!(s, Δ)

qp_solve!(qp, opts)

s.s = exp.((rand(10m) .* 10 .- 5) .* log(10))
s.z = exp.((rand(10m) .* 10 .- 5) .* log(10))

scaling(s)
