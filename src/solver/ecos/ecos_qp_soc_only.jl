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
    s0 = [1; 1e-2 *rand(m - 1)]
    x0 = ones(n)
    y0 = ones(p)
    z0 = [1; 1e-2 *rand(m - 1)]
    @assert strictly_positive(s0) #changed
    @assert strictly_positive(z0) #changed
    Solution15(n, p, m, s0, x0, y0, z0)
end

function strictly_positive(u::AbstractVector)
    return value(u) > 0.0
end

function value(u::AbstractVector)
    u0 = u[1]
    u1 = u[2:end]
    return (u0 - u1'*u1)
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
    v4 = - min(value(s.s), 0.0) #changed
    v5 = - min(value(s.z), 0.0) #changed
    v6 = norm(abs.(s.s .* s.z), Inf) #kept see Eq.(8) and above
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

function combined_direction!(Δ::Direction15, r::Residual15, s::Solution15,
        qp::QuadraticProgram15, σ::T; η::T = σ, γ::T = 1.0) where {T}
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

function scaling(s::Solution15)
    m = s.m
    e = [1; zeros(m - 1)]
    J = Diagonal([1; - ones(m - 1)])

    z̄ = 1 / sqrt(s.z' * J * s.z) * s.z
    s̄ = 1 / sqrt(s.s' * J * s.s) * s.s

    γ = sqrt((1 + z̄' * s̄) / 2)
    w̄ = 1 / (2γ) * (s̄ + J * z̄)
    v = 1 / sqrt(2 * (w̄[1] + 1)) * (w̄ + e)

    W̄ = 2 * v * v' - J
    W̄i = 2 * J * v * v' * J - J

    sca = exp(0.25 * log((s.s' * J * s.s) / (s.z' * J * s.z)))
    W = sca * W̄ #changed
    Wi = 1/sca * W̄i #changed

    λbar = [γ; 0.5 * (z̄[2:end] + s̄[2:end] + (z̄[1] - s̄[1]) / (w̄[1] + 1) * w̄[2:end])]
    λ = exp(0.25 * log((s.s' * J * s.s) * (s.z' * J * s.z))) * λbar

    # @show norm(λ - Wi * s.s, Inf)
    # @show norm(λ - W  * s.z, Inf)
    # @show norm(W - W', Inf)
    # @show norm(inv(W) - Wi, Inf)

    @assert norm(W - W', Inf) < 1e-10
    @assert norm(inv(W) - Wi, Inf) < 1e-8
    @assert norm(λ - W  * s.z, Inf) < 1e-10
    @assert norm(λ - Wi * s.s, Inf) < 1e-10
    return W, Wi, λ
end

function cone_prod(u, v)
    return [u' * v; u[1] * v[2:end] + v[1] * u[2:end]]
end

function cone_div(u, v)
    # check CVXOPT after equation 22b for inverse-free formula.
    m = length(u)
    A = [u[1]     u[2:end]';
         u[2:end] u[1]*I(m-1)]
    return A \ v
end

function step_size_centering(s::Solution15, Δ::Direction15)
    n = s.n
    p = s.p
    m = s.m
    W, Wi, λ = scaling(s) #changed

    # α = 1.0
    # for i = 1:m
    #     if Δ.Δs[i] < 0
    #         α = min(α, -(1 - 1e-10) * s.s[i] / Δ.Δs[i])
    #     end
    #     if Δ.Δz[i] < 0
    #         α = min(α, -(1 - 1e-10) * s.z[i] / Δ.Δz[i])
    #     end
    # end

    # J = Diagonal([1; - ones(m - 1)])
    # Δs̃ = Wi' * Δ.Δs
    # Δz̃ = W  *  Δ.Δz
    # λbar = λ ./ sqrt(λ' * J * λ)
    # ρα = 1 / sqrt(λ' * J * λ) * [λbar' * J * Δs̃; Δs̃[2:end] - (λbar' * J * Δs̃ + Δs̃[1]) / (λbar[1] + 1) .* λbar[2:end]]
    # σα = 1 / sqrt(λ' * J * λ) * [λbar' * J * Δz̃; Δz̃[2:end] - (λbar' * J * Δz̃ + Δz̃[1]) / (λbar[1] + 1) .* λbar[2:end]]
    #
    # α = 1 / max(0.0, norm(ρα[2:end]) - ρα[1], norm(σα[2:end]) - σα[1])
    # @show α
    # @show value(s.s + α * Δ.Δs)
    # @show value(s.z + α * Δ.Δz)
    # @assert strictly_positive(s.s + α * Δ.Δs)
    # @assert strictly_positive(s.z + α * Δ.Δz)

    @warn "changed step size"
    α = step_size(s, Δ, τ = 1.0)

    ρ_ = ((s.s + α * Δ.Δs)' * (s.z + α * Δ.Δz)) / (s.s' * s.z)
    ρ = 1 - α + α^2 * ((Wi' * Δ.Δs)' * (W * Δ.Δz)) / (λ' * λ) #changed
    @assert norm(ρ_ - ρ) < 1e-10
    σ = max(0, min(1, ρ))^3
    return α, ρ, σ
end

function step_size(s::Solution15, Δ::Direction15; τ::T = 1.0) where {T}
    α = 1.0
    val = min(value(s.s + α * (Δ.Δs ./ τ)),
              value(s.z + α * (Δ.Δz ./ τ)))
    while val < 0.0
        α *= 0.9
        val = min(value(s.s + α * (Δ.Δs ./ τ)),
                  value(s.z + α * (Δ.Δz ./ τ)))
    end
    # @show α
    # @show value(s.s + α * Δ.Δs)
    # @show value(s.z + α * Δ.Δz)
    @assert strictly_positive(s.s + α * Δ.Δs)
    @assert strictly_positive(s.z + α * Δ.Δz)
    return α
end

function identity_vector(qp::QuadraticProgram15)
    m = qp.m
    e = [1; zeros(m - 1)]
    return e
end

function update!(s::Solution15{T}, Δ::Direction15{T}; τ::T = 0.99) where {T}
    n = s.n
    p = s.p
    m = s.m

    W, Wi, λ = scaling(s) #changed
    # α = 1.0
    # for i = 1:m
    #     if (Wi' * Δ.Δs)[i] < 0
    #         α = min(α, - τ * λ[i] / (Wi' * Δ.Δs)[i])
    #     end
    #     if (W * Δ.Δz)[i] < 0
    #         α = min(α, - τ * λ[i] / (W * Δ.Δz)[i])
    #     end
    # end

    # J = Diagonal([1; - ones(m - 1)])
    # Δs̃ = (Wi' * Δ.Δs) ./ τ #TAU
    # Δz̃ = (W  *  Δ.Δz) ./ τ #TAU
    # λbar = λ ./ sqrt(λ' * J * λ)
    # ρα = 1 / sqrt(λ' * J * λ) * [λbar' * J * Δs̃; Δs̃[2:end] - (λbar' * J * Δs̃ + Δs̃[1]) / (λbar[1] + 1) .* λbar[2:end]]
    # σα = 1 / sqrt(λ' * J * λ) * [λbar' * J * Δz̃; Δz̃[2:end] - (λbar' * J * Δz̃ + Δz̃[1]) / (λbar[1] + 1) .* λbar[2:end]]
    #
    # α = 1 / max(0.0, norm(ρα[2:end]) - ρα[1], norm(σα[2:end]) - σα[1])
    # @show α
    # @show value(s.s + α * Δ.Δs)
    # @show value(s.z + α * Δ.Δz)
    # @assert strictly_positive(s.s + α * Δ.Δs)
    # @assert strictly_positive(s.z + α * Δ.Δz)

    @warn "changed step size"
    α = step_size(s, Δ, τ = τ)

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
# Test
################################################################################

function test_division(m::Int, seed::Int = 100)
    Random.seed!(100)
    u = [1; 1e-2 * rand(m - 1)]
    v = [1; 1e-2 * rand(m - 1)]
    @show norm(cone_prod(u, cone_div(u, v)) - v)
    @test norm(cone_prod(u, cone_div(u, v)) - v) < 1e-10
    return nothing
end

test_division(4)


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
