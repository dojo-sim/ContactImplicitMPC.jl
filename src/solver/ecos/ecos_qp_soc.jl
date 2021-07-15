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

mutable struct QuadraticProgram16{T}
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
    QuadraticProgram16(n, p, m, ix, iy, iz, c, h, b, P, G, A, M)
end

# interior-point solver options
@with_kw mutable struct QuadraticProgramOptions16{T}
    iter::Int = 10
    tol::T = 1.0e-8
    γ::T = 1.0
    τ::T = 1.0 - 1e-8 # recommend 0.99 however it is slow
end

mutable struct Solution16{T}
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
    Solution16(n, p, m, s0, x0, y0, z0)
end

function strictly_positive(u::AbstractVector)
    return value(u) > 0.0
end

function positive(u::AbstractVector)
    return value(u) >= 0.0
end

function value(u::AbstractVector)
    u0 = u[1]
    u1 = u[2:end]
    return (u0^2 - u1' * u1)
end

mutable struct Residual16{T}
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
    Residual16(n, p, m, rx0, ry0, rz0)
end

function residual!(r::Residual16, s::Solution16, qp::QuadraticProgram16)
    r.rx =       qp.P * s.x + qp.A' * s.y + qp.G' * s.z + qp.c
    r.ry =       qp.A * s.x                             - qp.b
    r.rz = s.s + qp.G * s.x                             - qp.h

    r_ = [zeros(n); zeros(p); s.s] + qp.M * [s.x; s.y; s.z] + [qp.c; - qp.b; - qp.h]
    @assert norm(r_[qp.ix] - r.rx, Inf) < 1e-10
    @assert norm(r_[qp.iy] - r.ry, Inf) < 1e-10
    @assert norm(r_[qp.iz] - r.rz, Inf) < 1e-10
    return nothing
end

function violation(r::Residual16, s::Solution16, qp::QuadraticProgram16)
    residual!(r, s, qp)
    v1 = norm(r.rx, Inf)
    v2 = norm(r.ry, Inf)
    v3 = norm(r.rz, Inf)
    v4 = - min(value(s.s), 0.0) #changed
    v5 = - min(value(s.z), 0.0) #changed
    W, Wi, λ = scaling(s) #changed
    v6 = norm(abs.(cone_prod(Wi' * s.s, W * s.z)), Inf) #changed
    # v6 = norm(abs.(cone_prod(s.s, s.z)), Inf) #changed
    v = [v1, v2, v3, v4, v5, v6]
    vmax, imax = findmax(v)
    @show scn(value(s.s))
    @show scn(value(s.z))
    @show imax
    @show cone_prod(s.s, s.z)
    @show scn.(v)
    return norm(v, Inf)
end

mutable struct Direction16{T}
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
    Direction16(n, p, m, Δx0, Δy0, Δz0, Δs0)
end

function affine_direction!(Δ::Direction16, r::Residual16, s::Solution16, qp::QuadraticProgram16)
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

function combined_direction!(Δ::Direction16, r::Residual16, s::Solution16,
        qp::QuadraticProgram16, σ::T; η::T = σ, γ::T = 1.0) where {T}
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

function scaling(s::Solution16)
    m = s.m
    e = [1; zeros(m - 1)]
    J = Diagonal([1; - ones(m - 1)])

    z̄ = 1 / sqrt(s.z' * J * s.z) * s.z
    s̄ = 1 / sqrt(s.s' * J * s.s) * s.s

    γ = sqrt((1 + z̄' * s̄) / 2)
    w̄ = (1 / 2γ) * (s̄ + J * z̄)
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
    # @assert norm(inv(W) - Wi, Inf) < 1e-8
    @assert norm(λ - W  * s.z, Inf) < 1e-8
    @assert norm(λ - Wi * s.s, Inf) < 1e-8
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

function step_size_centering(s::Solution16, Δ::Direction16)
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
    # ϵ = 1e-7
    # # α0 = 1 / (max(0.0, norm(ρα[2:end]) - ρα[1], norm(σα[2:end]) - σα[1]) + ϵ)
    # α0 = min(1 / (norm(ρα[2:end]) - ρα[1]), 1 / (norm(σα[2:end]) - σα[1]))
    # α0 = min(α0, 1.0 - ϵ)
    # # @show α0
    # # @show value(s.s + α0 * Δ.Δs)
    # # @show value(s.z + α0 * Δ.Δz)

    Δs̃ = (Wi' * Δ.Δs) # ./ τ #TAU
    Δz̃ = (W  *  Δ.Δz) # ./ τ #TAU
    α0 = analytical_step_size(λ, Δs̃, Δz̃)
    @assert value(s.s + α0 * Δ.Δs) > -1e-10
    @assert value(s.z + α0 * Δ.Δz) > -1e-10

    @warn "changed step size"
    # α = empirical_step_size(s, Δ, τ = 1.0)
    αs = analytical_step_size(s.s, Δ.Δs)
    αz = analytical_step_size(s.z, Δ.Δz)
    α = min(αs, αz)

    # @show scn(α0, digits = 2)
    # @show scn(α, digits = 2)

    ρ_ = ((s.s + α * Δ.Δs)' * (s.z + α * Δ.Δz)) / (s.s' * s.z)
    ρ = 1 - α + α^2 * ((Wi' * Δ.Δs)' * (W * Δ.Δz)) / (λ' * λ) #changed
    # @show norm(ρ_ - ρ)
    # @assert norm(ρ_ - ρ) < 1e-10
    σ = max(0, min(1, ρ))^3
    return α, ρ, σ
end

# function step_size(s::Solution16, Δ::Direction16; τ::T = 1.0) where {T}
#     α = 1.0
#     val = min(value(s.s + α * (Δ.Δs ./ τ)),
#               value(s.z + α * (Δ.Δz ./ τ)))
#     while val <= 1e-10
#         α *= 0.99
#         val = min(value(s.s + α * (Δ.Δs ./ τ)),
#                   value(s.z + α * (Δ.Δz ./ τ)))
#     end
#     # @show α
#     # @show value(s.s + α * Δ.Δs)
#     # @show value(s.z + α * Δ.Δz)
#     @assert strictly_positive(s.s + α * Δ.Δs)
#     @assert strictly_positive(s.z + α * Δ.Δz)
#     return α
# end

function identity_vector(qp::QuadraticProgram16)
    m = qp.m
    e = [1; zeros(m - 1)]
    return e
end

function update!(s::Solution16{T}, Δ::Direction16{T}; τ::T = 0.99) where {T}
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
    #
    # J = Diagonal([1; - ones(m - 1)])
    # Δs̃ = (Wi' * Δ.Δs) # ./ τ #TAU
    # Δz̃ = (W  *  Δ.Δz) # ./ τ #TAU
    # λbar = λ ./ sqrt(λ' * J * λ)
    # ρα = 1 / sqrt(λ' * J * λ) * [λbar' * J * Δs̃; Δs̃[2:end] - ((λbar' * J * Δs̃ + Δs̃[1]) / (λbar[1] + 1)) .* λbar[2:end]]
    # σα = 1 / sqrt(λ' * J * λ) * [λbar' * J * Δz̃; Δz̃[2:end] - ((λbar' * J * Δz̃ + Δz̃[1]) / (λbar[1] + 1)) .* λbar[2:end]]
    #
    # ϵ = 1e-7
    # α0 = min(1 / (norm(ρα[2:end]) - ρα[1]), 1 / (norm(σα[2:end]) - σα[1]))
    # α0 = max(0.0, norm(ρα[2:end]) - ρα[1], norm(σα[2:end]) - σα[1]) + ϵ
    # @show norm(ρα[2:end]) - ρα[1]
    # @show norm(σα[2:end]) - σα[1]
    # @show scn(max(0.0, norm(ρα[2:end]) - ρα[1], norm(σα[2:end]) - σα[1]) + ϵ, digits=2)
    # @show 1 / (max(0.0, norm(ρα[2:end]) - ρα[1], norm(σα[2:end]) - σα[1]) + ϵ)
    # α0 = min(α0, 1.0 - ϵ)
    # @show α0
    # @show value(s.s + α0 * Δ.Δs)
    # @show value(s.z + α0 * Δ.Δz)

    Δs̃ = (Wi' * Δ.Δs) # ./ τ #TAU
    Δz̃ = (W  *  Δ.Δz) # ./ τ #TAU
    α0 = analytical_step_size(λ, Δs̃, Δz̃)
    @assert value(s.s + α0 * Δ.Δs) > - 1e-10
    @assert value(s.z + α0 * Δ.Δz) > - 1e-10

    # @warn "changed step size"
    # α = empirical_step_size(s, Δ, τ = τ)
    αs = analytical_step_size(s.s, Δ.Δs)
    αz = analytical_step_size(s.z, Δ.Δz)
    α = min(αs*0.99, αz*0.99)

    # println("α0: ", scn(α0, digits = 2))
    # println("α:  ", scn(α, digits = 2))

    s.s += α * Δ.Δs
    s.x += α * Δ.Δx
    s.y += α * Δ.Δy
    s.z += α * Δ.Δz
    return nothing
end

function qp_solve!(qp::QuadraticProgram16, opts::QuadraticProgramOptions16)
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
    # @show norm(cone_prod(u, cone_div(u, v)) - v)
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

opts = QuadraticProgramOptions16(iter = 22, τ = 1.0)
qp = quadratic_program(n, p, m, seed = 200)
qp_solve!(qp, opts)

a = 10
a = 10
a = 10
a = 10


#
#
#
#
# e = identity_vector(qp)
# v = 10*ones(qp.m)
# α = 0.1
# value(e + α*v)
#
# α = norm(v[2:end]) - v[1]
#
# function test_step_size(qp::QuadraticProgram16)
#     e = identity_vector(qp)
#     v = 100 * (rand(qp.m) .- 0.5)
#     α = 1 / (norm(v[2:end]) - v[1])
#     @show value(e + α*v)
#     @test norm(value(e + α*v)) < 1e-10
#     return nothing
# end
#
# e = identity_vector(qp)
# v = 10*ones(qp.m)
# # α = 0.1
# α = 1/ (norm(v[2:end]) - v[1])
# value(e + α*v)
#
#
# for t = 1:100
#     test_step_size(qp)
# end
#




function empirical_step_size(s::Solution16, Δ::Direction16; τ::T = 1.0) where {T}
    α = 1.0
    val = min(value(s.s + α * (Δ.Δs ./ τ)),
              value(s.z + α * (Δ.Δz ./ τ)))
    while val <= 1e-10
        α *= 0.98
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

function analytical_step_size(λ::AbstractVector, Δs̃::AbstractVector, Δz̃::AbstractVector)
    m = length(λ)

    J = Diagonal([1; - ones(m - 1)])
    λbar = λ / sqrt(λ' * J * λ)

    ρ = 1 / sqrt(λ' * J * λ) * [λbar' * J * Δs̃; Δs̃[2:end] - (λbar' * J * Δs̃ + Δs̃[1]) / (λbar[1] + 1) * λbar[2:end]]
    σ = 1 / sqrt(λ' * J * λ) * [λbar' * J * Δz̃; Δz̃[2:end] - (λbar' * J * Δz̃ + Δz̃[1]) / (λbar[1] + 1) * λbar[2:end]]

    α = min(1.0, 1.0 / (norm(ρ[2:end]) - ρ[1]), 1.0 / (norm(σ[2:end]) - σ[1]))
    # @show value(λ + α * Δs̃)
    # @show value(λ + α * Δz̃)
    return α
end

function analytical_step_size(λ::AbstractVector, Δ::AbstractVector)
    m = length(λ)

    J = Diagonal([1; - ones(m - 1)])
    λbar = λ / sqrt(λ' * J * λ)

    ρ = 1 / sqrt(λ' * J * λ) * [λbar' * J * Δ; Δ[2:end] - (λbar' * J * Δ + Δ[1]) / (λbar[1] + 1) * λbar[2:end]]

    α = min(1.0, 1.0 / (norm(ρ[2:end]) - ρ[1]))
    @show value(λ + α * Δ)
    return α
end



Random.seed!(100)
λ = [3.0; rand(m - 1)]
Δs̃ = 100*rand(m)
Δz̃ = 100*rand(m)

s = solution(n, p, m)
Δ = direction(n, p, m)
s.s = λ
s.z = λ
Δ.Δs = Δs̃
Δ.Δz = Δz̃

αE = step_size(s, Δ)
αA = analytical_step_size(λ, Δs̃, Δz̃)

value(λ + αE*Δs̃)
value(λ + αE*Δz̃)

value(λ + αA*Δs̃)
value(λ + αA*Δz̃)

Random.seed!(10)
λ = [3.0; rand(m - 1)]
Δ = 100*rand(m)
analytical_step_size(λ, Δ)
