using Random
using LinearAlgebra

################################################################################
# Structures
################################################################################

# Filter
abstract type Filter{T}
end

mutable struct Filter12{T} <: Filter{T}
    nx::Int
    nu::Int
    nz::Int
    A::AbstractMatrix{T}
    B::AbstractMatrix{T}
    C::AbstractMatrix{T}
    Q::AbstractMatrix{T}
    R::AbstractMatrix{T}
end


################################################################################
# Methods
################################################################################

function process_model(filter::Filter{T}, x0::AbstractVector{T}, u0::AbstractVector{T}) where {T}
    nx = filter.nx
    A = filter.A
    B = filter.B
    Q = filter.Q
    ϵ = randn(nx)
    x1 = A * x0 + B * u0 + sqrt(Q) * rand(nx)
    return x1
end

function measurement_model(filter::Filter{T}, x0::AbstractVector{T}) where {T}
    nz = filter.nz
    C = filter.C
    R = filter.R
    z0 = C * x0 + sqrt(R) * rand(nz)
    return z0
end

################################################################################
# Script
################################################################################

Random.seed!(100)

# Dimensions
nx = 4 # state
nu = 2 # control
nz = 2 # observation

# Process
Δt = 0.1
A = Diagonal(ones(nx)) + Δt * [0 0 Δt 0;
                               0 0 0 Δt;
                               0 0 0  0;
                               0 0 0  0;]
B = [0  0;
     0  0;
     Δt 0;
     0 Δt;]

Q̄ = 0.1*rand(nx, nx)
Q = Q̄ * Q̄'

# Measurement
C = [1 0 0 0;
     0 1 0 0.;]
R̄ = 0.1*rand(nz, nz)
R = R̄ * R̄'

# Filter
kf = Filter12(nx, nu, nz, A, B, C, Q, R)

function step!(filter::Filter{T}, x0, u0, z1, Σ0)
    A = filter.A
    xbar = process_model(filter, x0, u0)
    Σbar = A * Σ0 * A' + Q
    K = Σbar * C' * inv(C * Σbar * C' + Qt)

    K =

    zbar = measurement_model(filter, x1)

    return x1, Σ1
end

x0 = rand(nx)
u0 = rand(nu)
x1 = process_model(kf, x0, u0)
z1 = measurement_model(kf, x1)





#
# plt = plot()
# for i = 1:400
#     x = sqrt(R) * randn(2)
#     scatter!(plt, aspectratio = 1.0, legend = false, [x[1]], [x[2]])
# end
# display(plt)
#
