include("interior_point.jl")
include("simulator.jl")
using InteractiveUtils, BenchmarkTools, ModelingToolkit

"""
    particle dynamics
    - 3D particle subject to contact forces

    - configuration: q = (x, y, z) ∈ R³
    - impact force (magnitude): γ ∈ R₊
    - friction force: β ∈ R⁴₊
        - friction coefficient: μ ∈ R₊

    Discrete Mechanics and Variational Integrators
        pg. 363
"""

struct Particle
      m # mass
      g # gravity
      μ # friction coefficient

      h # time step

      nq # configuration dimension
      nu # control dimesion
      nc # contact dimension
      nb # friction dimension
      nw # disturbance dimension
end

# mass matrix
function mass_matrix(model, q)
    m = model.m

    Diagonal(@SVector [m, m, m])
end

# gravity
function gravity(model, q)
    m = model.m
    g = model.g

    @SVector [0.0, 0.0, m * g]
end

# signed distance function
function signed_distance(model, q)
    q[3:3]
end

# control Jacobian
function control_jacobian(model, q)
    SMatrix{3, 3}([1.0 0.0 0.0;
                   0.0 1.0 0.0;
                   0.0 0.0 1.0])
end

# normal Jacobian
function normal_jacobian(model, q)
    SMatrix{1, 3}([0.0, 0.0, 1.0])
end

# tangent Jacobian
function tangent_jacobian(model, q)
    SMatrix{4, 3}([1.0 0.0 0.0;
                   0.0 1.0 0.0;
                   -1.0 0.0 0.0;
                   0.0 -1.0 0.0])
end

# dynamics
function dynamics(model, q0, q1, q2, u, γ, b, w, h)
    nq = model.nq
    SVector{nq}(mass_matrix(model, q0) * (q1 - q0) ./ h
              - mass_matrix(model, q1) * (q2 - q1) ./ h
              - h * gravity(model, q1)
              + control_jacobian(model, q2) * u
              + transpose(normal_jacobian(model, q2)) * γ
              + transpose(tangent_jacobian(model, q2)) * b)
end

function num_var(model)
    model.nq + model.nc + model.nb + model.nc + model.nb + 2 * model.nc
end

function num_data(model)
    model.nq + model.nq + model.nu + model.nw + 1
end
num_data(model)

function z_initialize!(z, model, q1)
    z .= 1.0
    z[1:3] = q1
end

function θ_initialize!(θ, model, q0, q1, u, w, h)
    θ[1:3] = q0
    θ[4:6] = q1
    θ[7:9] = u
    θ[10:10] = w
    θ[11:11] = h
end

function unpack_z(z, model)
    q2 = view(z, 1:3)  # configuration
    γ = view(z, 4:4)   # normal impulse
    b = view(z, 5:8)   # friction impulse (double parameterized)
    ψ = view(z, 9:9)   # friction cone dual (magnitude of tangential velocity)
    η = view(z, 10:13) # friction impulse duals
    s1 = view(z, 14:14)# signed-distance slack
    s2 = view(z, 15:15)# friction cone slack

    return q2, γ, b, ψ, η, s1, s2
end

function unpack_θ(θ, model)
    q0 = view(θ, 1:3)  # configuration
    q1 = view(θ, 4:6)  # configuration
    u = view(θ, 7:9)
    w = view(θ, 10:10)
    h = view(θ, 11:11)
    return q0, q1, u, w, h
end

inequality_indices(model) = collect(model.nq .+ (1:(num_var(model) - model.nq)))

# Model
h = 0.01
T = 500
model = Particle(1.0, 9.81, 1.0, h, 3, 3, 1, 4, 1)

function _r!(r, z, θ, κ)
    q2, γ, b, ψ, η, s1, s2 = unpack_z(z, model)
    q0, q1, u, w, h = unpack_θ(θ, model)

    ϕ = signed_distance(model, q2)        # signed-distance function
    vT = (tangent_jacobian(model, q2) * q2 - tangent_jacobian(model, q1) * q1) / h # tangent velocity

    r[1:3] = dynamics(model, q0, q1, q2, u, γ, b, w, h)
    r[4:4] = s1 - ϕ
    r[5:5] = γ .* s1 - κ
    r[6:9] = vT - η
    r[6:9] .+= ψ
    r[10:10] = s2 - (model.μ * γ .- sum(b))
    r[11:11] = ψ .* s2 - κ
    r[12:15] = b .* η .- κ

    nothing
end

r = zeros(num_var(model))
z = ones(num_var(model))
θ = ones(num_data(model))
κ = 1.0
@code_warntype _r!(r, z, θ, κ)
# @benchmark _r!($r, $z, $θ, $κ)

@variables r_sym[1:15]
@variables z_sym[1:num_var(model)]
@variables θ_sym[1:num_data(model)]
@variables κ_sym[1:1]

parallel = false # ModelingToolkit.MultithreadedForm()
_r!(r_sym, z_sym, θ_sym, κ_sym)
r_sym = simplify.(r_sym)
r! = eval(ModelingToolkit.build_function(r_sym, z_sym, θ_sym, κ_sym,
    parallel = parallel)[2])
rz_exp = ModelingToolkit.sparsejacobian(r_sym, z_sym, simplify = true)
rθ_exp = ModelingToolkit.jacobian(r_sym, θ_sym, simplify = true)
rz_sp = similar(rz_exp, Float64)
rθ_sp = similar(rθ_exp, Float64)
rz! = eval(ModelingToolkit.build_function(rz_exp, z_sym, θ_sym, κ_sym,
    parallel = parallel)[2])
rθ! = eval(ModelingToolkit.build_function(rθ_exp, z_sym, θ_sym, κ_sym,
    parallel = parallel)[2])

@code_warntype r!(r, z, θ, κ)
@benchmark r!($r, $z, $θ, $κ)

@code_warntype rz!(rz_sp, z, θ, κ)
rz!(rz_sp, z, θ, κ)
@benchmark rz!($rz_sp, $z, $θ, $κ)

@code_warntype rθ!(r, z, θ, κ)
rθ!(rθ_sp, z, θ, κ)
@benchmark rθ!($rθ_sp, $z, $θ, $κ)

q0 = @SVector [0.0, 0.0, 1.0]
q1 = @SVector [0.1, -.05, 1.0]
sim = simulator(model, q0, q1, h, T,
    rz = rz_sp,
    rθ = rθ_sp)

@time simulate!(sim)

using Plots
plot(hcat(Array.(sim.q)...)')
