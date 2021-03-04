include("interior_point.jl")
include("simulator.jl")

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
    q[3]
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
    model.nq + model.nq + model.nu + 2
end

function z_initialize!(z, model, q1)
    z .= 1.0
    z[1:3] = q1
end

function θ_initialize!(θ, model, q0, q1, u, w)
    θ[1:3] = q0
    θ[4:6] = q1
    θ[7:9] = u
end

function unpack_z(z, model)
    q2 = view(z, 1:3)  # configuration
    γ = view(z, 4:4)   # normal impulse
    b = view(z, 5:8)   # friction impulse (double parameterized)
    ψ = z[9]           # friction cone dual (magnitude of tangential velocity)
    η = view(z, 10:13) # friction impulse duals
    s1 = z[14]         # signed-distance slack
    s2 = z[15]         # friction cone slack

    return q2, γ, b, ψ, η, s1, s2
end

function unpack_θ(θ, model)
    q0 = view(θ, 1:3)  # configuration
    q1 = view(θ, 4:6)  # configuration
    u = view(θ, 7:9)
    w = nothing
    return q0, q1, u, w
end

inequality_indices(model) = collect(model.nq .+ (1:(model.nq - num_var(model))))

function r!(r, z, data)
    # unpack
    θ = data.θ
    κ = data.κ
    model = data.info[:model]

    q2, γ, b, ψ, η, s1, s2 = unpack_z(z, model)
    q0, q1, u, w = unpack_θ(θ, model)

    ϕ = signed_distance(model, q2)        # signed-distance function
    vT = (tangent_jacobian(model, q2) * q2 - tangent_jacobian(model, q1) * q1) / h # tangent velocity

    r[1:3] = dynamics(model, q0, q1, q2, u, γ, b, w, h)
    r[4] = s1 - ϕ
    r[5] = γ[1] * s1 - κ
    r[6:9] = vT - η
    r[6:9] .+= ψ
    r[10] = s2 - (model.μ * γ[1] - sum(b))
    r[11] = ψ * s2 - κ
    r[12:15] = b .* η .- κ

    nothing
end

function rz!(rz, z, data)
    _r(a, b) = r!(a, b, data)
    ForwardDiff.jacobian!(rz, _r, data.r, z)
end

function rθ!(rθ, z, data)
    _r(a, c) = r!(a, z, c, κ)
    ForwardDiff.jacobian!(rθ, _r, data.r, data.θ)
end

# Model
h = 0.01
T = 500
model = Particle(1.0, 9.81, 1.0, h, 3, 3, 1, 4, 0)
q0 = @SVector [0.0, 0.0, 1.0]
q1 = @SVector [0.1, 0.0, 1.0]
sim = simulator(model, q0, q1, h, T)
@time simulate!(sim)

sim.ip_data.δz
sim.q
sim.u
sim.γ
sim.dq2dq0
sim.dγdq0

SizedMatrix{100,100}(zeros(100,100))

xx = [@SVector ones(3) for t = 1:5]
