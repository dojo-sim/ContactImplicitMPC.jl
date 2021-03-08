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
struct Particle <: ContactDynamicsModel
    dim
    m # mass
    g # gravity
    μ_world # friction coefficient

	base
	dyn
	res
end

function lagrangian(model::Particle, q, q̇)
	L = 0.0

	L += 0.5 * model.m * transpose(q̇) * q̇
	L -= model.m * model.g * q[3]

	return L
end

# mass matrix
function M_func(model, q)
    m = model.m

    Diagonal(@SVector [m, m, m])
end

# gravity
function C_func(model, q, q̇)
    m = model.m
    g = model.g

    @SVector [0.0, 0.0, m * g]
end

# signed distance function
function ϕ_func(model, q)
    q[3:3]
end

# control Jacobian
function B_func(model, q)
    SMatrix{3, 3}([1.0 0.0 0.0;
                   0.0 1.0 0.0;
                   0.0 0.0 1.0])
end

# normal Jacobian
function N_func(model, q)
    SMatrix{1, 3}([0.0, 0.0, 1.0])
end

# tangent Jacobian
function P_func(model, q)
    SMatrix{4, 3}([1.0 0.0 0.0;
                   0.0 1.0 0.0;
                   -1.0 0.0 0.0;
                   0.0 -1.0 0.0])
end

function num_var(model)
    model.dim.q + model.dim.c + model.dim.b + model.dim.c + model.dim.b + 2 * model.dim.c
end

function num_data(model)
    model.dim.q + model.dim.q + model.dim.u + model.dim.w + 1
end

function z_initialize!(z, model, q1)
    z .= 1.0
    z[1:3] = q1
end

function θ_initialize!(θ, model, q0, q1, u, w, h)
    θ[1:3] = q0
    θ[4:6] = q1
    θ[7:9] = u
    θ[10:12] = w
    θ[13:13] .= h
end

function unpack_z(model, z)
    q2 = z[1:3]  # configuration
    γ = z[4:4]   # normal impulse
    b = z[5:8]   # friction impulse (double parameterized)
    ψ = z[9:9]   # friction cone dual (magnitude of tangential velocity)
    η = z[10:13] # friction impulse duals
    s1 = z[14:14]# signed-distance slack
    s2 = z[15:15]# friction cone slack

    return q2, γ, b, ψ, η, s1, s2
end

function unpack_θ(model, θ)
    q0 = θ[1:3]  # configuration
    q1 = θ[4:6]  # configuration
    u = θ[7:9]
    w = θ[10:12]
    h = θ[13:13]
    return q0, q1, u, w, h
end

# Model
particle = Particle(Dimensions(3, 3, 3, 1, 4), 1.0, 9.81, 1.0,
	BaseMethods(), DynamicsMethods(), ResidualMethods())

function dynamics(model::Particle, h, q0, q1, u1, w1, γ1, b1, q2)

	v1 = (q2 - q1) / h[1]

	return (1.0 / h[1] *
		  (M_func(model, q0) * (q1 - q0)
		- M_func(model, q1) * (q2 - q1))
		+ transpose(B_func(model, q2)) * u1
		+ transpose(N_func(model, q2)) * γ1
		+ transpose(P_func(model, q2)) * b1
		- h[1] * C_func(model, q2, v1))
end

function residual(model::Particle, z, θ, κ)
	nc = model.dim.c
	nb = model.dim.b

	q0, q1, u1, w1, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ, η, s1, s2 = unpack_z(model, z)

	ϕ = ϕ_func(model, q2)
	vT = (P_func(model, q2) * q2 - P_func(model, q1) * q1) / h[1]

	[dynamics(model, h, q0, q1, u1, w1, γ1, b1, q2);
	 s1 - ϕ;
	 γ1 .* s1 .- κ;
	 vT + ψ[1] .* ones(nb) - η;
	 s2 .- (model.μ_world * γ1
	   .- sum(b1));
	 ψ .* s2 .- κ;
	 b1 .* η .- κ]
end

# fast_expressions!(particle, "particle", generate = false)
