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
struct Particle{T} <: ContactDynamicsModel
    dim::Dimensions
    m::T # mass
    g::T # gravity
    μ_world::T # friction coefficient
	μ_joint::T

	base::BaseMethods
	dyn::DynamicsMethods
	res::ResidualMethods
	linearized::ResidualMethods

	spa::SparseStructure

	joint_friction::SVector

	env::Environment
end

function lagrangian(model::Particle, q, q̇)
	L = 0.0

	L += 0.5 * model.m * transpose(q̇) * q̇
	L -= model.m * model.g * q[3]

	return L
end

function kinematics(::Particle, q)
	return q
end

# mass matrix
function M_func(model::Particle, q)
    m = model.m

    Diagonal(@SVector [m, m, m])
end

# gravity
function C_func(model::Particle, q, q̇)
    m = model.m
    g = model.g

    @SVector [0.0, 0.0, m * g]
end

# signed distance function
function ϕ_func(model::Particle, q)
    # q[3:3] .- model.env.surf(q)
	SVector{model.dim.c}(q[3:3] - model.env.surf(q[1:2]))
end

# control Jacobian
function B_func(model::Particle, q)
    SMatrix{3, 3}([1.0 0.0 0.0;
                   0.0 1.0 0.0;
                   0.0 0.0 1.0])
end

# disturbance Jacobian
function A_func(model::Particle, q)
    SMatrix{3, 3}([1.0 0.0 0.0;
                   0.0 1.0 0.0;
                   0.0 0.0 1.0])
end

# contact Jacobian
function J_func(model::Particle, q)
    SMatrix{3, 3}([1.0 0.0 0.0;
                   0.0 1.0 0.0;
				   0.0 0.0 1.0])
end

# Model (flat surface)
particle = Particle(Dimensions(3, 3, 3, 1, 4), 1.0, 9.81, 1.0, 0.0,
	BaseMethods(), DynamicsMethods(), ResidualMethods(), ResidualMethods(),
	SparseStructure(spzeros(0,0),spzeros(0,0)),
	SVector{3}(zeros(3)),
	environment_3D_flat())

# Model (quadratic bowl)
particle_quadratic = Particle(Dimensions(3, 3, 3, 1, 4), 1.0, 9.81, 0.1, 0.0,
	BaseMethods(), DynamicsMethods(), ResidualMethods(), ResidualMethods(),
	SparseStructure(spzeros(0,0),spzeros(0,0)),
	SVector{3}(zeros(3)),
	environment_3D(x -> [transpose(x[1:2]) * x[1:2]]))

# no gravity
particle_sinusoidal = Particle(Dimensions(3, 3, 3, 1, 4), 1.0, 0.0 * 9.81, 1.0, 0.0,
	BaseMethods(), DynamicsMethods(), ResidualMethods(), ResidualMethods(),
	SparseStructure(spzeros(0,0),spzeros(0,0)),
	SVector{3}(zeros(3)),
	environment_3D(x -> [sin(x[1])+sin(x[2])]))
