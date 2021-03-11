"""
    particle dynamics
    - 2D particle subject to contact forces

    - configuration: q = (x, z) ∈ R²
    - impact force (magnitude): γ ∈ R₊
    - friction force: β ∈ R²₊
        - friction coefficient: μ ∈ R₊

    Discrete Mechanics and Variational Integrators
        pg. 363
"""
struct Particle2D{T} <: ContactDynamicsModel
    dim::Dimensions
    m::T # mass
    g::T # gravity
    μ_world::T # friction coefficient
	μ_joint::T

	base::BaseMethods
	dyn::DynamicsMethods
	res::ResidualMethods

	spa::SparseStructure

	joint_friction::SVector

	env::Environment
end

function lagrangian(model::Particle2D, q, q̇)
	L = 0.0

	L += 0.5 * model.m * transpose(q̇) * q̇
	L -= model.m * model.g * q[2]

	return L
end

# mass matrix
function M_func(model::Particle2D, q)
    m = model.m

    Diagonal(@SVector [m, m])
end

# gravity
function C_func(model::Particle2D, q, q̇)
    m = model.m
    g = model.g

    @SVector [0.0, m * g]
end

# signed distance function
function ϕ_func(model::Particle2D, q)
    q[2:2] .- model.env.surf(q)
end

# control Jacobian
function B_func(model::Particle2D, q)
    SMatrix{2, 2}([1.0 0.0;
                   0.0 1.0])
end

# disturbance Jacobian
function A_func(model::Particle2D, q)
	SMatrix{2, 2}([1.0 0.0;
                   0.0 1.0])
end

# contact Jacobian
function J_func(model::Particle2D, q)
	SMatrix{2, 2}([1.0 0.0;
                   0.0 1.0])
end

# Model (flat surface)
particle2D = Particle2D(Dimensions(2, 2, 2, 1, 2), 1.0, 9.81, 1.0, 0.0,
	BaseMethods(), DynamicsMethods(), ResidualMethods(),
	SparseStructure(spzeros(0,0),spzeros(0,0)),
	SVector{2}(zeros(2)),
	environment_2D_flat())

# Model (slope)
particle2D_slope = Particle2D(Dimensions(2, 2, 2, 1, 2), 1.0, 9.81, 0.1, 0.0,
	BaseMethods(), DynamicsMethods(), ResidualMethods(),
	SparseStructure(spzeros(0,0),spzeros(0,0)),
	SVector{2}(zeros(2)),
	environment_2D(x -> 0.5 * x[1]))
