mutable struct RigidBody{T} <: ContactDynamicsModel
    dim::Dimensions
    m::T # mass
	J::Vector{T} # inertia
    g::T # gravity
    μ_world::T # friction coefficient

	base::BaseMethods
	dyn::DynamicsMethods
	con::ContactMethods
	res::ResidualMethods
	linearized::ResidualMethods

	spa::SparseStructure

	joint_friction::SVector

	env::Environment
end

function kinematics(::RigidBody, q)
	return q[1:3]
end

function lagrangian(model::RigidBody, q, q̇)
	m = 1.0
	J = Diagonal([1.0, 1.0, 1.0])

	p = q[1:3]
	quat = q[4:7]
	v = q̇[1:3]
	ω = q̇[4:6]

	L = 0.0

	# linear
	L += 0.5 * m * transpose(v) * v
	L -= m * 9.81 * p[3]

	# angular
	L += 0.5 * transpose(ω) * J * ω

	return L
end

# mass matrix
function M_func(model::RigidBody, q)
    m = model.m
	J = model.J
    Diagonal(@SVector [m, m, m, J...])
end

# gravity
function C_func(model::RigidBody, q, q̇)
    m = model.m
    g = model.g

    @SVector [0.0, 0.0, m * g, 0.0, 0.0, 0.0]
end

# signed distance function
function ϕ_func(model::RigidBody, q)
	SVector{1}(q[3] - model.env.surf(q[1:2]))
end


# control Jacobian
function B_func(model::RigidBody, q)
	quat = q[4:7]
    SMatrix{6, 6}([1.0 0.0 0.0 0.0 0.0 0.0
                   0.0 1.0 0.0 0.0 0.0 0.0
                   0.0 0.0 1.0 0.0 0.0 0.0;
				   zeros(3, 3) quaternion_rotation_matrix(quat)])
end

# disturbance Jacobian
function A_func(model::RigidBody, q)
	SMatrix{6, 6}([1.0 0.0 0.0 0.0 0.0 0.0
				   0.0 1.0 0.0 0.0 0.0 0.0
				   0.0 0.0 1.0 0.0 0.0 0.0;
				   zeros(3, 3) quaternion_rotation_matrix(quat)])
end

# contact Jacobian
function J_func(model::RigidBody, q)
    SMatrix{3, 6}([1.0 0.0 0.0 0.0 0.0 0.0
				   0.0 1.0 0.0 0.0 0.0 0.0
				   0.0 0.0 1.0 0.0 0.0 0.0])
end

function contact_forces(model::RigidBody, γ1, b1, q2, k)
	m = friction_mapping(model.env)

	SVector{3}(transpose(rotation(model.env, k)) * [m * b1; γ1])
end

function velocity_stack(model::RigidBody, q1, q2, k, h)
	# k = kinematics(model, q2)
	v = J_func(model, q2) * (q2 - q1) / h[1]

	v1_surf = rotation(model.env, k) * v

	SVector{4}(friction_mapping(model.env)' * v1_surf[1:2])
end

# quaternion midpoint http://web.cse.ohio-state.edu/~parent.1/classes/682/Lectures/Lectures07/L05_Interpolation/Quaternions.pdf
function dynamics(model::RigidBody, h, q0, q1, u1, w1, λ1, q2)

	p0 = q0[1:3]
	quat0 = q0[4:7]

	p1 = q1[1:3]
	quat1 = q1[4:7]

	p2 = q2[1:3]
	quat2 = q2[4:7]

	# evalutate at midpoint
	qm1 = [0.5 * (p0 + p1); (quat0 + quat1) ./ norm(quat0 + quat1)]
    vm1 = [(p1 - p0) / h[1]; ω_finite_difference(quat0, quat1, h[1])]
    qm2 = [0.5 * (p1 + p2); (quat1 + quat2) ./ norm(quat1 + quat2)]
    vm2 = [(p2 - p1) / h[1]; ω_finite_difference(quat1, quat2, h[1])]

	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	return [(0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
		+ transpose(B_fast(model, qm2)) * u1
		+ transpose(A_fast(model, qm2)) * w1
		+ transpose(J_fast(model, q2)) * λ1
		- h[1] * model.joint_friction .* vm2);
		sqrt(quat2[1]^2.0 + quat2[2]^2.0 + quat2[3]^2.0 + quat2[4]^2.0) - 1.0]

end

# Model (flat surface)
model = RigidBody(Dimensions(3, 3, 3, 1, 4), 1.0, 9.81, 1.0, 0.0,
	BaseMethods(), DynamicsMethods(), ContactMethods(),
	ResidualMethods(), ResidualMethods(),
	SparseStructure(spzeros(0,0),spzeros(0,0)),
	SVector{3}(zeros(3)),
	environment_3D_flat())

nq = 7
nv = 6

lagrangian(model, rand(nq), rand(nv))

# Symbolics
@variables q[1:nq]
@variables q̇[1:nv]

G = [1.0 0.0 0.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 0.0 0.0;
	 0.0 0.0 1.0 0.0 0.0 0.0;
     0.0 0.0 0.0 -transpose(q[1:3]);
	 zeros(3, 3) Diagonal(q[4] * ones(3)) + skew(q[1:3])]

# Lagrangian
L = lagrangian(model, q, q̇)
L = Symbolics.simplify.(L)

_dLq = Symbolics.gradient(L, q, simplify=true)
dLq = G' * _dLq

dLq̇ = Symbolics.gradient(L, q̇, simplify=true)
ddL = Symbolics.hessian(L, [q; q̇], simplify=true)
ddLq̇q = ddL[nq .+ (1:nv), 1:nq] * G
M = ddL[nq .+ (1:nv), nq .+ (1:nv)]

# Coriolis and Centrifugal forces Jacobians
C = ddLq̇q * q̇ - dLq
C = Symbolics.simplify.(C)
