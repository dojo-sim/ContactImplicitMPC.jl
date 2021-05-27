# rigid body with quaternion representation (sphere w/ radius r)
mutable struct RigidBody{T} <: ContactModel
    dim::Dimensions
    m::T # mass
	J::Vector{T} # inertia
    g::T # gravity
    μ_world::T # friction coefficient
	r::T # radius

	base::BaseMethods
	dyn::DynamicsMethods

	joint_friction::SVector
end

function kinematics(model::RigidBody, q)
	return q[1:3] - [0.0; 0.0; model.r]
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
	L -= m * model.g * p[3]

	# angular
	L += 0.5 * transpose(ω) * J * ω

	return L
end

# mass matrix
function M_func(model::RigidBody, q)
    m = model.m
	J = model.J
    Diagonal(SVector{6}([m, m, m, J...]))
end

# gravity
function C_func(model::RigidBody, q, q̇)
    m = model.m
    g = model.g

	ω = q̇[4:6]

    SVector{6}([0.0, 0.0, m * g, cross(ω, Diagonal(model.J) * ω)...]) #TODO confirm cross product
end

# signed distance function
function ϕ_func(model::RigidBody, env::Environment, q)
	SVector{1}(q[3] - model.r - env.surf(q[1:2]))
end

# control Jacobian
function B_func(model::RigidBody, q)
    SMatrix{6, 6}([1.0 0.0 0.0 0.0 0.0 0.0;
                   0.0 1.0 0.0 0.0 0.0 0.0;
                   0.0 0.0 1.0 0.0 0.0 0.0;
				   0.0 0.0 0.0 1.0 0.0 0.0;
                   0.0 0.0 0.0 0.0 1.0 0.0;
                   0.0 0.0 0.0 0.0 0.0 1.0])
end

# disturbance Jacobian
function A_func(model::RigidBody, q)
	SMatrix{3, 6}([1.0 0.0 0.0 0.0 0.0 0.0
				   0.0 1.0 0.0 0.0 0.0 0.0
				   0.0 0.0 1.0 0.0 0.0 0.0])
end

# contact Jacobian
function J_func(model::RigidBody, q)
	r = [0.0; 0.0; -model.r]
    SMatrix{3, 6}([Diagonal(ones(3)) -1.0 * skew(r)])
end

function contact_forces(model::RigidBody, env::Environment{<:World, LinearizedCone}, γ1, b1, q2, k)
	m = friction_mapping(env)

	SVector{3}(transpose(rotation(env, k)) * [m * b1; γ1])
end

function contact_forces(model::RigidBody, env::Environment{<:World, NonlinearCone}, γ1, b1, q2, k)
	SVector{3}(transpose(rotation(env, k)) * [b1; γ1])
end

function velocity_stack(model::RigidBody, env::Environment{<:World, LinearizedCone}, q1, q2, k, h)
	p1 = q1[1:3]
	quat1 = q1[4:7]

	p2 = q2[1:3]
	quat2 = q2[4:7]

	v = J_func(model, q2) * [(p2 - p1) / h[1]; ω_finite_difference(quat1, quat2, h[1])]

	v1_surf = rotation(env, k) * v

	SVector{4}(friction_mapping(env)' * v1_surf[1:2])
end

function velocity_stack(model::RigidBody, env::Environment{<:World,NonlinearCone}, q1, q2, k, h)
	nc = model.dim.c
	ne = dim(env)

	p1 = q1[1:3]
	quat1 = q1[4:7]

	p2 = q2[1:3]
	quat2 = q2[4:7]

	v = J_func(model, q2) * [(p2 - p1) / h[1]; ω_finite_difference(quat1, quat2, h[1])]

	v_surf = [rotation(env, k[(i-1) * ne .+ (1:ne)]) * v[(i-1) * ne .+ (1:ne)] for i = 1:nc]
	vT_stack = vcat([v_surf[i][1:ne-1] for i = 1:nc]...)
end

function dynamics(model::RigidBody, h, q0, q1, u1, w1, λ1, q2)

	p0 = q0[1:3]
	quat0 = q0[4:7]

	p1 = q1[1:3]
	quat1 = q1[4:7]

	p2 = q2[1:3]
	quat2 = q2[4:7]

	# evalutate at midpoint
	qm1 = [0.5 * (p0 + p1); zeros(4)]
    vm1 = [(p1 - p0) / h[1]; zeros(3)]
    qm2 = [0.5 * (p1 + p2); zeros(4)]
    vm2 = [(p2 - p1) / h[1]; zeros(3)]

	ω1 = ω_finite_difference(quat0, quat1, h)
	ω2 = ω_finite_difference(quat1, quat2, h)

	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	d = [(0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2)[1:3];
		 -1.0 * (Diagonal(model.J) * ω2 * sqrt(4.0 / h[1]^2.0 - transpose(ω2) * ω2)
			+ cross(ω2, Diagonal(model.J) * ω2)
			- Diagonal(model.J) * ω1 * sqrt(4.0 / h[1]^2.0 - transpose(ω1) * ω1)
			+ cross(ω1, Diagonal(model.J) * ω1))]

	return (d
		+ transpose(B_fast(model, qm2)) * u1
		+ transpose(A_fast(model, qm2)) * w1
		+ transpose(J_fast(model, q2)) * λ1)
end

function G_func(::RigidBody, q)
	quat = q[4:7]
	[1.0 0.0 0.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 0.0 0.0;
	 0.0 0.0 1.0 0.0 0.0 0.0;
     zeros(4, 3) attitude_jacobian(quat)]
end

function Gz_func(model::RigidBody, env, z)
	nz = num_var(model, env)
	ny = nz - 1

	[G_func(model, z) zeros(7, ny - 6);
	 zeros(nz - 7, 6) I]
end

# Model
rigidbody = RigidBody(Dimensions(7, 6, 3, 1),
	10.0, [0.01, 0.01, 0.01], 9.81, 1.0, 0.25,
	BaseMethods(), DynamicsMethods(),
	SVector{6}(zeros(6)))

#
# @variables x
# @variables y
# @variables z
# @variables λ
#
# @variables px
# @variables py
# @variables pz
#
# L = (px - x)^2.0 + (py - y)^2.0 + (pz - z)^2.0 + λ * (z - x^2.0 - y^2.0)
#
# dL = Symbolics.gradient(L, [x; y; z; λ])
# dL = simplify.(dL)
#
# ddL = Symbolics.hessian(L, [x; y; z; λ])
