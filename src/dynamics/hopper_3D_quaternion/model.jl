"""
    hopper 3D
        orientation representation: quaternion
		similar to Raibert hopper, all mass is located at the body
		s = (px, py, pz, q0, q1, q2, q3, r) = (3d_position, quaternion, leg_length)
"""
struct Hopper3DQuaternion{T} <: ContactModel where T
    dim::Dimensions

	mb::T # mass of body
    ml::T # mass of leg
    Jb::T # inertia of body
    Jl::T # inertia of leg

    μ_world::T  # coefficient of friction
	μ_joint::T
    g::T # gravity

	base::BaseMethods
	dyn::DynamicsMethods

	joint_friction::SVector
end

function G_func(::Hopper3DQuaternion, q)
	quat = q[4:7]
	[1.0 0.0 0.0 0.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 0.0 0.0 0.0;
	 0.0 0.0 1.0 0.0 0.0 0.0 0.0;
     zeros(4, 3) attitude_jacobian(quat) zeros(4, 1);
	 zeros(1, 6) 1.0]
end

function Gz_func(model::Hopper3DQuaternion, env, z)
	nz = num_var(model, env)
	ny = nz - 1

	[G_func(model, z) zeros(8, ny - 7);
	 zeros(nz - 8, 7) I]
end

function lagrangian(model::Hopper3DQuaternion, q, q̇)
	0.0
end

# Kinematics
function kinematics(model::Hopper3DQuaternion, q)
	p = q[1:3]
	R = eval(model.orientation)(q[4:7]...)
	p + R * [0.0; 0.0; -1.0 * q[8]]
end

# Methods
function M_func(model::Hopper3DQuaternion, q)
	Diagonal(@SVector [model.mb + model.ml, model.mb + model.ml, model.mb + model.ml,
					   model.Jb + model.Jl, model.Jb + model.Jl, model.Jb + model.Jl,
					   model.ml])
end

function C_func(model::Hopper3DQuaternion, q, q̇)
	ω = q̇[4:6]
	e = cross(ω, Diagonal([model.Jb + model.Jl, model.Jb + model.Jl, model.Jb + model.Jl]) * ω)
	SVector{7}([0.0, 0.0, (model.mb + model.ml) * model.g,
				e[1], e[2], e[3],
			  0.0])
end

function ϕ_func(model::Hopper3DQuaternion, env::Environment, q)
	SVector{1}(kinematics(model, q)[3] - env.surf(kinematics(model, q)[1:2]))
end

function B_func(model::Hopper3DQuaternion, q)
    rot = view(q, 4:7)
    R = eval(model.orientation)(rot...)
    @SMatrix [0.0 0.0 0.0 R[1,1] R[2,1] R[3,1] 0.0;
              0.0 0.0 0.0 R[1,2] R[2,2] R[3,2] 0.0;
			  R[1,3] R[2,3] R[3,3] 0.0 0.0 0.0 1.0]
end

function A_func(::Hopper3DQuaternion, q)
    @SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0;
              0.0 1.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 1.0 0.0 0.0 0.0 0.0]
end

function J_func(model::Hopper3DQuaternion, q)
    k(z) = kinematics(model, z)
    ForwardDiff.jacobian(k, q) * G_func(model, q)
end

function contact_forces(model::Hopper3DQuaternion, env::Environment{R3, LinearizedCone}, γ1, b1, q2, k)
	m = friction_mapping(env)
	SVector{3}(transpose(rotation(env, k)) * [m * b1; γ1])
end

function contact_forces(model::Hopper3DQuaternion, env::Environment{R3, NonlinearCone}, γ1, b1, q2, k)
	SVector{3}(transpose(rotation(env, k)) * [b1; γ1])
end

function velocity_stack(model::Hopper3DQuaternion, env::Environment{R3, LinearizedCone}, q1, q2, k, h)
	p1 = q1[1:3]
	quat1 = q1[4:7]
	r1 = q1[8]

	p2 = q2[1:3]
	quat2 = q2[4:7]
	r2 = q2[8]

	v = J_func(model, q2) * [(p2 - p1) / h[1]; ω_finite_difference(quat1, quat2, h[1]); (r2 - r1) / h[1]]

	v1_surf = rotation(env, k) * v

	SVector{4}(transpose(friction_mapping(env)) * v1_surf[1:2])
end

function velocity_stack(model::Hopper3DQuaternion, env::Environment{R3, NonlinearCone}, q1, q2, k, h)
	p1 = q1[1:3]
	quat1 = q1[4:7]
	r1 = q1[8]

	p2 = q2[1:3]
	quat2 = q2[4:7]
	r2 = q2[8]

	v = J_func(model, q2) * [(p2 - p1) / h[1]; ω_finite_difference(quat1, quat2, h[1]); (r2 - r1) / h[1]]

	v1_surf = rotation(env, k) * v

	SVector{2}(v1_surf[1:2])
end

function dynamics(model::Hopper3DQuaternion, h, q0, q1, u1, w1, λ1, q2)

	p0 = q0[1:3]
	quat0 = q0[4:7]
	r0 = q0[8]

	p1 = q1[1:3]
	quat1 = q1[4:7]
	r1 = q1[8]

	p2 = q2[1:3]
	quat2 = q2[4:7]
	r2 = q2[8]

	ω1 = ω_finite_difference(quat0, quat1, h)
	ω2 = ω_finite_difference(quat1, quat2, h)

	# evalutate at midpoint
	qm1 = [0.5 * (p0 + p1); zeros(4); 0.5 * (r0 + r1)]
    vm1 = [(p1 - p0) / h[1]; ω1; (r1 - r0) / h[1]]
    qm2 = [0.5 * (p1 + p2); zeros(4); 0.5 * (r1 + r2)]
    vm2 = [(p2 - p1) / h[1]; ω2; (r2 - r1) / h[1]]

	ω1 = ω_finite_difference(quat0, quat1, h)
	ω2 = ω_finite_difference(quat1, quat2, h)
	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	J = Diagonal([model.Jb + model.Jl, model.Jb + model.Jl, model.Jb + model.Jl])
	d = [(0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2)[1:3];
		 -1.0 * (J * ω2 * sqrt(4.0 / h[1]^2.0 - transpose(ω2) * ω2)
			+ cross(ω2, J * ω2)
			- J * ω1 * sqrt(4.0 / h[1]^2.0 - transpose(ω1) * ω1)
			+ cross(ω1, J * ω1));
		 (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2)[7]]

	return (d
		+ transpose(B_fast(model, q2)) * u1
		+ transpose(A_fast(model, q2)) * w1
		+ transpose(J_fast(model, q2)) * λ1)
end

function get_stride(model::Hopper3DQuaternion, traj::ContactTraj)
    stride = zeros(SizedVector{model.dim.q})
    stride[1:2] = traj.q[end-1][1:2] - traj.q[1][1:2]
    return stride
end

# Dimensions
nq = 8 # configuration dimension
nu = 3 # control dimension
nw = 3 # disturbance dimension
nc = 1 # number of contact points

# Parameters
g = 9.81 # gravity
μ_world = 1.0 # coefficient of friction
μ_joint = 0.0

# TODO: change to Raibert parameters
mb = 3.0 # body mass
ml = 0.3  # leg mass
Jb = 0.75 # body inertia
Jl = 0.075 # leg inertia

hopper_3D_quaternion = Hopper3DQuaternion(Dimensions(nq, nu, nw, nc),
			mb, ml, Jb, Jl,
			μ_world, μ_joint, g,
			BaseMethods(), DynamicsMethods(),
			SVector{8}(zeros(8)))
