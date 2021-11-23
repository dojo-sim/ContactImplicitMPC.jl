"""
    hopper 3D
        orientation representation: modified rodrigues parameters
		similar to Raibert hopper, all mass is located at the body
		s = (px, py, pz, tx, ty, tz, r) = (3d_position, MRP(m1, m2, 0.0), leg_length)
"""
struct Hopper3D{T} <: Model{T}
    nq::Int 
	nu::Int
	nw::Int
	nc::Int

	mb::T # mass of body
    ml::T # mass of leg
    Jb::T # inertia of body
    Jl::T # inertia of leg

    μ_world::T  # coefficient of friction
	μ_joint::T
    g::T # gravity

	orientation::Symbol

	base::BaseMethods
	dyn::DynamicsMethods

	joint_friction::SVector
end

lagrangian(model::Hopper3D, q, q̇) = 0.0

# Kinematics
function kinematics(model::Hopper3D, q)
	p = view(q, 1:3)
	R = eval(model.orientation)(view(q, 4:6)...)
	p + R * [0.0; 0.0; -1.0 * q[7]]
end

# Methods
function M_func(model::Hopper3D, q)
	Diagonal(@SVector [model.mb + model.ml, model.mb + model.ml, model.mb + model.ml,
					   model.Jb + model.Jl, model.Jb + model.Jl, model.Jb + model.Jl,
					   model.ml])
end

function C_func(model::Hopper3D, q, q̇)
	@SVector [0.0, 0.0, (model.mb + model.ml) * model.g, 0.0, 0.0, 0.0, 0.0]
end

function ϕ_func(model::Hopper3D, env::Environment, q)
	SVector{1}(kinematics(model, q)[3] - env.surf(kinematics(model, q)[1:2]))
end

function B_func(::Hopper3D, q)
    rot = view(q, 4:6)
    R = eval(model.orientation)(rot...)
    @SMatrix [0.0 0.0 0.0 R[1,1] R[2,1] R[3,1] 0.0;
              0.0 0.0 0.0 R[1,2] R[2,2] R[3,2] 0.0;
			  R[1,3] R[2,3] R[3,3] 0.0 0.0 0.0 1.0]
end

function A_func(::Hopper3D, q)
    rot = view(q, 4:6)
    R = eval(model.orientation)(rot...)
    @SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0;
              0.0 1.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 1.0 0.0 0.0 0.0 0.0]
end

function J_func(model::Hopper3D, env::Environment, q)
    k(z) = kinematics(model, z)
    ForwardDiff.jacobian(k, q)
end

function contact_forces(model::Hopper3D, env::Environment{<:World, LinearizedCone}, γ1, b1, q2, k)
	m = friction_mapping(env)

	SVector{3}(transpose(rotation(env, k)) * [m * b1; γ1])
end

function velocity_stack(model::Hopper3D, env::Environment{<:World, LinearizedCone}, q1, q2, k, h)
	v = J_func(model, env, q2) * (q2 - q1) / h[1]

	v1_surf = rotation(env, k) * v

	SVector{4}(friction_mapping(env)' * v1_surf[1:2])
end
function get_stride(model::Hopper3D, traj)
    stride = zeros(model.nq)
    stride[1:2] = traj.q[end-1][1:2] - traj.q[1][1:2]
    return stride
end

# Dimensions
nq = 7 # configuration dimension
nu = 3 # control dimension
nw = 3 # disturbance dimension
nc = 1 # number of contact points
nquat = 0 # number of quaternions

# Parameters
g = 9.81 # gravity
μ_world = 1.5 # coefficient of friction
μ_joint = 0.0

# TODO: change to Raibert parameters
mb = 3.0 # body mass
ml = 0.3  # leg mass
Jb = 0.75 # body inertia
Jl = 0.075 # leg inertia

hopper_3D = Hopper3D(nq,nu,nw,nc,
			mb, ml, Jb, Jl,
			μ_world, μ_joint, g,
			:MRP,
			BaseMethods(), DynamicsMethods(),
			SVector{7}(zeros(7)))

function friction_coefficients(model::Hopper3D) 
	return [model.μ_world]
end

function initialize_z!(z, model::Hopper3D, idx::RoboDojo.IndicesZ, q)
    z .= 1.0
    z[idx.q] .= q
end