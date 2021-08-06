"""
    Hopper2D
    	model inspired by "Dynamically Stable Legged Locomotion"
		s = (x, z, t, r)
			x - lateral position
			z - vertical position
			t - body orientation
			r - leg length
"""
mutable struct Hopper2D{T} <: ContactModel
    dim::Dimensions

    mb::T # mass of body
    ml::T # mass of leg
    Jb::T # inertia of body
    Jl::T # inertia of leg

    μ_world::T  # coefficient of friction
    μ_joint::T  # gravity
	g::T

	base::BaseMethods
	dyn::DynamicsMethods

	joint_friction::SVector
end

lagrangian(model::Hopper2D, q, q̇) = 0.0

# Kinematics
function kinematics(::Hopper2D, q)
	[q[1] + q[4] * sin(q[3]),
	 q[2] - q[4] * cos(q[3])]
end

# Methods
function M_func(model::Hopper2D, q)
	Diagonal(@SVector [model.mb + model.ml,
					   model.mb + model.ml,
					   model.Jb + model.Jl,
					   model.ml])
 end

function C_func(model::Hopper2D, q, q̇)
	@SVector [0.0,
			  (model.mb + model.ml) * model.g,
			  0.0,
			  0.0]
end

function ϕ_func(model::Hopper2D, env::Environment, q)
    SVector{1}(q[2] - q[4] * cos(q[3]) - env.surf(q[1] + q[4] * sin(q[3])))
end

function J_func(::Hopper2D, env::Environment, q)
    @SMatrix [1.0 0.0 (q[4] * cos(q[3])) sin(q[3]);
		      0.0 1.0 (q[4] * sin(q[3])) (-1.0 * cos(q[3]))]
end

function B_func(::Hopper2D, q)
	@SMatrix [0.0 0.0 1.0 0.0;
             -sin(q[3]) cos(q[3]) 0.0 1.0]
end

function A_func(::Hopper2D, q)
	@SMatrix [1.0 0.0 0.0 0.0;
	          0.0 1.0 0.0 0.0]
end

function contact_forces(model::Hopper2D, env::Environment{<:World, LinearizedCone}, γ1, b1, q2, k)
	m = friction_mapping(env)

	SVector{2}(transpose(rotation(env, k)) * [m * b1; γ1])
end

function velocity_stack(model::Hopper2D, env::Environment{<:World, LinearizedCone}, q1, q2, k, h)
	v = J_func(model, env, q2) * (q2 - q1) / h[1]

	v1_surf = rotation(env, k) * v

	SVector{2}([v1_surf[1]; -v1_surf[1]])
end

# Working Parameters
gravity = 9.81 # gravity
μ_world = 0.8 # coefficient of friction
μ_joint = 0.0

# TODO: change to Raibert parameters
mb = 3.0 # body mass
ml = 0.3  # leg mass
Jb = 0.75 # body inertia
Jl = 0.075 # leg inertia

# Dimensions
nq = 4
nu = 2
nw = 2
nc = 1
nb = 2

hopper_2D = Hopper2D(Dimensions(nq, nu, nw, nc),
			   mb, ml, Jb, Jl,
			   μ_world, μ_joint, gravity,
			   BaseMethods(), DynamicsMethods(),
			   SVector{4}(zeros(4)))
