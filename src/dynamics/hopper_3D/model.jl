"""
    hopper 3D
        orientation representation: modified rodrigues parameters
		similar to Raibert hopper, all mass is located at the body
		s = (px, py, pz, tx, ty, tz, r) = (3d_position, 3d_orientation_MRP, leg_length)
"""
struct Hopper3D{T} <: ContactDynamicsModel where T
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
	con::ContactMethods
	res::ResidualMethods
	linearized::ResidualMethods

	spa::SparseStructure

	joint_friction::SVector

	env::Environment
end

lagrangian(model::Hopper3D, q, q̇) = 0.0

# Kinematics
function kinematics(::Hopper3D, q)
	p = view(q, 1:3)
	R = MRP(view(q, 4:6)...)
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

function ϕ_func(model::Hopper3D, q)
	SVector{1}(kinematics(model, q)[3] - model.env.surf(kinematics(model, q)[1:2]))
end

function B_func(::Hopper3D, q)
    rot = view(q, 4:6)
    R = MRP(rot...)
    @SMatrix [0.0 0.0 0.0 R[1,1] R[2,1] R[3,1] 0.0;
              0.0 0.0 0.0 R[1,2] R[2,2] R[3,2] 0.0;
			  R[1,3] R[2,3] R[3,3] 0.0 0.0 0.0 1.0]
end

function A_func(::Hopper3D, q)
    rot = view(q, 4:6)
    R = MRP(rot...)
    @SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0;
              0.0 1.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 1.0 0.0 0.0 0.0 0.0]
end

function J_func(model::Hopper3D, q)
    k(z) = kinematics(model, z)
    ForwardDiff.jacobian(k, q)
end

function contact_forces(model::Hopper3D, γ1, b1, q2, k)
	# k = kinematics(model, q2)
	m = friction_mapping(model.env)

	SVector{3}(transpose(rotation(model.env, k)) * [m * b1; γ1])
end

function velocity_stack(model::Hopper3D, q1, q2, k, h)
	# k = kinematics(model, q2)
	v = J_func(model, q2) * (q2 - q1) / h[1]

	v1_surf = rotation(model.env, k) * v

	SVector{4}([v1_surf[1:2]; -v1_surf[1:2]])
end

# Dimensions
nq = 7 # configuration dimension
nu = 3 # control dimension
nw = 3 # disturbance dimension
nc = 1 # number of contact points
nf = 4 # number of faces for friction cone pyramid
nb = nc * nf

# # Parameters
# g = 9.81 # gravity
# μ = 1.0  # coefficient of friction
# mb = 1.0 # body mass
# ml = 0.1  # leg mass
# Jb = 0.25 # body inertia
# Jl = 0.025 # leg inertia

# Parameters
g = 9.81 # gravity
μ_world = 0.8 # coefficient of friction
μ_joint = 0.0

# TODO: change to Raibert parameters
mb = 3.0 # body mass
ml = 0.3  # leg mass
Jb = 0.75 # body inertia
Jl = 0.075 # leg inertia


hopper_3D = Hopper3D(Dimensions(nq, nu, nw, nc, nb),
			mb, ml, Jb, Jl,
			μ_world, μ_joint, g,
			BaseMethods(), DynamicsMethods(), ContactMethods(),
			ResidualMethods(), ResidualMethods(),
			SparseStructure(spzeros(0, 0), spzeros(0, 0)),
			SVector{7}(zeros(7)),
			environment_3D_flat())

hopper_3D_sinusoidal = Hopper3D(Dimensions(nq, nu, nw, nc, nb),
			mb, ml, Jb, Jl,
			μ_world, μ_joint, g,
			BaseMethods(), DynamicsMethods(), ContactMethods(),
			ResidualMethods(), ResidualMethods(),
			SparseStructure(spzeros(0, 0), spzeros(0, 0)),
			SVector{7}(zeros(7)),
			environment_3D(x -> 0.10 * sin(2π * x[1])),
		    )
