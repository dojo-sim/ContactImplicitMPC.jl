"""
    Unicycle
"""
mutable struct Unicycle{T} <: ContactModel
    dim::Dimensions

	mb::T # body mass
    Jb::T # body inertia

	mw::T # wheel mass
	Jw::T # wheel inertia
	rw::T # wheel radius

	l0::T # suspension length
	k::T # suspension stiffness

    μ_world::T  # coefficient of friction
    μ_joint::T  # gravity
	g::T

	base::BaseMethods
	dyn::DynamicsMethods

	joint_friction::SVector
end


# Kinematics
function kinematics(model::Unicycle, q)
	x = q[1]
	z = q[2]
	l = q[3]
	return [x; z - l]
end

function lagrangian(model::Unicycle, q, q̇)
	L = 0.0

	x = q[1]
	z = q[2]
	l = q[3]
	θw = q[4]

	ẋ = q̇[1]
	ż = q̇[2]
	l̇ = q̇[3]
	θdw = q̇[4]

	L += 0.5 * model.mb * (ẋ^2 + ż^2)
	L += 0.5 * model.mw * l̇^2
	L += 0.5 * model.Jw * θdw^2

	L -= model.mb * model.g * z
	L -= model.mw * model.g * (z - l)
	L -= 0.5 * model.k  * (l - model.l0)^2
	return L
end

# Methods
function M_func(model::Unicycle, q)
	Diagonal(@SVector [model.mb + model.mw,
					   model.mb + model.mw,
					   model.mw,
					   model.Jw,
					   ])
 end

function C_func(model::Unicycle, q, q̇) #false
	@SVector [0.0,
			  (model.mb + model.mw) * model.g,
			  0.0,
			  0.0,
			  ]
end

function ϕ_func(model::Unicycle, env::Environment, q)
	# assumes point of contact is on the suspension axis (flat ground)
    SVector{1}(q[2] - q[3] - model.rw - env.surf(q[1]))
end

function J_func(::Unicycle, env::Environment, q)
    @SMatrix [-1.0/model.rw  0.0  0.0  1.0;
		       0.0           1.0 -1.0  0.0]
end

function B_func(model::Unicycle, q)
	@SMatrix [ 0.0  0.0  1.0  0.0; # /100 for scaling
			   0.0  0.0  0.0  1.0/100]
end

function A_func(::Unicycle, q)
	@SMatrix [ 1.0  0.0  0.0  0.0;
			   0.0  1.0  0.0  0.0]
end

# Parameters
g = 9.81 # gravity
μ_world = 0.1 # coefficient of friction
μ_joint = 10.0

mb = 1.0 # body mass
Jb = 0.1 # body inertia

mw = 1.0 # wheel mass
Jw = 0.1 # wheel inertia
rw = 0.1 # wheel radius

l0 = 0.5 # suspension length
k = 300.0 # suspension stiffness

# Dimensions
nq = 4
nu = 2
nw = 2
nc = 1

unicycle = Unicycle(Dimensions(nq, nu, nw, nc),
			   mb, Jb, mw, Jw, rw, l0, k,
			   μ_world, μ_joint, g,
			   BaseMethods(), DynamicsMethods(),
			   SVector{4}(μ_joint * [0.0; 0.0; 1.0; 0.0]))
