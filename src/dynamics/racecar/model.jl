"""
    RaceCar
"""
mutable struct RaceCar{T} <: ContactModel
    dim::Dimensions

    mb::T # body mass
	mw::T # wheel mass
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
function kinematics(model::RaceCar, q)
	x = q[1]
	z = q[2]
	l = q[3]
	return [x; z - l]
end

function lagrangian(model::RaceCar, q, q̇)
	L = 0.0

	x = q[1]
	z = q[2]
	l = q[3]
	ẋ = q̇[1]
	ż = q̇[2]
	l̇ = q̇[3]
	L += 0.5 * model.mb * (ẋ^2 + ż^2)
	L += 0.5 * model.mw * l̇^2

	L -= model.mb * model.g * z
	L -= model.mw * model.g * (z - l)
	L -= 0.5 * model.k  * (l - model.l0)^2
	return L
end

# Methods
function M_func(model::RaceCar, q)
	Diagonal(@SVector [model.mb + model.mw,
					   model.mb + model.mw,
					   model.mw,
					   ])
 end

function C_func(model::RaceCar, q, q̇)
	@SVector [0.0,
			  (model.mb + model.mw) * model.g,
			  0.0,
			  ]
end

function ϕ_func(model::RaceCar, env::Environment, q)
    SVector{1}(q[2] - q[3] - env.surf(q[1]))
end

function J_func(::RaceCar, env::Environment, q)
    @SMatrix [1.0 0.0  0.0;
		      0.0 1.0 -1.0]
end

function B_func(model::RaceCar, q)
	@SMatrix [0.0 0.0 1.0;]
end

function A_func(::RaceCar, q)
	@SMatrix [1.0 0.0 0.0;
			  0.0 1.0 0.0]
end

# Parameters
g = 9.81 # gravity
μ_world = 0.5  # coefficient of friction
μ_joint = 0.0

mb = 1.0 # body mass
mw = 0.1 # whell mass
l0 = 1.0 # suspension length
k = 30.0 # suspension stiffness

# Dimensions
nq = 3
nu = 1
nw = 2
nc = 1

racecar = RaceCar(Dimensions(nq, nu, nw, nc),
			   mb, mw, l0, k,
			   μ_world, μ_joint, g,
			   BaseMethods(), DynamicsMethods(),
			   SVector{3}(μ_joint * [1.0; 1.0; 1.0]))
