"""
    PushBot
"""
mutable struct PushBot{T} <: ContactDynamicsModel
    dim::Dimensions

    mb::T # mass
	ma::T
	l::T # length

    μ_world::T  # coefficient of friction
    μ_joint::T  # gravity
	g::T

	base::BaseMethods
	dyn::DynamicsMethods
	con::ContactMethods
	res::ResidualMethods
	linearized::ResidualMethods

	spa::SparseStructure

	joint_friction::SVector

	env::Environment
end


# Kinematics
function _kinematics(model::PushBot, q; mode = :com)
	θ, d  = q
	if mode == :d
		return [-model.l * sin(θ) + d * cos(θ);
				 model.l * cos(θ) + d * sin(θ)]
	elseif mode == :com
		return [-1.0 * model.l * sin(θ);
				 1.0 * model.l * cos(θ)]
	elseif mode == :ee
		return [-model.l * sin(θ);
				 model.l * cos(θ)]
	else
		@error "incorrect mode"
	 	return zeros(2)
	end
end

function _jacobian(model::PushBot, q; mode = :com)
	θ, d = q

	if mode == :d
 		return [(-model.l * cos(θ) - d * sin(θ)) cos(θ);
		        (-model.l * sin(θ) + d * cos(θ)) sin(θ)]
 	elseif mode == :com
 		return [-1.0 * model.l * cos(θ) 0.0;
		        -1.0 * model.l * sin(θ) 0.0]
	elseif mode == :ee
		return [-model.l * cos(θ) 0.0;
		        -model.l * sin(θ) 0.0]
 	else
 		@error "incorrect mode"
 	 	return
 	end
end

function kinematics(model::PushBot, q)
	[_kinematics(model, q, mode = :d);
	 _kinematics(model, q, mode = :d)]
end

function lagrangian(model::PushBot, q, q̇)
	L = 0.0

	vθ = _jacobian(model, q, mode = :com) * q̇
	L += 0.5 * model.mb * transpose(vθ) * vθ
	L -= model.mb * model.g * _kinematics(model, q, mode = :com)[2]

	vd1 = _jacobian(model, q, mode = :d) * q̇
	L += 0.5 * model.ma * transpose(vd1) * vd1
	L -= model.ma * model.g * _kinematics(model, q, mode = :d)[2]

	return L
end

function M_func(model::PushBot, q)
	Jθ = _jacobian(model, q, mode = :com)
	Jd = _jacobian(model, q, mode = :d)

	return model.mb * transpose(Jθ) * Jθ + model.ma * transpose(Jd) * Jd
end

function ϕ_func(model::PushBot, q)
	# walls at x = -0.5, x = 0.5
    SVector{2}([_kinematics(model, q, mode = :d)[1] + 0.5;
	            0.5 - _kinematics(model, q, mode = :d)[1]])
end

function rot2D(x)
	[cos(x) sin(x);
    -sin(x) cos(x)]
end

function J_func(model::PushBot, q)
	r1 = [0.0 -1.0; 1.0 0.0]
	r2 = [0.0 1.0; -1.0 0.0]

    SMatrix{4, 2}([r1 * _jacobian(model, q, mode = :d);
		           r2 * _jacobian(model, q, mode = :d)])
end

# function control_switch(model, x, ϵ = 1.0e-6)
# 	ϕ = ϕ_func(model, x)
# 	IfElse.ifelse(ϕ[1] > 1.0e-3, IfElse.ifelse(ϕ[2] > 1.0e-3, 0.0, 1.0), 1.0)
# end

function B_func(model::PushBot, q)
	@SMatrix [1.0 * model.l 1.0;
	          1.0 1.0 / model.l]
end

function A_func(::PushBot, q)
	@SMatrix [1.0 0.0;
			  0.0 1.0]
end

# Parameters
g = 9.81 # gravity
μ_world = 0.5  # coefficient of friction
μ_joint = 10.0

mb = 1.0 # body mass
# ma = 0.1 # mass of ee
ma = 0.01 # mass of ee # CHANGED
l = 1.0  # leg mass

# Dimensions
nq = 2
nu = 2
nw = nq
nc = 2
nf = 2
nb = nc * nf

pushbot = PushBot(Dimensions(nq, nu, nw, nc, nb),
					   mb, ma, l,
					   μ_world, μ_joint, g,
					   BaseMethods(), DynamicsMethods(), ContactMethods(),
					   ResidualMethods(), ResidualMethods(),
					   SparseStructure(spzeros(0, 0), spzeros(0, 0)),
					   SVector{2}(μ_joint * [1.0; 1.0]),
					   environment_2D_flat())
