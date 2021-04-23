"""
    InvertedPendulum
"""
mutable struct InvertedPendulum{T} <: ContactDynamicsModel
    dim::Dimensions

    mb::T # mass
	m1::T
	m2::T
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
function _kinematics(model::InvertedPendulum, q; mode = :com)
	θ, d1, d2 = q

	if mode == :d1
		return [model.l * sin(θ) + d1 * cos(θ);
				model.l * cos(θ) - d1 * sin(θ)]
	elseif mode == :d2
		return [model.l * sin(θ) + d2 * cos(θ);
				model.l * cos(θ) - d2 * sin(θ)]
	elseif mode == :com
		return [-0.5 * model.l * sin(θ);
				 0.5 * model.l * cos(θ)]
	else
		@error "incorrect mode"
	 	return zeros(2)
	end
end

function _jacobian(model::InvertedPendulum, q; mode = :com)
	θ, d1, d2 = q

	if mode == :d1
 		return [(-model.l * cos(θ) - d1 * sin(θ)) cos(θ) 0.0;
		        (-model.l * sin(θ) - d1 * cos(θ)) -sin(θ) 0.0]
	elseif mode == :d2
 		return [(-model.l * cos(θ) - d2 * sin(θ)) 0.0 cos(θ);
		        (-model.l * sin(θ) - d2 * cos(θ)) 0.0 -sin(θ)]
 	elseif mode == :com
 		return [-0.5 * model.l * cos(θ) 0.0 0.0;
		        -0.5 * model.l * sin(θ) 0.0 0.0]
 	else
 		@error "incorrect mode"
 	 	return zeros(2)
 	end
end

function kinematics(model::InvertedPendulum, q)
	[_kinematics(model, q, mode = :d1);
	 _kinematics(model, q, mode = :d2)]
end

function lagrangian(model::InvertedPendulum, q, q̇)
	L = 0.0

	vθ = _jacobian(model, q, mode = :com) * q̇[1]
	L += 0.5 * model.m * transpose(vθ) * vθ
	L -= model.mb * model.g * _kinematics(model, q, mode = :com)[2]

	vd1 = _jacobian(model, q, mode = :d1) * q̇[1]
	L += 0.5 * model.m1 * transpose(vd1) * vd1
	L -= model.m1 * model.g * _kinematics(model, q, mode = :d1)[2]

	vd2 = _jacobian(model, q, mode = :d2) * q̇[1]
	L += 0.5 * model.m2 * transpose(vd2) * vd2
	L -= model.m2 * model.g * _kinematics(model, q, mode = :d2)[2]

	return L
end

function M_func(model::InvertedPendulum, q)
	Jθ = _jacobian(model, q, mode = :com)
	Jd1 = _jacobian(model, q, model = :d1)
	Jd2 = _jacobian(model, q, model = :d2)

	return model.mb * transpose(Jθ) * Jθ + model.m1 * transpose(Jd1) * transpose(Jd2)
end

function ϕ_func(model::InvertedPendulum, q)
	# walls at x = -0.5, x = 0.5
    SVector{2}([_kinematics(model, q, mode = :d2)[1] + 0.5;
	            0.5 - _kinematics(model, q, mode = :d1)[1]])
end

function rot2D(x)
	[cos(x) sin(x);
    -sin(x) cos(x)]
end

function J_func(model::InvertedPendulum, q)
	r1 = rot2D(-0.5 * π)
	r2 = rot2D(0.5 * π)

    SMatrix{4, 1}([r1 * _jacobian(model, q, mode = :d1);
		           r2 * _jacobian(model, q, mode = :d2)])
end

# function control_switch(model, x, ϵ = 1.0e-6)
# 	ϕ = ϕ_func(model, x)
# 	IfElse.ifelse(ϕ[1] > 1.0e-3, IfElse.ifelse(ϕ[2] > 1.0e-3, 0.0, 1.0), 1.0)
# end

function B_func(model::InvertedPendulum, q)
	# @SMatrix [control_switch(model, q) * model.l * cos(q[1])]
	@SMatrix [model.l 1.0 0.0;
	 		  model.l 0.0 1.0]
end

function A_func(::InvertedPendulum, q)
	@SMatrix [1.0 0.0 0.0;
			  0.0 1.0 0.0;
			  0.0 0.0 1.0]
end

# Parameters
g = 9.81 # gravity
μ_world = 1.0  # coefficient of friction
μ_joint = 0.25

mb = 1.0 # body mass
m1 = 0.1 # mass of ee
m2 = 0.1 # mass of ee
l = 1.0  # leg mass

# Dimensions
nq = 3
nu = 2
nw = 3
nc = 2
nf = 2
nb = nc * nf

inverted_pendulum = InvertedPendulum(Dimensions(nq, nu, nw, nc, nb),
					   mb, m1, m2, l,
					   μ_world, μ_joint, g,
					   BaseMethods(), DynamicsMethods(), ContactMethods(),
					   ResidualMethods(), ResidualMethods(),
					   SparseStructure(spzeros(0, 0), spzeros(0, 0)),
					   SVector{3}(μ_joint * [1.0; 0.0; 0.0]),
					   environment_2D_flat())
