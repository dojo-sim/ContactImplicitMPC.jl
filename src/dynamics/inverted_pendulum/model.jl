"""
    InvertedPendulum
"""
mutable struct InvertedPendulum{T} <: ContactDynamicsModel
    dim::Dimensions

    m::T # mass
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
function _kinematics(model::InvertedPendulum, q; mode = :ee)
	if mode == :ee
		return [-model.l * sin(q[1]);
		         model.l * cos(q[1])]
	elseif mode == :com
		return [-0.5 * model.l * sin(q[1]);
				 0.5 * model.l * cos(q[1])]
	else
		@error "incorrect mode"
	 	return zeros(2)
	end
end

function _jacobian(model::InvertedPendulum, q; mode = :ee)
	[-model.l * cos(q[1]);
	 -model.l * sin(q[1])]

	 if mode == :ee
 		return [-model.l * cos(q[1]);
		        -model.l * sin(q[1])]
 	elseif mode == :com
 		return [-0.5 * model.l * cos(q[1]);
		        -0.5 * model.l * sin(q[1])]
 	else
 		@error "incorrect mode"
 	 	return zeros(2)
 	end
end

function kinematics(model::InvertedPendulum, q; mode = :ee)
	[_kinematics(model, q, mode = mode);
	 _kinematics(model, q, mode = mode)]
end

function lagrangian(model::InvertedPendulum, q, q̇)
	L = 0.0

	v = _jacobian(model, q, mode = :com) * q̇[1]
	L += 0.5 * model.m * transpose(v) * v
	L -= model.m * model.g * _kinematics(model, q, mode = :com)[2]

	return L
end

function M_func(model::InvertedPendulum, q)
	J = _jacobian(model, q, mode = :com)
	return [model.m * transpose(J) * J]
end

function ϕ_func(model::InvertedPendulum, q)
	# walls at x = -0.5, x = 0.5
    SVector{2}([_kinematics(model, q, mode = :ee)[1] + 0.5;
	            0.5 - _kinematics(model, q, mode = :ee)[1]])
end

function rot2D(x)
	[cos(x) sin(x);
    -sin(x) cos(x)]
end

function J_func(model::InvertedPendulum, q)
	r1 = rot2D(-0.5 * π)
	r2 = rot2D(0.5 * π)

    SMatrix{4, 1}([r1 * _jacobian(model, q, mode = :ee);
		           r2 * _jacobian(model, q, mode = :ee)])
end

function B_func(model::InvertedPendulum, q)
	@SMatrix [model.l * cos(q[1])]
end

function A_func(::InvertedPendulum, q)
	@SMatrix [1.0]
end

M_func(model, [1.0])

# Parameters
g = 9.81 # gravity
μ_world = 1.0  # coefficient of friction
μ_joint = 0.1

m = 1.0 # body mass
l = 1.0  # leg mass

# Dimensions
nq = 1
nu = 1
nw = 1
nc = 2
nf = 2
nb = nc * nf

inverted_pendulum = InvertedPendulum(Dimensions(nq, nu, nw, nc, nb),
					   m, l,
					   μ_world, μ_joint, g,
					   BaseMethods(), DynamicsMethods(), ContactMethods(),
					   ResidualMethods(), ResidualMethods(),
					   SparseStructure(spzeros(0, 0), spzeros(0, 0)),
					   SVector{1}(μ_joint * ones(1)),
					   environment_2D_flat())
