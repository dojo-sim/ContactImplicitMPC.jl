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

function control_switch(model, x, ϵ = 1.0e-6)
	ϕ = ϕ_func(model, x)
	IfElse.ifelse(ϕ[1] > 1.0e-3, IfElse.ifelse(ϕ[2] > 1.0e-3, 0.0, 1.0), 1.0)
end

function B_func(model::InvertedPendulum, q)
	# @SMatrix [control_switch(model, q) * model.l * cos(q[1])]
	@SMatrix [model.l * cos(q[1])]
end

function A_func(::InvertedPendulum, q)
	@SMatrix [1.0]
end

function dynamics(model::ContactDynamicsModel, h, q0, q1, u1, w1, λ1, q2)

	# evalutate at midpoint
	qm1 = 0.5 * (q0 + q1)
    vm1 = (q1 - q0) / h[1]
    qm2 = 0.5 * (q1 + q2)
    vm2 = (q2 - q1) / h[1]

	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	return (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
		+ transpose(B_func(model, qm2)) * u1 #NOTE: swap this to B_func for control flow
		+ transpose(A_fast(model, qm2)) * w1
		+ transpose(J_fast(model, q2)) * λ1
		- h[1] * model.joint_friction .* vm2)
end

function residual(model::InvertedPendulum, z, θ, κ)
	nc = model.dim.c
	nb = model.dim.b
	nf = Int(nb / nc)
	np = dim(model.env)

	q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
	q2, γ1, b1, ψ1, η1, s1, s2 = unpack_z(model, z)

	ϕ = ϕ_func(model, q2)

	k = kinematics(model, q2)
	λ1 = contact_forces(model, γ1, b1, q2, k)
	vT_stack = velocity_stack(model, q1, q2, k, h)
	ψ_stack = transpose(E_func(model)) * ψ1

	[dynamics(model, h, q0, q1, u1, w1, λ1, q2); # NOTE swap to dynamics for control flow
	 s1 - ϕ;
	 vT_stack + ψ_stack - η1;
	 s2 .- (μ[1] * γ1 .- E_func(model) * b1);
	 γ1 .* s1 .- κ;
	 b1 .* η1 .- κ;
	 ψ1 .* s2 .- κ]
end

# Parameters
g = 9.81 # gravity
μ_world = 1.0  # coefficient of friction
μ_joint = 0.25

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
