mutable struct QuadrupedLinear12{T} <: ContactDynamicsModel
	dim::Dimensions

	g::T
	μ_world::T
	μ_joint::T

	mb::T
	mf::T
	Ix::T
	Iy::T
	Iz::T
	l_torso::T
	w_torso::T

    # # torso
    # l_torso::T
    # d_torso::T
    # m_torso::T
    # J_torso::T
	#
    # # leg 1
    #     # thigh
    # l_thigh1::T
    # d_thigh1::T
    # m_thigh1::T
    # J_thigh1::T
    #     # calf
    # l_calf1::T
    # d_calf1::T
    # m_calf1::T
    # J_calf1::T
	#
    # # leg 2
    #     # thigh
    # l_thigh2::T
    # d_thigh2::T
    # m_thigh2::T
    # J_thigh2::T
    #     # calf
    # l_calf2::T
    # d_calf2::T
    # m_calf2::T
    # J_calf2::T
	#
	# # leg 3
    #     # thigh
    # l_thigh3::T
    # d_thigh3::T
    # m_thigh3::T
    # J_thigh3::T
    #     # calf
    # l_calf3::T
    # d_calf3::T
    # m_calf3::T
    # J_calf3::T
	#
	# # leg 4
    #     # thigh
    # l_thigh4::T
    # d_thigh4::T
    # m_thigh4::T
    # J_thigh4::T
    #     # calf
    # l_calf4::T
    # d_calf4::T
    # m_calf4::T
    # J_calf4::T

	alt

	# fast methods
	base
	dyn
	con
	res
	linearized

	spa::SparseStructure

	joint_friction

	env::Environment
end

# Trunk model
#                             middle   com
#                                267 mm
#             o---------------------|---X-----------------o
#                                   13mm
# BRhip o___________________________|___________________________o FRhip
#                    183 mm                    183mm
# kinematics
# function kinematics_1(model::QuadrupedLinear12, q; body = :torso, mode = :ee)
# 	x = q[1]
# 	z = q[2]
#
# 	if body == :torso
# 		l = model.l_torso
# 		d = model.d_torso
# 		θ = q[3]
# 	elseif body == :thigh_1
# 		l = model.l_thigh1
# 		d = model.d_thigh1
# 		θ = q[4]
# 	elseif body == :thigh_2
# 		l = model.l_thigh2
# 		d = model.d_thigh2
# 		θ = q[6]
# 	else
# 		@error "incorrect body specification"
# 	end
#
# 	if mode == :ee
# 		return [x + l * sin(θ); z - l * cos(θ)]
# 	elseif mode == :com
# 		return [x + d * sin(θ); z - d * cos(θ)]
# 	else
# 		@error "incorrect mode specification"
# 	end
# end
#
# function jacobian_1(model::QuadrupedLinear12, q; body = :torso, mode = :ee)
# 	jac = zeros(eltype(q), 2, model.dim.q)
# 	jac[1, 1] = 1.0
# 	jac[2, 2] = 1.0
# 	if body == :torso
# 		r = mode == :ee ? model.l_torso : model.d_torso
# 		θ = q[3]
# 		jac[1, 3] = r * cos(θ)
# 		jac[2, 3] = r * sin(θ)
# 	elseif body == :thigh_1
# 		r = mode == :ee ? model.l_thigh1 : model.d_thigh1
# 		θ = q[4]
# 		jac[1, 4] = r * cos(θ)
# 		jac[2, 4] = r * sin(θ)
# 	elseif body == :thigh_2
# 		r = mode == :ee ? model.l_thigh2 : model.d_thigh2
# 		θ = q[6]
# 		jac[1, 6] = r * cos(θ)
# 		jac[2, 6] = r * sin(θ)
# 	else
# 		@error "incorrect body specification"
# 	end
#
# 	return jac
# end
#
# function kinematics_2(model::QuadrupedLinear12, q; body = :calf_1, mode = :ee)
#
# 	if body == :calf_1
# 		p = kinematics_1(model, q, body = :thigh_1, mode = :ee)
# 		θb = q[5]
#
# 		lb = model.l_calf1
# 		db = model.d_calf1
# 	elseif body == :calf_2
# 		p = kinematics_1(model, q, body = :thigh_2, mode = :ee)
# 		θb = q[7]
#
# 		lb = model.l_calf2
# 		db = model.d_calf2
# 	elseif body == :thigh_3
# 		p = kinematics_1(model, q, body = :torso, mode = :ee)
# 		θb = q[8]
#
# 		lb = model.l_thigh3
# 		db = model.d_thigh3
# 	elseif body == :thigh_4
# 		p = kinematics_1(model, q, body = :torso, mode = :ee)
# 		θb = q[10]
#
# 		lb = model.l_thigh4
# 		db = model.d_thigh4
# 	else
# 		@error "incorrect body specification"
# 	end
#
# 	if mode == :ee
# 		return p + [lb * sin(θb); -1.0 * lb * cos(θb)]
# 	elseif mode == :com
# 		return p + [db * sin(θb); -1.0 * db * cos(θb)]
# 	else
# 		@error "incorrect mode specification"
# 	end
# end
#
# function jacobian_2(model::QuadrupedLinear12, q; body = :calf_1, mode = :ee)
#
# 	if body == :calf_1
# 		jac = jacobian_1(model, q, body = :thigh_1, mode = :ee)
#
# 		θb = q[5]
#
# 		r = mode == :ee ? model.l_calf1 : model.d_calf1
#
# 		jac[1, 5] += r * cos(θb)
# 		jac[2, 5] += r * sin(θb)
# 	elseif body == :calf_2
# 		jac = jacobian_1(model, q, body = :thigh_2, mode = :ee)
#
# 		θb = q[7]
#
# 		r = mode == :ee ? model.l_calf2 : model.d_calf2
#
# 		jac[1, 7] += r * cos(θb)
# 		jac[2, 7] += r * sin(θb)
# 	elseif body == :thigh_3
# 		jac = jacobian_1(model, q, body = :torso, mode = :ee)
# 		θb = q[8]
#
# 		r = mode == :ee ? model.l_thigh3 : model.d_thigh3
#
# 		jac[1, 8] += r * cos(θb)
# 		jac[2, 8] += r * sin(θb)
# 	elseif body == :thigh_4
# 		jac = jacobian_1(model, q, body = :torso, mode = :ee)
#
# 		θb = q[10]
#
# 		r = mode == :ee ? model.l_thigh4 : model.d_thigh4
#
# 		jac[1, 10] += r * cos(θb)
# 		jac[2, 10] += r * sin(θb)
# 	else
# 		@error "incorrect body specification"
# 	end
#
# 	return jac
# end
#
# function kinematics_3(model::QuadrupedLinear12, q; body = :calf_3, mode = :ee)
#
# 	if body == :calf_3
# 		p = kinematics_2(model, q, body = :thigh_3, mode = :ee)
# 		θc = q[9]
#
#
# 		lb = model.l_calf3
# 		db = model.d_calf3
# 	elseif body == :calf_4
# 		p = kinematics_2(model, q, body = :thigh_4, mode = :ee)
#
# 		θc = q[11]
#
# 		lb = model.l_calf4
# 		db = model.d_calf4
# 	else
# 		@error "incorrect body specification"
# 	end
#
# 	if mode == :ee
# 		return p + [lb * sin(θc); -1.0 * lb * cos(θc)]
# 	elseif mode == :com
# 		return p + [db * sin(θc); -1.0 * db * cos(θc)]
# 	else
# 		@error "incorrect mode specification"
# 	end
# end
#
# function jacobian_3(model::QuadrupedLinear12, q; body = :calf_3, mode = :ee)
#
# 	if body == :calf_3
# 		jac = jacobian_2(model, q, body = :thigh_3, mode = :ee)
#
# 		θc = q[9]
#
# 		r = mode == :ee ? model.l_calf3 : model.d_calf3
#
# 		jac[1, 9] += r * cos(θc)
# 		jac[2, 9] += r * sin(θc)
#
# 	elseif body == :calf_4
# 		jac = jacobian_2(model, q, body = :thigh_4, mode = :ee)
#
# 		θc = q[11]
#
# 		r = mode == :ee ? model.l_calf4 : model.d_calf4
#
# 		jac[1, 11] += r * cos(θc)
# 		jac[2, 11] += r * sin(θc)
# 	else
# 		@error "incorrect body specification"
# 	end
#
# 	return jac
# end

# Lagrangian
function lagrangian(model::QuadrupedLinear12, q, q̇)
	return 0.0
	# L = 0.0
	#
	# # # torso
	# p_torso = kinematics_1(model, q, body = :torso, mode = :com)
	# J_torso = jacobian_1(model, q, body = :torso, mode = :com)
	# v_torso = J_torso * q̇
	# #
	# L += 0.5 * model.m_torso * transpose(v_torso) * v_torso
	# L += 0.5 * model.J_torso * q̇[3]^2.0
	# L -= model.m_torso * model.g * p_torso[2]
	#
	# # thigh 1
	# p_thigh_1 = kinematics_1(model, q, body = :thigh_1, mode = :com)
	# J_thigh_1 = jacobian_1(model, q, body = :thigh_1, mode = :com)
	# v_thigh_1 = J_thigh_1 * q̇
	#
	# L += 0.5 * model.m_thigh1 * transpose(v_thigh_1) * v_thigh_1
	# L += 0.5 * model.J_thigh1 * q̇[4]^2.0
	# L -= model.m_thigh1 * model.g * p_thigh_1[2]
	#
	# # leg 1
	# p_calf_1 = kinematics_2(model, q, body = :calf_1, mode = :com)
	# J_calf_1 = jacobian_2(model, q, body = :calf_1, mode = :com)
	# v_calf_1 = J_calf_1 * q̇
	#
	# L += 0.5 * model.m_calf1 * transpose(v_calf_1) * v_calf_1
	# L += 0.5 * model.J_calf1 * q̇[5]^2.0
	# L -= model.m_calf1 * model.g * p_calf_1[2]
	#
	# # thigh 2
	# p_thigh_2 = kinematics_1(model, q, body = :thigh_2, mode = :com)
	# J_thigh_2 = jacobian_1(model, q, body = :thigh_2, mode = :com)
	# v_thigh_2 = J_thigh_2 * q̇
	#
	# L += 0.5 * model.m_thigh2 * transpose(v_thigh_2) * v_thigh_2
	# L += 0.5 * model.J_thigh2 * q̇[6]^2.0
	# L -= model.m_thigh2 * model.g * p_thigh_2[2]
	#
	# # leg 2
	# p_calf_2 = kinematics_2(model, q, body = :calf_2, mode = :com)
	# J_calf_2 = jacobian_2(model, q, body = :calf_2, mode = :com)
	# v_calf_2 = J_calf_2 * q̇
	#
	# L += 0.5 * model.m_calf2 * transpose(v_calf_2) * v_calf_2
	# L += 0.5 * model.J_calf2 * q̇[7]^2.0
	# L -= model.m_calf2 * model.g * p_calf_2[2]
	#
	# # thigh 3
	# p_thigh_3 = kinematics_2(model, q, body = :thigh_3, mode = :com)
	# J_thigh_3 = jacobian_2(model, q, body = :thigh_3, mode = :com)
	# v_thigh_3 = J_thigh_3 * q̇
	#
	# L += 0.5 * model.m_thigh3 * transpose(v_thigh_3) * v_thigh_3
	# L += 0.5 * model.J_thigh3 * q̇[8]^2.0
	# L -= model.m_thigh3 * model.g * p_thigh_3[2]
	#
	# # leg 3
	# p_calf_3 = kinematics_3(model, q, body = :calf_3, mode = :com)
	# J_calf_3 = jacobian_3(model, q, body = :calf_3, mode = :com)
	# v_calf_3 = J_calf_3 * q̇
	#
	# L += 0.5 * model.m_calf3 * transpose(v_calf_3) * v_calf_3
	# L += 0.5 * model.J_calf3 * q̇[9]^2.0
	# L -= model.m_calf3 * model.g * p_calf_3[2]
	#
	# # thigh 4
	# p_thigh_4 = kinematics_2(model, q, body = :thigh_4, mode = :com)
	# J_thigh_4 = jacobian_2(model, q, body = :thigh_4, mode = :com)
	# v_thigh_4 = J_thigh_4 * q̇
	#
	# L += 0.5 * model.m_thigh4 * transpose(v_thigh_4) * v_thigh_4
	# L += 0.5 * model.J_thigh4 * q̇[10]^2.0
	# L -= model.m_thigh4 * model.g * p_thigh_4[2]
	#
	# # leg 4
	# p_calf_4 = kinematics_3(model, q, body = :calf_4, mode = :com)
	# J_calf_4 = jacobian_3(model, q, body = :calf_4, mode = :com)
	# v_calf_4 = J_calf_4 * q̇
	#
	# L += 0.5 * model.m_calf4 * transpose(v_calf_4) * v_calf_4
	# L += 0.5 * model.J_calf4 * q̇[11]^2.0
	# L -= model.m_calf4 * model.g * p_calf_4[2]
	#
	# return L
end

function _dLdq(model::QuadrupedLinear12, q, q̇)
	Lq(x) = lagrangian(model, x, q̇)
	ForwardDiff.gradient(Lq, q)
end

function _dLdq̇(model::QuadrupedLinear12, q, q̇)
	Lq̇(x) = lagrangian(model, q, x)
	ForwardDiff.gradient(Lq̇, q̇)
end

function kinematics(model::QuadrupedLinear12, q)
	pw1 = q[6 .+ (1:3)]
	pw2 = q[9 .+ (1:3)]
	pw3 = q[12 .+ (1:3)]
	pw4 = q[15 .+ (1:3)]
	SVector{12}([pw1; pw2; pw3; pw4])
end

# Methods
function M_func(model::QuadrupedLinear12, q)
	Diagonal(@SVector [model.mb,
					   model.mb,
					   model.mb,
					   model.Ix, model.Iy, model.Iz,
					   model.mf, model.mf, model.mf,
					   model.mf, model.mf, model.mf,
					   model.mf, model.mf, model.mf,
					   model.mf, model.mf, model.mf])
end

function ϕ_func(model::QuadrupedLinear12, q)
	k = kinematics(model, q)
	@SVector [k[3], k[6], k[9], k[12]]
	# p_calf_1 = kinematics_2(model, q, body = :calf_1, mode = :ee)
	# p_calf_2 = kinematics_2(model, q, body = :calf_2, mode = :ee)
	# p_calf_3 = kinematics_3(model, q, body = :calf_3, mode = :ee)
	# p_calf_4 = kinematics_3(model, q, body = :calf_4, mode = :ee)
	#
	# SVector{4}([p_calf_1[2] - model.env.surf(p_calf_1[1:1]);
	# 			p_calf_2[2] - model.env.surf(p_calf_2[1:1]);
	# 			p_calf_3[2] - model.env.surf(p_calf_3[1:1]);
	# 			p_calf_4[2] - model.env.surf(p_calf_4[1:1])])
end

function B_func(model::QuadrupedLinear12, q)
	p_torso = q[1:3]
	rot = MRP(q[4:6]...)

	# r in world frame
	r1 = q[6 .+ (1:3)] - p_torso
	r2 = q[9 .+ (1:3)] - p_torso
	r3 = q[12 .+ (1:3)] - p_torso
	r4 = q[15 .+ (1:3)] - p_torso

	z3 = zeros(3, 3)

	SMatrix{18, 12}([I I I I;
	                transpose(rot) * skew(r1) transpose(rot) * skew(r2) transpose(rot) * skew(r3) transpose(rot) * skew(r4);
					-I z3 z3 z3;
					z3 -I z3 z3;
					z3 z3 -I z3;
					z3 z3 z3 -I])
	# @SMatrix [0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
	# 		  0.0  0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0;
	# 		  0.0  0.0 -1.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0;
	# 		  0.0  0.0  0.0  0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0;
	# 		  0.0  0.0 -1.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0;
	# 		  0.0  0.0  0.0  0.0  0.0  0.0  0.0 -1.0  1.0  0.0  0.0;
	# 		  0.0  0.0 -1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0;
	# 		  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 -1.0  1.0]
end

function A_func(model::QuadrupedLinear12, q)
	@SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
end

function J_func(model::QuadrupedLinear12, q)
	p_torso = q[1:3]
	rot = MRP(q[4:6]...)

	# r in world frame
	r1 = q[6 .+ (1:3)] - p_torso
	r2 = q[9 .+ (1:3)] - p_torso
	r3 = q[12 .+ (1:3)] - p_torso
	r4 = q[15 .+ (1:3)] - p_torso

	z3 = zeros(3, 3)
	J = transpose([
		 zeros(6, 12);
	     I z3 z3 z3;
	     z3 I z3 z3;
	     z3 z3 I z3;
	     z3 z3 z3 I])
	# idx = collect([3, 6, 9, 12])
	# J_calf_1 = jacobian_2(model, q, body = :calf_1, mode = :ee)
	# J_calf_2 = jacobian_2(model, q, body = :calf_2, mode = :ee)
	# J_calf_3 = jacobian_3(model, q, body = :calf_3, mode = :ee)
	# J_calf_4 = jacobian_3(model, q, body = :calf_4, mode = :ee)
	# # map = [1.0; -1.0]
	#
	# return [J_calf_1;
	# 		J_calf_2;
	# 		J_calf_3;
	# 		J_calf_4]
end

function C_func(model::QuadrupedLinear12, q, q̇)
	SVector{18}([0.0, 0.0, model.g * model.mb,
				 cross(q̇[4:6], Diagonal([model.Ix, model.Iy, model.Iz]) * q̇[4:6])...,
				 0.0, 0.0, model.g * model.mf,
				 0.0, 0.0, model.g * model.mf,
				 0.0, 0.0, model.g * model.mf,
				 0.0, 0.0, model.g * model.mf])
	# tmp_q(z) = _dLdq̇(model, z, q̇)
	# tmp_q̇(z) = _dLdq̇(model, q, z)
	#
	# ForwardDiff.jacobian(tmp_q, q) * q̇ - _dLdq(model, q, q̇)
end

function contact_forces(model::QuadrupedLinear12, γ1, b1, q2, k)
	k = kinematics(model, q2)
	m = friction_mapping(model.env)

	SVector{12}([transpose(rotation(model.env, k[1:2]))  * [m * b1[1:4];   γ1[1]];
				transpose(rotation(model.env, k[4:5]))   * [m * b1[5:8];   γ1[2]];
				transpose(rotation(model.env, k[7:8]))   * [m * b1[9:12];  γ1[3]];
				transpose(rotation(model.env, k[10:11])) * [m * b1[13:16]; γ1[4]]])
end

function velocity_stack(model::QuadrupedLinear12, q1, q2, k, h)
	k = kinematics(model, q2)
	v = J_func(model, q2) * (q2 - q1) / h[1]

	v1_surf = rotation(model.env, k[1:2]) * v[1:3]
	v2_surf = rotation(model.env, k[4:5]) * v[4:6]
	v3_surf = rotation(model.env, k[7:8]) * v[7:9]
	v4_surf = rotation(model.env, k[10:11]) * v[10:12]

	# v1_surf = v[1:2]
	# v2_surf = v[3:4]
	# v3_surf = v[5:6]
	# v4_surf = v[7:8]

	SVector{16}([friction_mapping(model.env)' * v1_surf[1:2];
				friction_mapping(model.env)' * v2_surf[1:2];
				friction_mapping(model.env)' * v3_surf[1:2];
				friction_mapping(model.env)' * v4_surf[1:2]])
end

################################################################################
# Instantiation
################################################################################
# Dimensions
nq = 3 + 3 + 3 * 4        # configuration dimension
nu = 4 * 3                # control dimension
nc = 4                    # number of contact points
nf = 4                    # number of parameters for friction cone
nb = nc * nf
nw = 3

# World parameters
g = 9.81      # gravity
μ_world = 0.5 # coefficient of friction
μ_joint = 0.0 # coefficient of torque friction at the joints

# ~Unitree A1
mf = 0.01

# TRUCK ONLY TODO: parallel axis theorem to add shoulders
mb = 4.713
Ix = 0.01683993
Iy = 0.056579028
Iz = 0.064713601

l_torso = 0.5 * 0.267 # dimension from com
w_torso = 0.5 * 0.194 # dimension from com

# # ~Unitree A1
# # Model parameters
# m_torso = 4.713 + 4 * 0.696
# m_thigh = 1.013
# m_leg = 0.166
#
# J_torso = 0.01683 + 4 * 0.696 * 0.183^2.0
# J_thigh = 0.00552
# J_leg = 0.00299
#
# l_torso = 0.183*2
# l_thigh = 0.2
# l_leg = 0.2
#
# d_torso = 0.5 * l_torso + 0.0127
# d_thigh = 0.5 * l_thigh - 0.00323
# d_leg = 0.5 * l_leg - 0.006435
#
# m_payload = 5.0
# J_payload = 0.05

quadrupedlinear = QuadrupedLinear12(Dimensions(nq, nu, nw, nc, nb),
				g, μ_world, μ_joint,
				mb, mf, Ix, Iy, Iz, l_torso, w_torso,
				# l_torso, d_torso, m_torso, J_torso,
				# l_thigh, d_thigh, m_thigh, J_thigh,
				# l_leg, d_leg, m_leg, J_leg,
				# l_thigh, d_thigh, m_thigh, J_thigh,
				# l_leg, d_leg, m_leg, J_leg,
				# l_thigh, d_thigh, m_thigh, J_thigh,
				# l_leg, d_leg, m_leg, J_leg,
				# l_thigh, d_thigh, m_thigh, J_thigh,
				# l_leg, d_leg, m_leg, J_leg,
				zeros(nc),
				BaseMethods(), DynamicsMethods(), ContactMethods(),
				ResidualMethods(), ResidualMethods(),
				SparseStructure(spzeros(0, 0), spzeros(0, 0)),
				SVector{nq}([zeros(3); μ_joint * ones(nq - 3)]),
				environment_3D_flat())
