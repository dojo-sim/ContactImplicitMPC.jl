mutable struct Quadruped3DAlt{T} <: ContactModel
	dim::Dimensions

	g::T
	μ_world::T
	μ_joint::T

    # torso
    l_torso::T
	d_torso::T
    w_torso::T
    m_torso::T
    J_torso

    # leg 1
        # thigh
    l_thigh1::T
    d_thigh1::T
    m_thigh1::T
    J_thigh1::T
        # calf
    l_calf1::T
    d_calf1::T
    m_calf1::T
    J_calf1::T

    # leg 2
        # thigh
    l_thigh2::T
    d_thigh2::T
    m_thigh2::T
    J_thigh2::T
        # calf
    l_calf2::T
    d_calf2::T
    m_calf2::T
    J_calf2::T

	# leg 3
        # thigh
    l_thigh3::T
    d_thigh3::T
    m_thigh3::T
    J_thigh3::T
        # calf
    l_calf3::T
    d_calf3::T
    m_calf3::T
    J_calf3::T

	# leg 4
        # thigh
    l_thigh4::T
    d_thigh4::T
    m_thigh4::T
    J_thigh4::T
        # calf
    l_calf4::T
    d_calf4::T
    m_calf4::T
    J_calf4::T

	# fast methods
	base
	dyn

	joint_friction
end

# Trunk model
#                             middle   com
#                                267 mm
#             o---------------------|---X-----------------o
#                                   13mm
# BRhip o___________________________|___________________________o FRhip
#                    183 mm                    183mm
# kinematics
function kinematics(model::Quadruped3DAlt, q; body = :leg1, mode = :ee)
	x, y, z = q[1:3]
	R = euler_rotation_matrix(q[4:6])

	if body == :torso
		return R[:, 1] * model.d_torso + q[1:3]
	end

	if body == :leg1
		l_torso = model.l_torso
		w_torso = model.w_torso
		θ = q[6 .+ (1:3)]

		if mode == :ee
			px = model.l_thigh1 * sin(θ[1]) + model.l_calf1 * sin(θ[2])
			pz = -model.l_thigh1 * cos(θ[1]) - model.l_calf1 * cos(θ[2])
		elseif mode == :calf
			px = model.l_thigh1 * sin(θ[1]) + model.d_calf1 * sin(θ[2])
			pz = -model.l_thigh1 * cos(θ[1]) - model.d_calf1 * cos(θ[2])
		elseif mode == :thigh
			px = model.d_thigh1 * sin(θ[1])
			pz = -model.d_thigh1 * cos(θ[1])
		else
			@error "incorrect mode specification"
		end

		p̄x = px + l_torso
		p̄y = -1.0 * sin(θ[3]) * pz + w_torso
		p̄z = cos(θ[3]) * pz

		return R * [p̄x; p̄y; p̄z] + q[1:3]
	end

	if body == :leg2
		l_torso = model.l_torso
		w_torso = -model.w_torso
		θ = q[6 + 3 .+ (1:3)]

		if mode == :ee
			px = model.l_thigh2 * sin(θ[1]) + model.l_calf2 * sin(θ[2])
			pz = -model.l_thigh2 * cos(θ[1]) - model.l_calf2 * cos(θ[2])
		elseif mode == :calf
			px = model.l_thigh2 * sin(θ[1]) + model.d_calf2 * sin(θ[2])
			pz = -model.l_thigh2 * cos(θ[1]) - model.d_calf2 * cos(θ[2])
		elseif mode == :thigh
			px = model.d_thigh2 * sin(θ[1])
			pz = -model.d_thigh2 * cos(θ[1])
		else
			@error "incorrect mode specification"
		end

		p̄x = px + l_torso
		p̄y = -1.0 * sin(θ[3]) * pz + w_torso
		p̄z = cos(θ[3]) * pz

		return R * [p̄x; p̄y; p̄z] + q[1:3]
	end

	if body == :leg3
		l_torso = -model.l_torso
		w_torso = model.w_torso
		θ = q[6 + 6 .+ (1:3)]

		if mode == :ee
			px = model.l_thigh3 * sin(θ[1]) + model.l_calf3 * sin(θ[2])
			pz = -model.l_thigh3 * cos(θ[1]) - model.l_calf3 * cos(θ[2])
		elseif mode == :calf
			px = model.l_thigh3 * sin(θ[1]) + model.d_calf3 * sin(θ[2])
			pz = -model.l_thigh3 * cos(θ[1]) - model.d_calf3 * cos(θ[2])
		elseif mode == :thigh
			px = model.d_thigh3 * sin(θ[1])
			pz = -model.d_thigh3 * cos(θ[1])
		else
			@error "incorrect mode specification"
		end

		p̄x = px + l_torso
		p̄y = -1.0 * sin(θ[3]) * pz + w_torso
		p̄z = cos(θ[3]) * pz

		return R * [p̄x; p̄y; p̄z] + q[1:3]
	end

	if body == :leg4
		l_torso = -model.l_torso
		w_torso = -model.w_torso
		θ = q[6 + 9 .+ (1:3)]

		if mode == :ee
			px = model.l_thigh4 * sin(θ[1]) + model.l_calf4 * sin(θ[2])
			pz = -model.l_thigh4 * cos(θ[1]) - model.l_calf4 * cos(θ[2])
		elseif mode == :calf
			px = model.l_thigh4 * sin(θ[1]) + model.d_calf4 * sin(θ[2])
			pz = -model.l_thigh4 * cos(θ[1]) - model.d_calf4 * cos(θ[2])
		elseif mode == :thigh
			px = model.d_thigh4 * sin(θ[1])
			pz = -model.d_thigh4 * cos(θ[1])
		else
			@error "incorrect mode specification"
		end

		p̄x = px + l_torso
		p̄y = -1.0 * sin(θ[3]) * pz + w_torso
		p̄z = cos(θ[3]) * pz

		return R * [p̄x; p̄y; p̄z] + q[1:3]
	end

end

function jacobian(model::Quadruped3DAlt, q; body = :torso, mode = :ee)
	k(z) = kinematics(model, z, body = body, mode = mode)
	return ForwardDiff.jacobian(k, q)
end

# Lagrangian
function lagrangian(model::Quadruped3DAlt, q, q̇)
	L = 0.0

	## BODY
	# torso
	p_torso = kinematics(model, q, body = :torso)
	J_torso = jacobian(model, q, body = :torso)
	v_torso = J_torso * q̇

	L += 0.5 * model.m_torso * transpose(v_torso) * v_torso
	L += 0.5 * transpose(q̇[4:6]) * Diagonal(model.J_torso) * q̇[4:6]
	L -= model.m_torso * model.g * p_torso[3]

	## LEG 1
	# thigh 1
	p_thigh_1 = kinematics(model, q, body = :leg1, mode = :thigh)
	J_thigh_1 = jacobian(model, q, body = :leg1, mode = :thigh)
	v_thigh_1 = J_thigh_1 * q̇

	L += 0.5 * model.m_thigh1 * transpose(v_thigh_1) * v_thigh_1
	L -= model.m_thigh1 * model.g * p_thigh_1[3]

	# # leg 1
	# p_calf_1 = kinematics(model, q, body = :leg1, mode = :calf)
	# J_calf_1 = jacobian(model, q, body = :leg1, mode = :calf)
	# v_calf_1 = J_calf_1 * q̇
	#
	# L += 0.5 * model.m_calf1 * transpose(v_calf_1) * v_calf_1
	# L -= model.m_calf1 * model.g * p_calf_1[3]


	# # LEG 2
	# # thigh 2
	# p_thigh_2 = kinematics(model, q, body = :leg2, mode = :thigh)
	# J_thigh_2 = jacobian(model, q, body = :leg2, mode = :thigh)
	# v_thigh_2 = J_thigh_2 * q̇
	#
	# L += 0.5 * model.m_thigh2 * transpose(v_thigh_2) * v_thigh_2
	# L -= model.m_thigh2 * model.g * p_thigh_2[3]
	#
	# # leg 2
	# p_calf_2 = kinematics(model, q, body = :leg2, mode = :calf)
	# J_calf_2 = jacobian(model, q, body = :leg2, mode = :calf)
	# v_calf_2 = J_calf_2 * q̇
	#
	# L += 0.5 * model.m_calf2 * transpose(v_calf_2) * v_calf_2
	# L -= model.m_calf2 * model.g * p_calf_2[3]
	#
	#
	# # LEG 3
	# # thigh 3
	# p_thigh_3 = kinematics(model, q, body = :leg3, mode = :thigh)
	# J_thigh_3 = jacobian(model, q, body = :leg3, mode = :thigh)
	# v_thigh_3 = J_thigh_3 * q̇
	#
	# L += 0.5 * model.m_thigh3 * transpose(v_thigh_3) * v_thigh_3
	# L -= model.m_thigh3 * model.g * p_thigh_3[3]
	#
	# # leg 3
	# p_calf_3 = kinematics(model, q, body = :leg3, mode = :calf)
	# J_calf_3 = jacobian(model, q, body = :leg3, mode = :calf)
	# v_calf_3 = J_calf_3 * q̇
	#
	# L += 0.5 * model.m_calf3 * transpose(v_calf_3) * v_calf_3
	# L -= model.m_calf3 * model.g * p_calf_3[3]
	#
	# # LEG 4
	# # thigh 4
	# p_thigh_4 = kinematics(model, q, body = :leg4, mode = :thigh)
	# J_thigh_4 = jacobian(model, q, body = :leg4, mode = :thigh)
	# v_thigh_4 = J_thigh_4 * q̇
	#
	# L += 0.5 * model.m_thigh4 * transpose(v_thigh_4) * v_thigh_4
	# L -= model.m_thigh4 * model.g * p_thigh_4[3]
	#
	# # leg 4
	# p_calf_4 = kinematics(model, q, body = :leg4, mode = :calf)
	# J_calf_4 = jacobian(model, q, body = :leg4, mode = :calf)
	# v_calf_4 = J_calf_4 * q̇
	#
	# L += 0.5 * model.m_calf4 * transpose(v_calf_4) * v_calf_4
	# L -= model.m_calf4 * model.g * p_calf_4[3]

	return L
end

function _dLdq(model::Quadruped3DAlt, q, q̇)
	Lq(x) = lagrangian(model, x, q̇)
	ForwardDiff.gradient(Lq, q)
end

function _dLdq̇(model::Quadruped3DAlt, q, q̇)
	Lq̇(x) = lagrangian(model, q, x)
	ForwardDiff.gradient(Lq̇, q̇)
end

function kinematics(model::Quadruped3DAlt, q)
	p_leg_1 = kinematics(model, q, body = :leg1, mode = :ee)
	p_leg_2 = kinematics(model, q, body = :leg2, mode = :ee)
	p_leg_3 = kinematics(model, q, body = :leg3, mode = :ee)
	p_leg_4 = kinematics(model, q, body = :leg4, mode = :ee)

	SVector{12}([p_leg_1; p_leg_2; p_leg_3; p_leg_4])
end

function ϕ_func(model::Quadruped3DAlt, env::Environment, q)
	p_leg_1 = kinematics(model, q, body = :leg1, mode = :ee)
	p_leg_2 = kinematics(model, q, body = :leg2, mode = :ee)
	p_leg_3 = kinematics(model, q, body = :leg3, mode = :ee)
	p_leg_4 = kinematics(model, q, body = :leg4, mode = :ee)

	SVector{4}([p_leg_1[3] - env.surf(p_leg_1[1:2]);
				p_leg_2[3] - env.surf(p_leg_2[1:2]);
				p_leg_3[3] - env.surf(p_leg_3[1:2]);
				p_leg_4[3] - env.surf(p_leg_4[1:2])])
end

function B_func(model::Quadruped3DAlt, q)
	# @SMatrix [0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
	# 		  0.0  0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0;
	# 		  0.0  0.0 -1.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0;
	# 		  0.0  0.0  0.0  0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0;
	# 		  0.0  0.0 -1.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0;
	# 		  0.0  0.0  0.0  0.0  0.0  0.0  0.0 -1.0  1.0  0.0  0.0;
	# 		  0.0  0.0 -1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0;
	# 		  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 -1.0  1.0]

	zeros(model.dim.u, model.dim.q)
end

function A_func(model::Quadruped3DAlt, q)
	# @SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
	# 		  0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]

	zeros(model.dim.w, model.dim.q)
end

function J_func(model::Quadruped3DAlt, q)
	k(z) = kinematics(model, z)
	return ForwardDiff.jacobian(k, q)
end

function C_func(model::Quadruped3DAlt, q, q̇)
	tmp_q(z) = _dLdq̇(model, z, q̇)
	tmp_q̇(z) = _dLdq̇(model, q, z)

	ForwardDiff.jacobian(tmp_q, q) * q̇ - _dLdq(model, q, q̇)
end

function contact_forces(model::Quadruped3DAlt, env::Environment{R3, LinearizedCone}, γ1, b1, q2, k)
	m = friction_mapping(env)

	SVector{12}([transpose(rotation(env, k[1:2]))   * [m * b1[1:4];   γ1[1]];
				 transpose(rotation(env, k[4:5]))   * [m * b1[5:8];   γ1[2]];
				 transpose(rotation(env, k[7:8]))   * [m * b1[9:12];  γ1[3]];
				 transpose(rotation(env, k[10:11])) * [m * b1[13:16]; γ1[4]]])
end

function velocity_stack(model::Quadruped3DAlt, env::Environment{R3, LinearizedCone}, q1, q2, k, h)
	v = J_func(model, q2) * (q2 - q1) / h[1]

	v1_surf = rotation(env, k[1:2]) * v[1:3]
	v2_surf = rotation(env, k[4:5]) * v[4:6]
	v3_surf = rotation(env, k[7:8]) * v[7:9]
	v4_surf = rotation(env, k[10:11]) * v[10:12]

	SVector{16}([transpose(friction_mapping(env)) * v1_surf[1:2];
				 transpose(friction_mapping(env)) * v2_surf[1:2];
				 transpose(friction_mapping(env)) * v3_surf[1:2];
				 transpose(friction_mapping(env)) * v4_surf[1:2]])
end


################################################################################
# Instantiation
################################################################################
# Dimensions
nq = 6 + 4 * 3           # configuration dimension
nu = 4 * 3                # control dimension
nc = 4                    # number of contact points
nf = 4                    # number of parameters for friction cone
nb = nc * nf
nw = 3

# World parameters
g = 9.81      # gravity
μ_world = 1.0 # coefficient of friction
μ_joint = 0.1 # coefficient of torque friction at the joints

# ~Unitree A1
# Model parameters
m_torso = 4.713 + 4 * 0.696
m_thigh = 1.013
m_leg = 0.166

J_torso = 0.01683 + 4 * 0.696 * 0.183^2.0
J_thigh = 0.00552
J_leg = 0.00299

l_torso = 0.183*2
l_thigh = 0.2
l_leg = 0.2

d_torso = 0.5 * l_torso + 0.0127
d_thigh = 0.5 * l_thigh - 0.00323
d_leg = 0.5 * l_leg - 0.006435

w_torso = 0.5 * 0.203

quadruped = Quadruped3DAlt(Dimensions(nq, nu, nw, nc),
				g, μ_world, μ_joint,
				l_torso, d_torso, w_torso, m_torso, [J_torso, J_torso, J_torso],
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				BaseMethods(), DynamicsMethods(),
				SVector{nq}([zeros(3); μ_joint * ones(nq - 3)]))


lagrangian(quadruped, rand(quadruped.dim.q), rand(quadruped.dim.q))
model = quadruped



nv = model.dim.q

nq = model.dim.q
nu = model.dim.u
nw = model.dim.w
nc = model.dim.c

# Declare variables
@variables q[1:nq]
@variables q̇[1:nv]

L(z) = lagrangian(model, z[1:nq], z[nq .+ (1:nq)])
ForwardDiff.gradient(L, [q; q̇])

# Lagrangian
L = lagrangian(model, q, q̇)
# L = Symbolics.simplify.(L, expand = true)

dLq = Symbolics.gradient(L, q) # including mapping for orientation (e.g., attitude Jacobian)
# dLq = simplify.(dLq, expand = true)
dLq̇ = Symbolics.gradient(L, q̇)
# dLq̇ = simplify.(dLq, expand = true)
ddL = Symbolics.hessian(L, [q; q̇])
# ddL = simplify.(ddL, expand = true)
# ddLq̇q = ddL[nq .+ (1:nv), 1:nq]

M = ddL[nq .+ (1:nv), nq .+ (1:nv)]

# Coriolis and Centrifugal forces Jacobians
C = ddLq̇q * q̇ - dLq
C = Symbolics.simplify.(C, expand = true)

# Control input Jacobian
B = B_func(model, q)
B = reshape(B, (nu, nv))
B = Symbolics.simplify.(B, expand = true)

# Disturbance input Jacobian
A = A_func(model, q)
A = reshape(A, (nw, nv))
A = Symbolics.simplify.(A, expand = true)

# Contact Jacobian
J = J_func(model, q)
J = reshape(J, size(J_func(model, zeros(nq))))
J = Symbolics.simplify.(J, expand = true)

# Kinematics
k = kinematics(model, q)
k = simplify.(k, expand = true)

# Build function
expr = Dict{Symbol, Expr}()
expr[:L]    = build_function([L], q, q̇)[1] # need to transpose to get a line vector
expr[:M]    = build_function(M, q)[1]
expr[:B]    = build_function(B, q)[1]
expr[:A]    = build_function(A, q)[1]
expr[:J]    = build_function(J, q)[1]
expr[:C]    = build_function(C, q, q̇)[1]
expr[:k]    = build_function(k, q)[1]

################################################################################
dir = joinpath(pwd(), "src/dynamics/quadruped_3D_alt")
model = deepcopy(quadruped)
path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

# expr_base = generate_base_expressions(model, M_analytical = false)
save_expressions(expr, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

# model.base.C(rand(model.dim.q), rand(model.dim.q))
# M_fast(model, rand(model.dim.q))
# J_fast(model, rand(model.dim.q))
# C_fast(model, rand(model.dim.q), rand(model.dim.q))
# nq = model.dim.q
# nu = model.dim.u
# nw = model.dim.w
# nc = model.dim.c
# ncf = size(J_func(model, zeros(nq)))[1]
#
# # Declare variables
# @variables q0[1:nq]
# @variables q1[1:nq]
# @variables u1[1:nu]
# @variables w1[1:nw]
# @variables λ1[1:ncf]
# @variables q2[1:nq]
# @variables h[1:1]
#
# # Expressions
# expr = Dict{Symbol, Expr}()
#
# # Dynamics
# d = dynamics(model, h, q0, q1, u1, w1, λ1, q2)
# d = Symbolics.simplify.(d)
