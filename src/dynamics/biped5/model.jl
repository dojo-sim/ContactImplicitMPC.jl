struct Biped5{T} <: ContactDynamicsModel
	dim::Dimensions

    g::T
    μ_world::T
	μ_joint::T

    # torso
    l_torso
    d_torso
    m_torso
    J_torso

    # leg 1
        # thigh
    l_thigh1
    d_thigh1
    m_thigh1
    J_thigh1

        # calf
    l_calf1
    d_calf1
    m_calf1
    J_calf1

    # leg 2
        # thigh
    l_thigh2
    d_thigh2
    m_thigh2
    J_thigh2

        # calf
    l_calf2
    d_calf2
    m_calf2
    J_calf2

    # joint limits
    qL
    qU

    # torque limits
    uL
    uU

    nq
    nu
    nc
    nf
    nb
    ns

    idx_u
    idx_λ
    idx_b
    idx_ψ
    idx_η
    idx_s

	alt

	# fast methods
	base
	dyn
	res
	linearized

	spa::SparseStructure

	joint_friction

	env::Environment
end

function kinematics_1(model::Biped5, q; body = :torso, mode = :ee)
	x = q[1]
	z = q[2]

	if body == :torso
		l = model.l_torso
		d = model.d_torso
		θ = q[3]
		if mode == :ee
			return [x - l * sin(θ); z + l * cos(θ)]
		elseif mode == :com
			return [x - d * sin(θ); z + d * cos(θ)]
		end
	elseif body == :thigh_1
		l = model.l_thigh1
		d = model.d_thigh1
		θ = q[4]
	elseif body == :thigh_2
		l = model.l_thigh2
		d = model.d_thigh2
		θ = q[6]
	else
		@error "incorrect body specification"
	end

	if mode == :ee
		return [x + l * sin(θ); z - l * cos(θ)]
	elseif mode == :com
		return [x + d * sin(θ); z - d * cos(θ)]
	else
		@error "incorrect mode specification"
	end
end

function jacobian_1(model::Biped5, q; body = :torso, mode = :ee)
	jac = zeros(eltype(q), 2, model.dim.q)
	jac[1, 1] = 1.0
	jac[2, 2] = 1.0
	if body == :torso
		r = mode == :ee ? model.l_torso : model.d_torso
		θ = q[3]
		jac[1, 3] = -r * cos(θ)
		jac[2, 3] = -r * sin(θ)
	elseif body == :thigh_1
		r = mode == :ee ? model.l_thigh1 : model.d_thigh1
		θ = q[4]
		jac[1, 4] = r * cos(θ)
		jac[2, 4] = r * sin(θ)
	elseif body == :thigh_2
		r = mode == :ee ? model.l_thigh2 : model.d_thigh2
		θ = q[6]
		jac[1, 6] = r * cos(θ)
		jac[2, 6] = r * sin(θ)
	else
		@error "incorrect body specification"
	end

	return jac
end

function kinematics_2(model::Biped5, q; body = :calf_1, mode = :ee)

	if body == :calf_1
		p = kinematics_1(model, q, body = :thigh_1, mode = :ee)

		θb = q[5]

		lb = model.l_calf1
		db = model.d_calf1
	elseif body == :calf_2
		p = kinematics_1(model, q, body = :thigh_2, mode = :ee)

		θb = q[7]

		lb = model.l_calf2
		db = model.d_calf2
	else
		@error "incorrect body specification"
	end

	if mode == :ee
		return p + [lb * sin(θb); -1.0 * lb * cos(θb)]
	elseif mode == :com
		return p + [db * sin(θb); -1.0 * db * cos(θb)]
	else
		@error "incorrect mode specification"
	end
end

function jacobian_2(model::Biped5, q; body = :calf_1, mode = :ee)

	if body == :calf_1
		jac = jacobian_1(model, q, body = :thigh_1, mode = :ee)

		θb = q[5]

		r = mode == :ee ? model.l_calf1 : model.d_calf1

		jac[1, 5] += r * cos(θb)
		jac[2, 5] += r * sin(θb)
	elseif body == :calf_2
		jac = jacobian_1(model, q, body = :thigh_2, mode = :ee)

		θb = q[7]

		r = mode == :ee ? model.l_calf2 : model.d_calf2

		jac[1, 7] += r * cos(θb)
		jac[2, 7] += r * sin(θb)
	else
		@error "incorrect body specification"
	end

	return jac
end

# Lagrangian

function lagrangian(model::Biped5, q, q̇)
	L = 0.0

	# torso
	p_torso = kinematics_1(model, q, body = :torso, mode = :com)
	J_torso = jacobian_1(model, q, body = :torso, mode = :com)
	v_torso = J_torso * q̇

	L += 0.5 * model.m_torso * transpose(v_torso) * v_torso
	L += 0.5 * model.J_torso * q̇[3]^2.0
	L -= model.m_torso * model.g * p_torso[2]

	# thigh 1
	p_thigh_1 = kinematics_1(model, q, body = :thigh_1, mode = :com)
	J_thigh_1 = jacobian_1(model, q, body = :thigh_1, mode = :com)
	v_thigh_1 = J_thigh_1 * q̇

	L += 0.5 * model.m_thigh1 * transpose(v_thigh_1) * v_thigh_1
	L += 0.5 * model.J_thigh1 * q̇[4]^2.0
	L -= model.m_thigh1 * model.g * p_thigh_1[2]

	# leg 1
	p_calf_1 = kinematics_2(model, q, body = :calf_1, mode = :com)
	J_calf_1 = jacobian_2(model, q, body = :calf_1, mode = :com)
	v_calf_1 = J_calf_1 * q̇

	L += 0.5 * model.m_calf1 * transpose(v_calf_1) * v_calf_1
	L += 0.5 * model.J_calf1 * q̇[5]^2.0
	L -= model.m_calf1 * model.g * p_calf_1[2]

	# thigh 2
	p_thigh_2 = kinematics_1(model, q, body = :thigh_2, mode = :com)
	J_thigh_2 = jacobian_1(model, q, body = :thigh_2, mode = :com)
	v_thigh_2 = J_thigh_2 * q̇

	L += 0.5 * model.m_thigh2 * transpose(v_thigh_2) * v_thigh_2
	L += 0.5 * model.J_thigh2 * q̇[6]^2.0
	L -= model.m_thigh2 * model.g * p_thigh_2[2]

	# leg 2
	p_calf_2 = kinematics_2(model, q, body = :calf_2, mode = :com)
	J_calf_2 = jacobian_2(model, q, body = :calf_2, mode = :com)
	v_calf_2 = J_calf_2 * q̇

	L += 0.5 * model.m_calf2 * transpose(v_calf_2) * v_calf_2
	L += 0.5 * model.J_calf2 * q̇[7]^2.0
	L -= model.m_calf2 * model.g * p_calf_2[2]

	return L
end

function _dLdq(model::Biped5, q, q̇)
	Lq(x) = lagrangian(model, x, q̇)
	ForwardDiff.gradient(Lq, q)
end

function _dLdq̇(model::Biped5, q, q̇)
	Lq̇(x) = lagrangian(model, q, x)
	ForwardDiff.gradient(Lq̇, q̇)
end

function _C_func(model::Biped5, q, q̇)
	tmp_q(z) = _dLdq̇(model, z, q̇)
	tmp_q̇(z) = _dLdq̇(model, q, z)

	ForwardDiff.jacobian(tmp_q, q) * q̇ - _dLdq(model, q, q̇)
end

function kinematics(model::Biped5, q)
	p_calf_1 = kinematics_2(model, q, body = :calf_1, mode = :ee)
	p_calf_2 = kinematics_2(model, q, body = :calf_2, mode = :ee)

	SVector{4}([p_calf_1; p_calf_2])
end

# Methods
function M_func(model::Biped5, q)
	M = Diagonal([0.0, 0.0, model.J_torso, model.J_thigh1, model.J_calf1, model.J_thigh2, model.J_calf2])

	# torso
	J_torso = jacobian_1(model, q, body = :torso, mode = :com)
	M += model.m_torso * transpose(J_torso) * J_torso

	# thigh 1
	J_thigh_1 = jacobian_1(model, q, body = :thigh_1, mode = :com)
	M += model.m_thigh1 * transpose(J_thigh_1) * J_thigh_1

	# leg 1
	J_calf_1 = jacobian_2(model, q, body = :calf_1, mode = :com)
	M += model.m_calf1 * transpose(J_calf_1) * J_calf_1

	# thigh 2
	J_thigh_2 = jacobian_1(model, q, body = :thigh_2, mode = :com)
	M += model.m_thigh2 * transpose(J_thigh_2) * J_thigh_2

	# leg 2
	J_calf_2 = jacobian_2(model, q, body = :calf_2, mode = :com)
	M += model.m_calf2 * transpose(J_calf_2) * J_calf_2

	return M
end

function _C_func(model::Biped5, q, q̇)
	tmp_q(z) = _dLdq̇(model, z, q̇)
	tmp_q̇(z) = _dLdq̇(model, q, z)

	ForwardDiff.jacobian(tmp_q, q) * q̇ - _dLdq(model, q, q̇)
end

function ϕ_func(model::Biped5, q)
	p_calf_1 = kinematics_2(model, q, body = :calf_1, mode = :ee)
	p_calf_2 = kinematics_2(model, q, body = :calf_2, mode = :ee)

	@SVector [p_calf_1[2], p_calf_2[2]]
end

function B_func(model::Biped5, q)
	@SMatrix [0.0 0.0 1.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 0.0 1.0 0.0 0.0 0.0;
			  0.0 0.0 0.0 0.0 1.0 0.0 0.0;
			  0.0 0.0 0.0 0.0 0.0 1.0 0.0;
			  0.0 0.0 0.0 0.0 0.0 0.0 1.0]
end

function A_func(model::Biped5, q)
	@SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 1.0 0.0 0.0 0.0 0.0 0.0]
end

function J_func(model::Biped5, q)
	J_calf_1 = jacobian_2(model, q, body = :calf_1, mode = :ee)
	J_calf_2 = jacobian_2(model, q, body = :calf_2, mode = :ee)

	return [J_calf_1; J_calf_2]
end

# Dimensions
nq = 2 + 5                # configuration dimension
nu = 5                    # control dimension
nc = 2                    # number of contact points
nf = 2                    # number of parameters for friction cone
nb = nc * nf
nw = 2

# World parameters
μ_world = 1.0      # coefficient of friction
μ_joint = 0.1
g = 9.81     # gravity

# Model parameters
m_torso = 0.5 + 0.48 * 2.0
m_thigh = 0.8112
m_calf = 0.3037

J_torso = 0.0029
J_thigh = 0.00709
J_calf = 0.00398

l_torso = 0.15 + 0.15
l_thigh = 0.2755
l_calf = 0.308

d_torso = 0.0342
d_thigh = 0.2176
d_calf = 0.1445

model = Biped5{Dimensions(nq, nu, nw, nc, nb),
			  g, μ_world, μ_joint,
			  l_torso, d_torso, m_torso, J_torso,
			  l_thigh, d_thigh, m_thigh, J_thigh,
			  l_calf, d_calf, m_calf, J_calf,
			  l_foot, d_foot, m_foot, J_foot,
			  l_thigh, d_thigh, m_thigh, J_thigh,
			  l_calf, d_calf, m_calf, J_calf,
			  l_foot, d_foot, m_foot, J_foot,
			  zeros(nc),
			  BaseMethods(), DynamicsMethods(), ResidualMethods(), ResidualMethods(),
			  SparseStructure(spzeros(0, 0), spzeros(0, 0)),
			  SVector{nq}([zeros(3); 0.0 * μ_joint * ones(nq - 3)]),
			  environment_2D_flat())

function initial_configuration(model, θ_torso, θ_thigh_1, θ_leg_1, θ_thigh_2)
    q1 = zeros(model.dim.q)
    q1[3] = θ_torso
    q1[4] = θ_thigh_1
    q1[5] = θ_leg_1
    z1 = model.l_thigh1 * cos(q1[4]) + model.l_calf1 * cos(q1[5])
    q1[6] = θ_thigh_2
    q1[7] = -1.0 * acos((z1 - model.l_thigh2 * cos(q1[6])) / model.l_calf2)
    q1[2] = z1

    p1 = kinematics_2(model, q1, body = :calf_1, mode = :ee)
    p2 = kinematics_2(model, q1, body = :calf_2, mode = :ee)
    @show stride = abs(p1[1] - p2[1])

    q1[1] = -1.0 * p1[1]

    qM = copy(q1)
    qM[4] = q1[6]
    qM[5] = q1[7]
    qM[6] = q1[4]
    qM[7] = q1[5]
    qM[1] = abs(p2[1])

    pM_1 = kinematics_2(model, qM, body = :calf_1, mode = :ee)
    pM_2 = kinematics_2(model, qM, body = :calf_2, mode = :ee)

    qT = copy(q1)
    qT[1] = 2 * stride

    pT_1 = kinematics_2(model, qT, body = :calf_1, mode = :ee)
    pT_2 = kinematics_2(model, qT, body = :calf_2, mode = :ee)

    return q1, qM, qT
end
