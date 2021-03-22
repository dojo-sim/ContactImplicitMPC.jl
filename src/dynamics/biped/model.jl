struct Biped{T} <: ContactDynamicsModel
    dim::Dimensions

    g::T
    μ_world::T
	μ_joint::T

    # torso
    l_torso::T
    d_torso::T
    m_torso::T
    J_torso::T

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

		# foot
	l_foot1::T # toe length
	d_foot1::T # heel length
	m_foot1::T
	J_foot1::T

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

		# foot
	l_foot2::T # toe length
	d_foot2::T # heel length
	m_foot2::T
	J_foot2::T

	alt

	# fast methods
	base
	dyn
	res
	approx

	spa::SparseStructure

	joint_friction

	env::Environment
end

function kinematics_1(model::Biped, q; body = :torso, mode = :ee)
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

function jacobian_1(model::Biped, q; body = :torso, mode = :ee)
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

function kinematics_2(model::Biped, q; body = :calf_1, mode = :ee)

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

function jacobian_2(model::Biped, q; body = :calf_1, mode = :ee)

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

function kinematics_3(model::Biped, q; body = :foot_1, mode = :ee)

	if body == :foot_1
		p = kinematics_2(model, q, body = :calf_1, mode = :ee)

		θb = q[8]

		lb = model.l_foot1
		db = model.d_foot1
		cb = 0.5 * (model.l_foot1 - model.d_foot1)
	elseif body == :foot_2
		p = kinematics_2(model, q, body = :calf_2, mode = :ee)

		θb = q[9]

		lb = model.l_foot2
		db = model.d_foot2
		cb = 0.5 * (model.l_foot2 - model.d_foot2)
	else
		@error "incorrect body specification"
	end

	if mode == :toe
		return p + [lb * sin(θb); -1.0 * lb * cos(θb)]
	elseif mode == :heel
		return p + [-db * sin(θb); 1.0 * db * cos(θb)]
	elseif mode == :com
		return p + [cb * sin(θb); -1.0 * cb * cos(θb)]
	else
		@error "incorrect mode specification"
	end
end

function jacobian_3(model::Biped, q; body = :foot_1, mode = :ee)

	if body == :foot_1
		jac = jacobian_2(model, q, body = :calf_1, mode = :ee)

		θb = q[8]

		if mode == :toe
			r = model.l_foot1
		elseif mode == :heel
			r = -1.0 * model.d_foot1
		elseif mode == :com
			r = 0.5 * (model.l_foot1 - model.d_foot1)
		else
			@error "incorrect mode specification"
		end

		jac[1, 8] += r * cos(θb)
		jac[2, 8] += r * sin(θb)

	elseif body == :foot_2
		jac = jacobian_2(model, q, body = :calf_2, mode = :ee)

		θb = q[9]

		if mode == :toe
			r = model.l_foot2
		elseif mode == :heel
			r = -1.0 * model.d_foot2
		elseif mode == :com
			r = 0.5 * (model.l_foot2 - model.d_foot2)
		else
			@error "incorrect mode specification"
		end
		jac[1, 9] += r * cos(θb)
		jac[2, 9] += r * sin(θb)

	else
		@error "incorrect body specification"
	end

	return jac
end

# Lagrangian

function lagrangian(model::Biped, q, q̇)
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

	# foot 1
	p_foot_1 = kinematics_3(model, q, body = :foot_1, mode = :com)
	J_foot_1 = jacobian_3(model, q, body = :foot_1, mode = :com)
	v_foot_1 = J_foot_1 * q̇

	L += 0.5 * model.m_foot1 * transpose(v_foot_1) * v_foot_1
	L += 0.5 * model.J_foot1 * q̇[8]^2.0
	L -= model.m_foot1 * model.g * p_foot_1[2]

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

	# foot 2
	p_foot_2 = kinematics_3(model, q, body = :foot_2, mode = :com)
	J_foot_2 = jacobian_3(model, q, body = :foot_2, mode = :com)
	v_foot_2 = J_foot_2 * q̇

	L += 0.5 * model.m_foot2 * transpose(v_foot_2) * v_foot_2
	L += 0.5 * model.J_foot2 * q̇[9]^2.0
	L -= model.m_foot2 * model.g * p_foot_2[2]

	return L
end

function _dLdq(model::Biped, q, q̇)
	Lq(x) = lagrangian(model, x, q̇)
	ForwardDiff.gradient(Lq, q)
end

function _dLdq̇(model::Biped, q, q̇)
	Lq̇(x) = lagrangian(model, q, x)
	ForwardDiff.gradient(Lq̇, q̇)
end

function _C_func(model::Biped, q, q̇)
	tmp_q(z) = _dLdq̇(model, z, q̇)
	tmp_q̇(z) = _dLdq̇(model, q, z)

	ForwardDiff.jacobian(tmp_q, q) * q̇ - _dLdq(model, q, q̇)
end

# Methods
function M_func(model::Biped, q)
	M = Diagonal([0.0, 0.0,
		model.J_torso,
		model.J_thigh1, model.J_calf1,
		model.J_thigh2, model.J_calf2,
		model.J_foot1, model.J_foot2])

	# torso
	J_torso = jacobian_1(model, q, body = :torso, mode = :com)
	M += model.m_torso * transpose(J_torso) * J_torso

	# thigh 1
	J_thigh_1 = jacobian_1(model, q, body = :thigh_1, mode = :com)
	M += model.m_thigh1 * transpose(J_thigh_1) * J_thigh_1

	# leg 1
	J_calf_1 = jacobian_2(model, q, body = :calf_1, mode = :com)
	M += model.m_calf1 * transpose(J_calf_1) * J_calf_1

	# foot 1
	J_foot_1 = jacobian_3(model, q, body = :foot_1, mode = :com)
	M += model.m_foot1 * transpose(J_foot_1) * J_foot_1

	# thigh 2
	J_thigh_2 = jacobian_1(model, q, body = :thigh_2, mode = :com)
	M += model.m_thigh2 * transpose(J_thigh_2) * J_thigh_2

	# leg 2
	J_calf_2 = jacobian_2(model, q, body = :calf_2, mode = :com)
	M += model.m_calf2 * transpose(J_calf_2) * J_calf_2

	# foot 2
	J_foot_2 = jacobian_3(model, q, body = :foot_2, mode = :com)
	M += model.m_foot2 * transpose(J_foot_2) * J_foot_2

	return M
end

function ϕ_func(model::Biped, q)
	p_toe_1 = kinematics_3(model, q, body = :foot_1, mode = :toe)
	p_heel_1 = kinematics_3(model, q, body = :foot_1, mode = :heel)
	p_toe_2 = kinematics_3(model, q, body = :foot_2, mode = :toe)
	p_heel_2 = kinematics_3(model, q, body = :foot_2, mode = :heel)

	@SVector [p_toe_1[2], p_heel_1[2], p_toe_2[2], p_heel_2[2]]
end

function B_func(model::Biped, q)
	@SMatrix [0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0;
			  0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0;
			  0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0;
			  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]
end

function A_func(model::Biped, q)
	@SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
end

function J_func(model::Biped, q)
	J_toe_1 = jacobian_3(model, q, body = :foot_1, mode = :toe)
	J_heel_1 = jacobian_3(model, q, body = :foot_1, mode = :heel)
	J_toe_2 = jacobian_3(model, q, body = :foot_2, mode = :toe)
	J_heel_2 = jacobian_3(model, q, body = :foot_2, mode = :heel)

	return [J_toe_1;
			J_heel_1;
			J_toe_2;
			J_heel_2]
end

# Dimensions
nq = 2 + 5 + 2            # configuration dimension
nu = 7                    # control dimension
nc = 4                    # number of contact points
nf = 2                    # number of parameters for friction cone
nb = nc * nf              # number of friction parameters
nw = 2                    # disturbance dimension

# World parameters
μ_world = 0.5      # coefficient of friction
μ_joint = 0.1
g = 9.81     # gravity

# Model parameters
m_torso = 0.5 + 0.48 * 2.0
m_thigh = 0.8112
m_calf = 0.3037
m_foot = 0.4 * m_calf

l_torso = 0.15 + 0.15
l_thigh = 0.2755
l_calf = 0.308
l_foot = 0.4 * l_calf

d_torso = 0.0342
d_thigh = 0.2176
d_calf = 0.1445
d_foot = 0.4 * d_calf

J_torso = 0.0029
J_thigh = 0.00709
J_calf = 0.00398
J_foot = 0.4 * J_calf # 1.0 / 12.0 * m_foot * (l_foot + d_foot)^2.0

biped = Biped(Dimensions(nq, nu, nw, nc, nb),
			  g, μ_world, μ_joint,
			  l_torso, d_torso, m_torso, J_torso,
			  l_thigh, d_thigh, m_thigh, J_thigh,
			  l_calf, d_calf, m_calf, J_calf,
			  l_foot, d_foot, m_foot, J_foot,
			  l_thigh, d_thigh, m_thigh, J_thigh,
			  l_calf, d_calf, m_calf, J_calf,
			  l_foot, d_foot, m_foot, J_foot,
			  zeros(nc),
			  BaseMethods(), DynamicsMethods(), ResidualMethods(), ApproximateMethods(),
			  SparseStructure(spzeros(0, 0), spzeros(0, 0)),
			  SVector{nq}([zeros(3); μ_joint * ones(nq - 3)]),
			  environment_2D_flat())
