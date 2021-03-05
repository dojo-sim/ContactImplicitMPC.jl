mutable struct QuadrupedBasic20{T} <: QuadrupedModel
	dim::Dimensions12
	ind::Indices13
	bas::BaseMethods12
	dyn::DynamicsMethods13
	res::ResidualMethods14
	alt::AbstractVector{T}

	h::T
	g::T
	μ::T
	joint_friction::T

    # n::Int
    # m::Int
    # d::Int
	#
	# # dim
	# alt1
	# alt2
	# alt3
	# alt4
	#
	# h
    # g
    # μ
	# joint_friction

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

	# leg 3
        # thigh
    l_thigh3
    d_thigh3
    m_thigh3
    J_thigh3
        # calf
    l_calf3
    d_calf3
    m_calf3
    J_calf3

	# leg 4
        # thigh
    l_thigh4
    d_thigh4
    m_thigh4
    J_thigh4
        # calf
    l_calf4
    d_calf4
    m_calf4
    J_calf4

    # # joint limits
    # qL
    # qU
	#
    # # torque limits
    # uL
    # uU
	#
    # nq
    # nu
    # nc
    # nf
    # nb
    # ns
	#
    # idx_u
    # idx_λ
    # idx_b
    # idx_ψ
    # idx_η
    # idx_s
end

function kinematics_1(model::QuadrupedBasic20, q; body = :torso, mode = :ee)
	x = q[1]
	z = q[2]

	if body == :torso
		l = model.l_torso
		d = model.d_torso
		θ = q[3]
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

function jacobian_1(model::QuadrupedBasic20, q; body = :torso, mode = :ee)
	jac = zeros(eltype(q), 2, model.dim.q)
	jac[1, 1] = 1.0
	jac[2, 2] = 1.0
	if body == :torso
		r = mode == :ee ? model.l_torso : model.d_torso
		θ = q[3]
		jac[1, 3] = r * cos(θ)
		jac[2, 3] = r * sin(θ)
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

function kinematics_2(model::QuadrupedBasic20, q; body = :calf_1, mode = :ee)

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
	elseif body == :thigh_3
		p = kinematics_1(model, q, body = :torso, mode = :ee)
		θb = q[8]

		lb = model.l_thigh3
		db = model.d_thigh3
	elseif body == :thigh_4
		p = kinematics_1(model, q, body = :torso, mode = :ee)
		θb = q[10]

		lb = model.l_thigh4
		db = model.d_thigh4
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

function jacobian_2(model::QuadrupedBasic20, q; body = :calf_1, mode = :ee)

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
	elseif body == :thigh_3
		jac = jacobian_1(model, q, body = :torso, mode = :ee)
		θb = q[8]

		r = mode == :ee ? model.l_thigh3 : model.d_thigh3

		jac[1, 8] += r * cos(θb)
		jac[2, 8] += r * sin(θb)
	elseif body == :thigh_4
		jac = jacobian_1(model, q, body = :torso, mode = :ee)

		θb = q[10]

		r = mode == :ee ? model.l_thigh4 : model.d_thigh4

		jac[1, 10] += r * cos(θb)
		jac[2, 10] += r * sin(θb)
	else
		@error "incorrect body specification"
	end

	return jac
end

function kinematics_3(model::QuadrupedBasic20, q; body = :calf_3, mode = :ee)

	if body == :calf_3
		p = kinematics_2(model, q, body = :thigh_3, mode = :ee)
		θc = q[9]


		lb = model.l_calf3
		db = model.d_calf3
	elseif body == :calf_4
		p = kinematics_2(model, q, body = :thigh_4, mode = :ee)

		θc = q[11]

		lb = model.l_calf4
		db = model.d_calf4
	else
		@error "incorrect body specification"
	end

	if mode == :ee
		return p + [lb * sin(θc); -1.0 * lb * cos(θc)]
	elseif mode == :com
		return p + [db * sin(θc); -1.0 * db * cos(θc)]
	else
		@error "incorrect mode specification"
	end
end

function jacobian_3(model::QuadrupedBasic20, q; body = :calf_3, mode = :ee)

	if body == :calf_3
		jac = jacobian_2(model, q, body = :thigh_3, mode = :ee)

		θc = q[9]

		r = mode == :ee ? model.l_calf3 : model.d_calf3

		jac[1, 9] += r * cos(θc)
		jac[2, 9] += r * sin(θc)

	elseif body == :calf_4
		jac = jacobian_2(model, q, body = :thigh_4, mode = :ee)

		θc = q[11]

		r = mode == :ee ? model.l_calf4 : model.d_calf4

		jac[1, 11] += r * cos(θc)
		jac[2, 11] += r * sin(θc)
	else
		@error "incorrect body specification"
	end

	return jac
end

# Lagrangian

function lagrangian(model::QuadrupedBasic20, q, q̇)
	L = 0.0

	# # torso
	p_torso = kinematics_1(model, q, body = :torso, mode = :com)
	J_torso = jacobian_1(model, q, body = :torso, mode = :com)
	v_torso = J_torso * q̇
	#
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

	# thigh 3
	p_thigh_3 = kinematics_2(model, q, body = :thigh_3, mode = :com)
	J_thigh_3 = jacobian_2(model, q, body = :thigh_3, mode = :com)
	v_thigh_3 = J_thigh_3 * q̇

	L += 0.5 * model.m_thigh3 * transpose(v_thigh_3) * v_thigh_3
	L += 0.5 * model.J_thigh3 * q̇[8]^2.0
	L -= model.m_thigh3 * model.g * p_thigh_3[2]

	# leg 3
	p_calf_3 = kinematics_3(model, q, body = :calf_3, mode = :com)
	J_calf_3 = jacobian_3(model, q, body = :calf_3, mode = :com)
	v_calf_3 = J_calf_3 * q̇

	L += 0.5 * model.m_calf3 * transpose(v_calf_3) * v_calf_3
	L += 0.5 * model.J_calf3 * q̇[9]^2.0
	L -= model.m_calf3 * model.g * p_calf_3[2]

	# thigh 4
	p_thigh_4 = kinematics_2(model, q, body = :thigh_4, mode = :com)
	J_thigh_4 = jacobian_2(model, q, body = :thigh_4, mode = :com)
	v_thigh_4 = J_thigh_4 * q̇

	L += 0.5 * model.m_thigh4 * transpose(v_thigh_4) * v_thigh_4
	L += 0.5 * model.J_thigh4 * q̇[10]^2.0
	L -= model.m_thigh4 * model.g * p_thigh_4[2]

	# leg 4
	p_calf_4 = kinematics_3(model, q, body = :calf_4, mode = :com)
	J_calf_4 = jacobian_3(model, q, body = :calf_4, mode = :com)
	v_calf_4 = J_calf_4 * q̇

	L += 0.5 * model.m_calf4 * transpose(v_calf_4) * v_calf_4
	L += 0.5 * model.J_calf4 * q̇[11]^2.0
	L -= model.m_calf4 * model.g * p_calf_4[2]

	return L
end

function _dLdq(model::QuadrupedBasic20, q, q̇)
	Lq(x) = lagrangian(model, x, q̇)
	ForwardDiff.gradient(Lq, q)
end

function _dLdq̇(model::QuadrupedBasic20, q, q̇)
	Lq̇(x) = lagrangian(model, q, x)
	ForwardDiff.gradient(Lq̇, q̇)
end


# Methods
function M_func(model::QuadrupedBasic20, q)
	M = Diagonal([0.0, 0.0, model.J_torso, model.J_thigh1, model.J_calf1, model.J_thigh2, model.J_calf2, model.J_thigh3, model.J_calf3, model.J_thigh4, model.J_calf4])

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

	# thigh 3
	J_thigh_3 = jacobian_2(model, q, body = :thigh_3, mode = :com)
	M += model.m_thigh3 * transpose(J_thigh_3) * J_thigh_3

	# leg 3
	J_calf_3 = jacobian_3(model, q, body = :calf_3, mode = :com)
	M += model.m_calf3 * transpose(J_calf_3) * J_calf_3

	# thigh 4
	J_thigh_4 = jacobian_2(model, q, body = :thigh_4, mode = :com)
	M += model.m_thigh4 * transpose(J_thigh_4) * J_thigh_4

	# leg 4
	J_calf_4 = jacobian_3(model, q, body = :calf_4, mode = :com)
	M += model.m_calf4 * transpose(J_calf_4) * J_calf_4

	return M
end

function _C_func(model::QuadrupedBasic20, q, q̇)
	tmp_q(z) = _dLdq̇(model, z, q̇)
	tmp_q̇(z) = _dLdq̇(model, q, z)

	ForwardDiff.jacobian(tmp_q, q) * q̇ - _dLdq(model, q, q̇)
end


function ϕ_func(model::QuadrupedBasic20, q)
	p_calf_1 = kinematics_2(model, q, body = :calf_1, mode = :ee)
	p_calf_2 = kinematics_2(model, q, body = :calf_2, mode = :ee)
	p_calf_3 = kinematics_3(model, q, body = :calf_3, mode = :ee)
	p_calf_4 = kinematics_3(model, q, body = :calf_4, mode = :ee)
	alt = model.alt

	@SVector [p_calf_1[2]-alt[1], p_calf_2[2]-alt[2], p_calf_3[2]-alt[3], p_calf_4[2]-alt[4]]
end

# function ϕ(model::QuadrupedBasic20, q)
# 	return ϕ_func(model, q)
# end
#
# function ϕ_no_alt(model::QuadrupedBasic20, q)
# 	p_calf_1 = kinematics_2(model, q, body = :calf_1, mode = :ee)
# 	p_calf_2 = kinematics_2(model, q, body = :calf_2, mode = :ee)
# 	p_calf_3 = kinematics_3(model, q, body = :calf_3, mode = :ee)
# 	p_calf_4 = kinematics_3(model, q, body = :calf_4, mode = :ee)
#
# 	@SVector [p_calf_1[2], p_calf_2[2], p_calf_3[2], p_calf_4[2]]
# end

# function tangential_contact_pos(model::QuadrupedBasic20, q)
# 	p_calf_1 = kinematics_2(model, q, body = :calf_1, mode = :ee)
# 	p_calf_2 = kinematics_2(model, q, body = :calf_2, mode = :ee)
# 	p_calf_3 = kinematics_3(model, q, body = :calf_3, mode = :ee)
# 	p_calf_4 = kinematics_3(model, q, body = :calf_4, mode = :ee)
#
# 	@SVector [p_calf_1[1], 0., p_calf_2[1], 0., p_calf_3[1], 0., p_calf_4[1], 0.]
# end

function B_func(model::QuadrupedBasic20, q)
	@SMatrix [0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0;
			  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0;
			  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0;
			  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 q[1]]
end

function N_func(model::QuadrupedBasic20, q)
	J_calf_1 = jacobian_2(model, q, body = :calf_1, mode = :ee)
	J_calf_2 = jacobian_2(model, q, body = :calf_2, mode = :ee)
	J_calf_3 = jacobian_3(model, q, body = :calf_3, mode = :ee)
	J_calf_4 = jacobian_3(model, q, body = :calf_4, mode = :ee)

	return [view(J_calf_1, 2:2, :);
			view(J_calf_2, 2:2, :);
			view(J_calf_3, 2:2, :);
			view(J_calf_4, 2:2, :)]
end

# function Nϕ(model::QuadrupedBasic20, q)
# 	return N_func(model, q)
# end

# function D(model::QuadrupedBasic20, q)
# 	return _P_func(model, q)
# end

# function _P_func(model::QuadrupedBasic20, q)
# 	J_calf_1 = jacobian_2(model, q, body = :calf_1, mode = :ee)
# 	J_calf_2 = jacobian_2(model, q, body = :calf_2, mode = :ee)
# 	J_calf_3 = jacobian_3(model, q, body = :calf_3, mode = :ee)
# 	J_calf_4 = jacobian_3(model, q, body = :calf_4, mode = :ee)
#
# 	return [view(J_calf_1, 1:1, :);
# 			zeros(SVector{11})';
# 			view(J_calf_2, 1:1, :);
# 			zeros(SVector{11})';
# 			view(J_calf_3, 1:1, :);
# 			zeros(SVector{11})';
# 			view(J_calf_4, 1:1, :);
# 			zeros(SVector{11})'
# 			]
# end

function P_func(model::QuadrupedBasic20, q)
	J_calf_1 = jacobian_2(model, q, body = :calf_1, mode = :ee)
	J_calf_2 = jacobian_2(model, q, body = :calf_2, mode = :ee)
	J_calf_3 = jacobian_3(model, q, body = :calf_3, mode = :ee)
	J_calf_4 = jacobian_3(model, q, body = :calf_4, mode = :ee)
	map = [1.0; -1.0]

	return [map * view(J_calf_1, 1:1, :);
			map * view(J_calf_2, 1:1, :);
			map * view(J_calf_3, 1:1, :);
			map * view(J_calf_4, 1:1, :)]
end

function friction_cone(model::QuadrupedBasic20, u)
	λ = view(u, model.idx_λ)
	b = view(u, model.idx_b)
	return @SVector [model.μ * λ[1] - sum(view(b, 1:2)),
					 model.μ * λ[2] - sum(view(b, 3:4)),
					 model.μ * λ[3] - sum(view(b, 5:6)),
					 model.μ * λ[4] - sum(view(b, 7:8))]
end

function maximum_dissipation(model::QuadrupedBasic20, x⁺, u, h)
	q3 = view(x⁺, model.dim.q .+ (1:model.dim.q))
	q2 = view(x⁺, 1:model.dim.q)
	ψ = view(u, model.idx_ψ)
	ψ_stack = [ψ[1] * ones(2); ψ[2] * ones(2); ψ[3] * ones(2); ψ[4] * ones(2)]
	η = view(u, model.idx_η)
	return P_func(model, q3) * (q3 - q2) / h + ψ_stack - η
end

function no_slip(model::QuadrupedBasic20, x⁺, u, h)
	q3 = view(x⁺, model.dim.q .+ (1:model.dim.q))
	q2 = view(x⁺, 1:model.dim.q)
	λ = view(u, model.idx_λ)
	s = view(u, model.idx_s)
	λ_stack = [λ[1]; λ[2]; λ[3]; λ[4]]
	return s[1] - (λ_stack' * _P_func(model, q3) * (q3 - q2) / h)[1]
end

function fd(model::QuadrupedBasic20, x⁺, x, u, w, h, t)
	q3 = view(x⁺, model.dim.q .+ (1:model.dim.q))
	q2⁺ = view(x⁺, 1:model.dim.q)
	q2⁻ = view(x, model.dim.q .+ (1:model.dim.q))
	q1 = view(x, 1:model.dim.q)
	u_ctrl = view(u, model.idx_u)
	λ = view(u, model.idx_λ)
	b = view(u, model.idx_b)

	v = (q3 - q2⁺) / h
	joint_fric = model.joint_friction * v
	joint_fric[1:2] .= 0
    [q2⁺ - q2⁻;
    ((1.0 / h) * (M_func(model, q1) * (SVector{11}(q2⁺) - SVector{11}(q1))
    - M_func(model, q2⁺) * (SVector{11}(q3) - SVector{11}(q2⁺)))
    + h * (transpose(B_func(model, q3)) * SVector{8}(u_ctrl)
	- joint_fric
    + transpose(N_func(model, q3)) * SVector{4}(λ)
    + transpose(P_func(model, q3)) * SVector{8}(b)
    - C_func(model, q3, (q3 - q2⁺) / h)))]
end

# function fd(model::QuadrupedBasic20, x⁺, x, u, w, h, t)
# 	q3 = view(x⁺, model.dim.q .+ (1:model.dim.q))
# 	q2⁺ = view(x⁺, 1:model.dim.q)
# 	q2⁻ = view(x, model.dim.q .+ (1:model.dim.q))
# 	q1 = view(x, 1:model.dim.q)
# 	u_ctrl = view(u, model.idx_u)
# 	λ = view(u, model.idx_λ)
# 	b = view(u, model.idx_b)
# 	h = u[end]
# 	v = (q3 - q2⁺) / h
# 	joint_fric = model.joint_friction * v
# 	joint_fric[1:2] .= 0
#     [q2⁺ - q2⁻;
#     ((1.0 / h) * (M_func(model, q1) * (SVector{11}(q2⁺) - SVector{11}(q1))
#     - M_func(model, q2⁺) * (SVector{11}(q3) - SVector{11}(q2⁺)))
#     + h * (transpose(B_func(model, q3)) * SVector{8}(u_ctrl)
# 	- joint_fric
#     + transpose(N_func(model, q3)) * SVector{4}(λ)
#     + transpose(P_func(model, q3)) * SVector{8}(b)
#     - C_func(model, q3, (q3 - q2⁺) / h)))]
# end

function maximum_dissipation(model::QuadrupedBasic20, x⁺, u, h)
	q3 = view(x⁺, model.dim.q .+ (1:model.dim.q))
	q2 = view(x⁺, 1:model.dim.q)
	ψ = view(u, model.idx_ψ)
	ψ_stack = [ψ[1] * ones(2); ψ[2] * ones(2); ψ[3] * ones(2); ψ[4] * ones(2)]
	η = view(u, model.idx_η)
	h = u[end]
	return P_func(model, q3) * (q3 - q2) / h + ψ_stack - η
end



################################################################################
# Instantiation
################################################################################
# Dimensions
nq = 2 + 5 + 4            # configuration dimension
nu = 4 + 4                # control dimension
# nγ = 4
nc = 4                    # number of contact points
nf = 2                    # number of parameters for friction cone
nb = nc * nf
# ns = 1

# World parameters
h = 0.025
μ = 0.5      # coefficient of friction
g = 9.81     # gravity
joint_friction = 0.1 # coefficient of torque friction at the joints

# ~Unitree A1
# Model parameters
m_torso = 4.713
m_thigh = 1.013
m_leg = 0.166

J_torso = 0.01683
J_thigh = 0.00552
J_leg = 0.00299

l_torso = 0.267
l_thigh = 0.2
l_leg = 0.2

d_torso = 0.0127
d_thigh = 0.00323
d_leg = 0.006435

# n = 2 * nq
# m = nu + nc + nb + nc + nb + ns
# d = 0
#
# idx_u = (1:nu)
# idx_λ = nu .+ (1:nc)
# idx_b = nu + nc .+ (1:nb)
# idx_ψ = nu + nc + nb .+ (1:nc)
# idx_η = nu + nc + nb + nc .+ (1:nb)
# idx_s = nu + nc + nb + nc + nb .+ (1:ns)
#
# qL = -Inf * ones(nq)
# qU = Inf * ones(nq)
#
# uL = -33.5 * ones(nu)
# uU = 33.5 * ones(nu)

# alt2 = 0.0
# alt3 = 0.0
# alt4 = 0.0
# alt1 = 0.0

dim = Dimensions12(nq, nu, nc, nb)
ind = Indices13(nq, nu, nc, nb)
dummy_bas = BaseMethods12()
dummy_dyn = DynamicsMethods13()
dummy_res = ResidualMethods14()
alt = zeros(SizedVector{nc})
quadruped = QuadrupedBasic20(dim, ind, dummy_bas, dummy_dyn, dummy_res, alt,
				h, g, μ, joint_friction,
				# n, m, d,
				# alt1, alt2, alt3, alt4,
				l_torso, d_torso, m_torso, J_torso,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				l_thigh, d_thigh, m_thigh, J_thigh,
				l_leg, d_leg, m_leg, J_leg,
				# qL, qU,
				# uL, uU,
				# nq,
				# nu,
				# nc,
				# nf,
				# nb,
				# ns,
				# idx_u,
				# idx_λ,
				# idx_b,
				# idx_ψ,
				# idx_η,
				# idx_s,
				)


# visualization
function visualize!(vis, model::QuadrupedBasic20, q;
	r = 0.025, Δt = 0.1, name::String="quadruped")
	# default_background!(vis)

	torso = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0, 0.0, 0.0), GeometryBasics.Point3f0(0.0, 0.0, model.l_torso),
		convert(Float32, 0.035))
	setobject!(vis[name*"/torso"], torso,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	thigh_1 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0), GeometryBasics.Point3f0(0.0, 0.0, model.l_thigh1),
		convert(Float32, 0.0175))
	setobject!(vis[name*"/thigh1"], thigh_1,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	calf_1 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0), GeometryBasics.Point3f0(0.0, 0.0, model.l_calf1),
		convert(Float32, 0.0125))
	setobject!(vis[name*"/leg1"], calf_1,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	thigh_2 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0), GeometryBasics.Point3f0(0.0, 0.0, model.l_thigh2),
		convert(Float32, 0.0175))
	setobject!(vis[name*"/thigh2"], thigh_2,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	calf_2 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0), GeometryBasics.Point3f0(0.0, 0.0, model.l_calf2),
		convert(Float32, 0.0125))
	setobject!(vis[name*"/leg2"], calf_2,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	thigh_3 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0), GeometryBasics.Point3f0(0.0, 0.0, model.l_thigh3),
		convert(Float32, 0.0175))
	setobject!(vis[name*"/thigh3"], thigh_3,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	calf_3 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0), GeometryBasics.Point3f0(0.0, 0.0, model.l_calf3),
		convert(Float32, 0.0125))
	setobject!(vis[name*"/leg3"], calf_3,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	thigh_4 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0), GeometryBasics.Point3f0(0.0, 0.0, model.l_thigh4),
		convert(Float32, 0.0175))
	setobject!(vis[name*"/thigh4"], thigh_4,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	calf_4 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0), GeometryBasics.Point3f0(0.0, 0.0, model.l_calf4),
		convert(Float32, 0.0125))
	setobject!(vis[name*"/leg4"], calf_4,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

	hip1 = setobject!(vis[name*"/hip1"], GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
        convert(Float32, 0.035)),
        MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))

	hip2 = setobject!(vis[name*"/hip2"], GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
        convert(Float32, 0.035)),
        MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))

	knee1 = setobject!(vis[name*"/knee1"], GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
        convert(Float32, 0.025)),
        MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))

	knee2 = setobject!(vis[name*"/knee2"], GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		convert(Float32, 0.025)),
		MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))

	knee3 = setobject!(vis[name*"/knee3"], GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		convert(Float32, 0.025)),
		MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))

	knee4 = setobject!(vis[name*"/knee4"], GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
        convert(Float32, 0.025)),
        MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))

	feet1 = setobject!(vis[name*"/feet1"], GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
        convert(Float32, r)),
        MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	feet2 = setobject!(vis[name*"/feet2"], GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	feet3 = setobject!(vis[name*"/feet3"], GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	feet4 = setobject!(vis[name*"/feet4"], GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
        convert(Float32, r)),
        MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	T = length(q)
	p_shift = [0.0, 0.0, r]
	for t = 1:T
		MeshCat.atframe(anim, t) do
			p = [q[t][1]; 0.0; q[t][2]] + p_shift

			k_torso = kinematics_1(model, q[t], body = :torso, mode = :ee)
			p_torso = [k_torso[1], 0.0, k_torso[2]] + p_shift

			k_thigh_1 = kinematics_1(model, q[t], body = :thigh_1, mode = :ee)
			p_thigh_1 = [k_thigh_1[1], 0.0, k_thigh_1[2]] + p_shift

			k_calf_1 = kinematics_2(model, q[t], body = :calf_1, mode = :ee)
			p_calf_1 = [k_calf_1[1], 0.0, k_calf_1[2]] + p_shift

			k_thigh_2 = kinematics_1(model, q[t], body = :thigh_2, mode = :ee)
			p_thigh_2 = [k_thigh_2[1], 0.0, k_thigh_2[2]] + p_shift

			k_calf_2 = kinematics_2(model, q[t], body = :calf_2, mode = :ee)
			p_calf_2 = [k_calf_2[1], 0.0, k_calf_2[2]] + p_shift


			k_thigh_3 = kinematics_2(model, q[t], body = :thigh_3, mode = :ee)
			p_thigh_3 = [k_thigh_3[1], 0.0, k_thigh_3[2]] + p_shift

			k_calf_3 = kinematics_3(model, q[t], body = :calf_3, mode = :ee)
			p_calf_3 = [k_calf_3[1], 0.0, k_calf_3[2]] + p_shift

			k_thigh_4 = kinematics_2(model, q[t], body = :thigh_4, mode = :ee)
			p_thigh_4 = [k_thigh_4[1], 0.0, k_thigh_4[2]] + p_shift

			k_calf_4 = kinematics_3(model, q[t], body = :calf_4, mode = :ee)
			p_calf_4 = [k_calf_4[1], 0.0, k_calf_4[2]] + p_shift

			settransform!(vis[name*"/thigh1"], cable_transform(p, p_thigh_1))
			settransform!(vis[name*"/leg1"], cable_transform(p_thigh_1, p_calf_1))
			settransform!(vis[name*"/thigh2"], cable_transform(p, p_thigh_2))
			settransform!(vis[name*"/leg2"], cable_transform(p_thigh_2, p_calf_2))
			settransform!(vis[name*"/thigh3"], cable_transform(p_torso, p_thigh_3))
			settransform!(vis[name*"/leg3"], cable_transform(p_thigh_3, p_calf_3))
			settransform!(vis[name*"/thigh4"], cable_transform(p_torso, p_thigh_4))
			settransform!(vis[name*"/leg4"], cable_transform(p_thigh_4, p_calf_4))
			settransform!(vis[name*"/torso"], cable_transform(p, p_torso))
			settransform!(vis[name*"/hip1"], MeshCat.Translation(p))
			settransform!(vis[name*"/hip2"], MeshCat.Translation(p_torso))
			settransform!(vis[name*"/knee1"], MeshCat.Translation(p_thigh_1))
			settransform!(vis[name*"/knee2"], MeshCat.Translation(p_thigh_2))
			settransform!(vis[name*"/knee3"], MeshCat.Translation(p_thigh_3))
			settransform!(vis[name*"/knee4"], MeshCat.Translation(p_thigh_4))
			settransform!(vis[name*"/feet1"], MeshCat.Translation(p_calf_1))
			settransform!(vis[name*"/feet2"], MeshCat.Translation(p_calf_2))
			settransform!(vis[name*"/feet3"], MeshCat.Translation(p_calf_3))
			settransform!(vis[name*"/feet4"], MeshCat.Translation(p_calf_4))
		end
	end

	MeshCat.setanimation!(vis, anim)
end

function initial_configuration(model::QuadrupedBasic20, θ)
    q1 = zeros(model.dim.q)
    q1[3] = pi / 2.0
    q1[4] = -θ
    q1[5] = θ
    q1[6] = -θ
    q1[7] = θ
    q1[8] = -θ
    q1[9] = θ
    q1[10] = -θ
    q1[11] = θ
    q1[2] = model.l_thigh1 * cos(q1[4]) + model.l_calf1 * cos(q1[5])
    return q1
end
