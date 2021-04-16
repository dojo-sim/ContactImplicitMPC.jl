function com_func(model::Biped5, q)
	m = model.m_torso + model.m_thigh1 + model.m_thigh2 + model.m_calf1 + model.m_calf2
	p_torso  = kinematics_1(model, q, body=:torso, mode=:com)
	p_thigh1 = kinematics_1(model, q, body=:thigh_1, mode=:com)
	p_thigh2 = kinematics_1(model, q, body=:thigh_2, mode=:com)
	p_calf1  = kinematics_2(model, q, body=:calf_1, mode=:com)
	p_calf2  = kinematics_2(model, q, body=:calf_2, mode=:com)
	p_com = model.m_torso  * p_torso +
			model.m_thigh1 * p_thigh1 +
			model.m_thigh2 * p_thigh2 +
			model.m_calf1  * p_thigh1 +
			model.m_calf2  * p_thigh2
	p_com ./= m
	return p_com
end

function com_func(model::PinnedBiped513, q)
	m = model.m_torso + model.m_thigh1 + model.m_thigh2 + model.m_calf1 + model.m_calf2
	p_calf1  = kinematics__1(model, q, body=:calf_1, mode=:com)
	p_thigh1 = kinematics__2(model, q, body=:thigh_1, mode=:com)
	p_torso  = kinematics__3(model, q, body=:torso, mode=:com)
	p_thigh2 = kinematics__3(model, q, body=:thigh_2, mode=:com)
	p_calf2  = kinematics__4(model, q, body=:calf_2, mode=:com)
	p_com = model.m_torso  * p_torso +
			model.m_thigh1 * p_thigh1 +
			model.m_thigh2 * p_thigh2 +
			model.m_calf1  * p_thigh1 +
			model.m_calf2  * p_thigh2
	p_com ./= m
	return p_com
end

function h_func(p_com, p_stance, p_swing, p_hip, p_torso)

	v_swing  = cast3d(p_swing  - p_com)
	v_stance = cast3d(p_stance - p_com)
	v_trunk  = cast3d(p_torso - p_hip)

	r_com = norm(p_com - p_stance)
	r_foot = norm(p_com - p_swing)
	θ_foot = atan(cross(v_swing, v_stance)[2], v_swing'*v_stance)
	θ_trunk = atan(cross(v_trunk, -v_stance)[2], v_trunk'*(-v_stance))
	h = [r_com, r_foot, θ_foot, θ_trunk]
	return h
end

function h_func(model::Biped5, q)
	p_com = com_func(model, q)
	p_stance  = kinematics_2(model, q, body=:calf_1, mode=:ee)
	p_swing   = kinematics_2(model, q, body=:calf_2, mode=:ee)
	p_hip     = [q[1], q[2]]
	p_torso   = kinematics_1(model, q, body=:torso, mode=:ee)
	return h_func(p_com, p_stance, p_swing, p_hip, p_torso)
end

function h_func(model::PinnedBiped513, q)
	p_com = com_func(model, q)
	p_stance  = [0.0, 0.0]
	p_swing   = kinematics__4(model, q, body=:calf_2, mode=:ee)
	p_hip     = kinematics__2(model, q, body=:thigh_1, mode=:ee)
	p_torso   = kinematics__3(model, q, body=:torso, mode=:ee)
	return h_func(p_com, p_stance, p_swing, p_hip, p_torso)
end

function hd_func(model::PinnedBiped513, q, qd)
	h_(q) = h_func(model, q)
	∇qh = ForwardDiff.jacobian(h_, q)
	hd = ∇qh*qd
	return hd
end

function θcom_func(p_com, p_stance)
	vx = cast3d([1.0, 0.0])
	v_stance = cast3d(p_com - p_stance)
	θcom = atan(cross(v_stance, vx)[2], v_stance'*(vx))
	return θcom
end

function θcom_func(model::Biped5, q)
	p_com = com_func(model, q)
	p_stance  = kinematics_2(model, q, body=:calf_1, mode=:ee)
	@show p_com
	@show p_stance

	return θcom_func(p_com, p_stance)
end

function θcom_func(model::PinnedBiped513, q)
	p_com = com_func(model, q)
	p_stance  = [0.0, 0.0]
	return θcom_func(p_com, p_stance)
end





function f_func(model::PinnedBiped513, x)
	# M(q)*qdd + C(q, qd) = B'*u + A*w + J*λ
	# M(q)*qdd + C(q, qd) = B'*u
	# xd = [qd ] = [qd              ] = [qd       ] + [0       ] * [u]
	#    = [qdd] = [inv(M)*(B'*u - C)] = [-inv(M)*C] + [inv(M)*B']   [ ]
	# xd = f(x) + g(x)*u
	# f(x) = [qd       ] ∈ R{nx}
	#        [-inv(M)*C]
	nq = model.dim.q
	nx = 2nq
	q  = x[1:nq]
	qd = x[nq .+ (1:nq)]

	M = M_func(model, q)
	C = _C_func(model, q, qd)
	f = [qd; - M\C]
	return f
end

function g_func(model::PinnedBiped513, x)
	# M(q)*qdd + C(q, qd) = B*u + A*w + J*λ
	# M(q)*qdd + C(q, qd) = B*u
	# xd = [qd ] = [qd              ] = [qd       ] + [0       ] * [u]
	#    = [qdd] = [inv(M)*(B*u - C)] = [-inv(M)*C] + [inv(M)*B]   [ ]
	# xd = f(x) + g(x)*u
	# g(x) = [0{nq,nu}] ∈ R{nx,nu}
	#        [inv(M)*B]
	nq = model.dim.q
	nx = 2nq
	nh = 4
	q  = x[1:nq]
	qd = x[nq .+ (1:nq)]

	M = M_func(model, q)
	B = B_func(model, q)
	g = zeros(nx, nh)
	# g[nq .+ (1:nq), 1:nh] = (M\transpose(B)#)[:,1 .+ (1:nh)]
	g[nq .+ (1:nq), 1:nh] = M\transpose(B[[1,2,4,5], :])
	return g
end

function ẋ_func(model::PinnedBiped513, x, u)
	ẋ = f_func(model, x) + g_func(model, x)*u
	return ẋ
end

function ∇fu_func(model::PinnedBiped513, x)
	nh = 4

	f_(x) = f_func(model, x)
	h_(x) = h_func(model, x)
	hi_ = [x -> h_func(model, x)[i] for i=1:nh]

	f = f_func(model, x)
	dfdx = ForwardDiff.jacobian(f_, x)
	dhdx = ForwardDiff.jacobian(h_, x)
	∇fu = dhdx * dfdx
	# @show size(dhdx)
	# @show size(dfdx)
	# @show size(∇fu)
	# @show size(f)
	for i = 1:nh
		# @show size(∇fu[i:1,:])
		# @show size(ForwardDiff.hessian(hi_[i], x))
		# @show size(f'*ForwardDiff.hessian(hi_[i], x))
		∇fu[i:i,:] += f'*ForwardDiff.hessian(hi_[i], x)
	end
	return ∇fu
end

function Afu_func(model::PinnedBiped513, x; ∇fu=∇fu_func(model, x))
	g = g_func(model, x)
	A = ∇fu * g
	return A
end

function bfu_func(model::PinnedBiped513, x; ∇fu=∇fu_func(model, x))
	f = f_func(model, x)
	b = ∇fu * f
	return b
end




function dynamics(model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)

	# evalutate at midpoint
	qm1 = 0.5 * (q0 + q1)
    vm1 = (q1 - q0) / h[1]
    qm2 = 0.5 * (q1 + q2)
    vm2 = (q2 - q1) / h[1]

	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	nc = model.dim.c
	nb = model.dim.b
	nf = Int(nb / nc)
	ne = dim(model.env)
	k = kinematics(model, q2)
	# λ1 = vcat([transpose(rotation(model.env, q2)) * [friction_mapping(model.env) * b1[(i-1) * nf .+ (1:nf)]; γ1[i]] for i = 1:nc]...) # TODO: make efficient
	λ1 = vcat([transpose(rotation(model.env, k[(i-1) * ne .+ (1:ne)])) * [friction_mapping(model.env) * b1[(i-1) * nf .+ (1:nf)]; γ1[i]] for i = 1:nc]...) # TODO: make efficient

	return (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
		+ transpose(B_fast(model, qm2)) * u1
		+ transpose(A_fast(model, qm2)) * w1
		+ transpose(J_fast(model, q2)) * λ1
		- h[1] * model.joint_friction .* vm2)
end


function pin_state(q)
	qp = zeros(5)
	qp[1] = q[5]
	qp[2] = q[4]
	qp[3] = q[3]
	qp[4] = q[6]
	qp[5] = q[7]
	return qp
end

function unpin_state(qp)
	q = zeros(7)
	q[1] = 0.0
	q[2] = 0.0
	q[3] = qp[3]
	q[4] = qp[2]
	q[5] = qp[1]
	q[6] = qp[4]
	q[7] = qp[5]
	return q
end



function s_func(model::PinnedBiped513, q;
		q_ini=[-0.228, 0.228, -0.1, -0.1, -0.3],
		q_mid=[-0.35,  0.2,   -0.1,  0.3, -0.4],
		q_end=[-0.3, -0.1, -0.1, 0.228, -0.228],
		)
	θcom_ini = θcom_func(model, q_ini)
	θcom_end = θcom_func(model, q_end)

	θcom = θcom_func(model, q)
	s = (θcom - θcom_ini)/(θcom_end - θcom_ini)
	# s = clamp(s, 0.0, 1.0)
	return s
end

function sd_func(model::PinnedBiped513, q, qd;
		q_ini=[-0.228, 0.228, -0.1, -0.1, -0.3],
		q_mid=[-0.35,  0.2,   -0.1,  0.3, -0.4],
		q_end=[-0.3, -0.1, -0.1, 0.228, -0.228],
		)

	θcom_(q) = θcom_func(model, q)
	∇qθcom = ForwardDiff.gradient(θcom_, q)
	θcomd = ∇qθcom'*qd

	θcom_ini = θcom_func(model, q_ini)
	θcom_end = θcom_func(model, q_end)
	sd = θcomd/(θcom_end - θcom_ini)

	return sd
end

function p_func(model::PinnedBiped513, s;
		q_ini=[-0.228, 0.228, -0.1, -0.1, -0.3],
		q_mid=[-0.35,  0.2,   -0.1,  0.3, -0.4],
		q_end=[-0.3, -0.1, -0.1, 0.228, -0.228],
		)
	nh = 4
	θcom_ini = θcom_func(model, q_ini)
	θcom_mid = θcom_func(model, q_mid)
	θcom_end = θcom_func(model, q_end)

	h_ini = h_func(model, q_ini)
	h_mid = h_func(model, q_mid)
	h_end = h_func(model, q_end)

	s_ini = 0.0
	s_mid = (θcom_mid - θcom_ini)/(θcom_end - θcom_ini)
	s_end = 1.0

	p = zeros(nh)
	if s < s_mid
		α = (s - s_ini)/(s_mid - s_ini)
		p = h_mid*α + (1-α)*h_ini
	elseif s >= 0*s+s_mid
		α = (s - s_mid)/(s_end - s_mid)
		p = h_end*α + (1-α)*h_mid
	end
	return p
end

function pd_func(model::PinnedBiped513, s;
		q_ini=[-0.228, 0.228, -0.1, -0.1, -0.3],
		q_mid=[-0.35,  0.2,   -0.1,  0.3, -0.4],
		q_end=[-0.3, -0.1, -0.1, 0.228, -0.228],
		)
	nh = 4
	θcom_ini = θcom_func(model, q_ini)
	θcom_mid = θcom_func(model, q_mid)
	θcom_end = θcom_func(model, q_end)

	h_ini = h_func(model, q_ini)
	h_mid = h_func(model, q_mid)
	h_end = h_func(model, q_end)

	s_ini = 0.0
	s_mid = (θcom_mid - θcom_ini)/(θcom_end - θcom_ini)
	s_end = 1.0

	pd = zeros(nh)
	if s < s_mid
		pd = (h_mid - h_ini)/(s_mid - s_ini)
	elseif s >= s_mid
		pd = (h_end - h_mid)/(s_end - s_mid)
	end
	return pd
end

function pdd_func(model::PinnedBiped513, s;
		q_ini=[-0.228, 0.228, -0.1, -0.1, -0.3],
		q_mid=[-0.35,  0.2,   -0.1,  0.3, -0.4],
		q_end=[-0.3, -0.1, -0.1, 0.228, -0.228],
		)
	nh = 4
	pdd = zeros(nh)
	return pdd
end


function pin_state(q)
	qp = zeros(5)
	qp[1] = q[5]
	qp[2] = q[4]
	qp[3] = q[3]
	qp[4] = q[6]
	qp[5] = q[7]
	return qp
end

function unpin_state(qp)
	q = zeros(7)
	q[1] = 0.0
	q[2] = 0.0
	q[3] = qp[3]
	q[4] = qp[2]
	q[5] = qp[1]
	q[6] = qp[4]
	q[7] = qp[5]
	return q
end






#
#
# function r_com_func(model::PinnedBiped513, q; body=:toe_1)
# 	com = com_func(model, q)
# 	p_toe = toe_func(model, q, body=body)
# 	r_com = norm(com - p_toe)
# 	return r_com
# end
#
# function opposite_toe(body)
# 	if body == :toe_1
# 		return :toe_2
# 	elseif body == :toe_2
# 		return :toe_1
# 	else
# 		error("incorrect body specification")
# 	end
# end
#
# function θ_foot_func(model::PinnedBiped513, q; stance_body=:toe_1)
# 	com = com_func(model, q)
# 	p_swing_toe = toe_func(model, q, body=opposite_toe(stance_body))
# 	p_stance_toe = toe_func(model, q, body=stance_body)
# 	v_swing  = cast3d(p_swing_toe  - com)
# 	v_stance = cast3d(p_stance_toe - com)
# 	θ_foot = atan(cross(v_swing, v_stance)[2], v_swing'*v_stance)
# 	return θ_foot
# end
#
# function θ_trunk_func(model::PinnedBiped513, q; stance_body=:toe_1)
# 	com = com_func(model, q)
# 	p_stance_toe = toe_func(model, q, body=stance_body)
# 	v_stance = cast3d(com - p_stance_toe)
# 	p_hip = kinematics__2(model, q, body=:thigh_1, mode=:ee)
# 	p_torso = kinematics__3(model, q, body=:torso, mode=:ee)
# 	v_trunk  = cast3d(p_torso - p_hip)
# 	θ_foot = atan(cross(v_trunk, v_stance)[2], v_trunk'*v_stance)
# 	return θ_foot
# end
#
# function h_func(model::PinnedBiped513, q; stance_body=:toe_1)
# 	swing_body = opposite_toe(stance_body)
# 	r_com = r_com_func(model, q, body=stance_body)
# 	r_foot = r_com_func(model, q, body=swing_body)
# 	θ_foot = θ_foot_func(model, q, stance_body=stance_body)
# 	θ_trunk = θ_trunk_func(model, q, stance_body=stance_body)
# 	return [r_com, r_foot, θ_foot, θ_trunk]
# end
#
# function toe_func(model::Biped5, q; body=:toe_1)
# 	if body == :toe_1
# 		p_toe = kinematics_2(model, q, body=:calf_1, mode=:ee)
# 	elseif body == :toe_2
# 		p_toe = kinematics_2(model, q, body=:calf_2, mode=:ee)
# 	else
# 		error("incorrect body specification")
# 	end
# 	return p_toe
# end
#
