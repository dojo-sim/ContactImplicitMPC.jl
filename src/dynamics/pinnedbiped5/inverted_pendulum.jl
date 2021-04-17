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


function unpin_state(model::Biped5, qp; body=:toe_1)
	# qp = [q_calf_front, q_thigh_front, q_torso, q_thigh_rear, q_calf_rear]
	# q  = [x, y, q_torso, q_thigh_1, q_calf_1, q_thigh_2, q_calf_2]
	q = zeros(7)
	if body == :toe_1
		q[1] = 0.0
		q[2] = 0.0
		q[3] = qp[3]
		q[4] = qp[2]
		q[5] = qp[1]
		q[6] = qp[4]
		q[7] = qp[5]
		q[1:2] = -kinematics_2(model, q, body=:calf_1, mode=:ee)
	elseif body == :toe_2
		q[1] = 0.0
		q[2] = 0.0
		q[3] = qp[3]
		q[4] = qp[4]
		q[5] = qp[5]
		q[6] = qp[2]
		q[7] = qp[1]
		q[1:2] = -kinematics_2(model, q, body=:calf_2, mode=:ee)
	else
		@error "incorrect body specification"
	end
	return q
end

function pin_state(q; body=:toe_1)
	qp = zeros(5)
	if body == :toe_1
		qp[1] = q[5]
		qp[2] = q[4]
		qp[3] = q[3]
		qp[4] = q[6]
		qp[5] = q[7]
	elseif body == :toe_2
		qp[1] = q[7]
		qp[2] = q[6]
		qp[3] = q[3]
		qp[4] = q[4]
		qp[5] = q[5]
	else
		@error "incorrect body specification"
	end
	return qp
end

function unpin_control(up; body=:toe_1)
	u = zeros(5)
	if body == :toe_1
		u[1] = 0.0   # u_torso
		u[2] = up[2] # u_thigh_1
		u[3] = up[1] # u_calf_1
		u[4] = up[3] # u_thigh_2
		u[5] = up[4] # u_calf_2
	elseif body == :toe_2
		u[1] = 0.0   # u_torso
		u[2] = up[3] # u_thigh_1
		u[3] = up[4] # u_calf_1
		u[4] = up[2] # u_thigh_2
		u[5] = up[1] # u_calf_2
	else
		@error "incorrect body specification"
	end
	return u
end





qp = [-0.228, 0.228, -0.1, -0.1, -0.3]
q1 = unpin_state(model, qp, body=:toe_1)
q2 = unpin_state(model, qp, body=:toe_2)

q1p = pin_state(unpin_state(model, qp, body=:toe_1), body=:toe_1)
q2p = pin_state(unpin_state(model, qp, body=:toe_2), body=:toe_2)

build_robot!(vis, pinnedmodel, r=0.004)
build_robot!(vis, model, r=0.004)

set_robot!(vis, pinnedmodel, qp, r=0.004)
set_robot!(vis, pinnedmodel, q1p, r=0.004)
set_robot!(vis, pinnedmodel, q2p, r=0.004)
set_robot!(vis, model, q1, r=0.004)
set_robot!(vis, model, q2, r=0.004)





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
