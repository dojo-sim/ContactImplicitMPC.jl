"""
    Racecar12
	y
	|              3-------------2 ∡α
	0----> x   rear|             |front
	z              4-------------1 ∡α
"""
mutable struct Racecar12{T} <: ContactModel
    dim::Dimensions

	mb::T # body mass
	Jb::T # body inertia
	lb::T # body length
	wb::T # body width

	mw::T # wheel mass
	Jw::T # wheel inertia
	rw::T # wheel radius

	l0::T # nominal suspension length
	k::T  # suspension stiffness

    μ_world::T  # coefficient of friction
    μ_joint::T  # gravity
	g::T

	base::BaseMethods
	dyn::DynamicsMethods

	joint_friction::SVector
end


# Kinematics
function kinematics_1(model::Racecar12, q; body::Symbol = :wheel1)
	x = q[1]
	y = q[2]
	z = q[3]
	rot = MRP(q[4:6]...)
	p = [x, y, z]

	lb = model.lb
	wb = model.wb

	if body == :wheel1
		return p + rot * [+lb/2, -wb/2, 0.0]
	elseif body == :wheel2
		return p + rot * [+lb/2, +wb/2, 0.0]
	elseif body == :wheel3
		return p + rot * [-lb/2, +wb/2, 0.0]
	elseif body == :wheel4
		return p + rot * [-lb/2, -wb/2, 0.0]
	else
		@error "unknown body."
	end
end

function kinematics_2(model::Racecar12, q; body::Symbol = :wheel_hub1)
	# rot = MRP(q[4:6]...)
	rot = I(3)
	l1 = q[8]
	l2 = q[9]
	l3 = q[10]
	l4 = q[11]

	if body == :wheel_hub1
		p = kinematics_1(model, q; body = :wheel1)
		return p + rot * [0.0, 0.0, -l1]
	elseif body == :wheel_hub2
		p = kinematics_1(model, q; body = :wheel2)
		return p + rot * [0.0, 0.0, -l2]
	elseif body == :wheel_hub3
		p = kinematics_1(model, q; body = :wheel3)
		return p + rot * [0.0, 0.0, -l3]
	elseif body == :wheel_hub4
		p = kinematics_1(model, q; body = :wheel4)
		return p + rot * [0.0, 0.0, -l4]
	else
		@error "unknown body."
	end
end

function kinematics_3(model::Racecar12, q; body::Symbol = :wheel_contact1)
	rw = model.rw

	if body == :wheel_contact1
		p = kinematics_2(model, q; body = :wheel_hub1)
		return p + [0.0, 0.0, -rw]
	elseif body == :wheel_contact2
		p = kinematics_2(model, q; body = :wheel_hub2)
		return p + [0.0, 0.0, -rw]
	elseif body == :wheel_contact3
		p = kinematics_2(model, q; body = :wheel_hub3)
		return p + [0.0, 0.0, -rw]
	elseif body == :wheel_contact4
		p = kinematics_2(model, q; body = :wheel_hub4)
		return p + [0.0, 0.0, -rw]
	else
		@error "unknown body."
	end
end

function kinematics(model::Racecar12, q)
	p1 = kinematics_3(model, q, body = :wheel_contact1)
	p2 = kinematics_3(model, q, body = :wheel_contact2)
	p3 = kinematics_3(model, q, body = :wheel_contact3)
	p4 = kinematics_3(model, q, body = :wheel_contact4)
	SVector{12}([p1; p2; p3; p4])
end

function jacobian_1(model::Racecar12, q; body = :wheel1)
	rot = MRP(q[4:6]...)
	lb = model.lb
	wb = model.wb

	jac = zeros(eltype(q), 3, model.dim.q)
	jac[1, 1] = 1.0
	jac[2, 2] = 1.0
	jac[3, 3] = 1.0

	if body == :wheel1
		jac[:, 4:6] += Rotations.∇rotate(rot, [+lb/2, -wb/2, 0.0])
	elseif body == :wheel2
		jac[:, 4:6] += Rotations.∇rotate(rot, [+lb/2, +wb/2, 0.0])
	elseif body == :wheel3
		jac[:, 4:6] += Rotations.∇rotate(rot, [-lb/2, +wb/2, 0.0])
	elseif body == :wheel4
		jac[:, 4:6] += Rotations.∇rotate(rot, [-lb/2, -wb/2, 0.0])
	else
		@error "incorrect mode specification"
	end
	return jac
end

function jacobian_2(model::Racecar12, q; body = :wheel_hub1)
	# rot = MRP(q[4:6]...)
	rot = I(3)
	l1 = q[8]
	l2 = q[9]
	l3 = q[10]
	l4 = q[11]

	if body == :wheel_hub1
		jac = jacobian_1(model, q; body = :wheel1)
		# jac[:, 4:6] += Rotations.∇rotate(rot, [0.0, 0.0, -l1])
		jac[:, 8] += rot * [0.0, 0.0, -1.0]
	elseif body == :wheel_hub2
		jac = jacobian_1(model, q; body = :wheel2)
		# jac[:, 4:6] += Rotations.∇rotate(rot, [0.0, 0.0, -l2])
		jac[:, 9] += rot * [0.0, 0.0, -1.0]
	elseif body == :wheel_hub3
		jac = jacobian_1(model, q; body = :wheel3)
		# jac[:, 4:6] += Rotations.∇rotate(rot, [0.0, 0.0, -l3])
		jac[:, 10] += rot * [0.0, 0.0, -1.0]
	elseif body == :wheel_hub4
		jac = jacobian_1(model, q; body = :wheel4)
		# jac[:, 4:6] += Rotations.∇rotate(rot, [0.0, 0.0, -l4])
		jac[:, 11] += rot * [0.0, 0.0, -1.0]
	else
		@error "incorrect mode specification"
	end
	return jac
end

function jacobian_3(model::Racecar12, q; body = :wheel_contact1)
	rot = MRP(q[4:6]...)
	# rot = I(3)
	α = q[7]
	θ1 = q[12]
	θ2 = q[13]
	θ3 = q[14]
	θ4 = q[15]
	rw = model.rw
	# the velocity of the contact point attached to front wheel i is
	# rot * rw * ωi [-cos α, -sin α, 0.0]
	# the velocity of the contact point attached to rear wheel i is
	# rot * rw * ωi [-1.0, 0.0, 0.0]
	if body == :wheel_contact1
		jac = jacobian_2(model, q; body = :wheel_hub1)
		jac[:, 12] += rot * rw * [-cos(α), -sin(α), 0.0]
	elseif body == :wheel_contact2
		jac = jacobian_2(model, q; body = :wheel_hub2)
		jac[:, 13] += rot * rw * [-cos(α), -sin(α), 0.0]
	elseif body == :wheel_contact3
		jac = jacobian_2(model, q; body = :wheel_hub3)
		jac[:, 14] += rot * rw * [-1.0, -0.0, 0.0]
	elseif body == :wheel_contact4
		jac = jacobian_2(model, q; body = :wheel_hub4)
		jac[:, 15] += rot * rw * [-1.0, -0.0, 0.0]
	else
		@error "incorrect mode specification"
	end
	return jac
end

function lagrangian(model::Racecar12, q, q̇)
	L = 0.0

	x = q[1]
	y = q[2]
	z = q[3]
	rot = MRP(q[4:6]...)
	α = q[7]
	l1 = q[8]
	l2 = q[9]
	l3 = q[10]
	l4 = q[11]
	θ1 = q[12]
	θ2 = q[13]
	θ3 = q[14]
	θ4 = q[15]

	ẋ = q[1]
	ẏ = q[2]
	ż = q[3]
	ω = q̇[4:6]
	αdot = q̇[7]
	l̇1 = q[8]
	l̇2 = q[9]
	l̇3 = q[10]
	l̇4 = q[11]
	ω1 = q[12]
	ω2 = q[13]
	ω3 = q[14]
	ω4 = q[15]

	# Kinetic energy
		# body
	L += 0.5 * model.mb * (ẋ^2 + ẏ^2 + ż^2)
	L += 0.5 * ω' * Diagonal(model.Jb*ones(3)) * ω
		# wheels
	v1 = jacobian_2(model, q, body = :wheel_hub1) * q̇
	v2 = jacobian_2(model, q, body = :wheel_hub2) * q̇
	v3 = jacobian_2(model, q, body = :wheel_hub3) * q̇
	v4 = jacobian_2(model, q, body = :wheel_hub4) * q̇
	L += 0.5 * model.mw * v1' * v1
	L += 0.5 * model.mw * v2' * v2
	L += 0.5 * model.mw * v3' * v3
	L += 0.5 * model.mw * v4' * v4
	L += 0.5 * model.Jw * (ω1^2 + ω2^2 + ω3^2 + ω4^2)
	L += 0.5 * 2model.Jw * αdot^2

	# Potential energy
	L -= model.mb * model.g * z
	L -= model.mw * model.g * kinematics_2(model, q, body = :wheel_hub1)[3]
	L -= model.mw * model.g * kinematics_2(model, q, body = :wheel_hub2)[3]
	L -= model.mw * model.g * kinematics_2(model, q, body = :wheel_hub3)[3]
	L -= model.mw * model.g * kinematics_2(model, q, body = :wheel_hub4)[3]
	L -= 0.5 * model.k  * (l1 - model.l0)^2
	L -= 0.5 * model.k  * (l2 - model.l0)^2
	L -= 0.5 * model.k  * (l3 - model.l0)^2
	L -= 0.5 * model.k  * (l4 - model.l0)^2
	return L
end

function _dLdq(model::Racecar12, q, q̇)
	# rot = MRP(q[4:6]...)

	dLdq = zeros(eltype(q), model.dim.q)
	dLdq[3] += - (model.mb + 4model.mw) * model.g
	# useless since the opposite wheels compensate each others contribution to potnetial enrgy during rotation.
	# indeed the barycenter of the 4 wheels is always located at the center of the vehicle (using some approx)
	# dLdq[4:6] = Rotations.∇rotate(rot, [0.0, 0.0, 0.0])[3,:]
	dLdq[8]   += - model.mw * model.g
	dLdq[9]   += - model.mw * model.g
	dLdq[10]  += - model.mw * model.g
	dLdq[11]  += - model.mw * model.g
	dLdq[8]  += model.k * (model.l0 - q[8])
	dLdq[9]  += model.k * (model.l0 - q[9])
	dLdq[10] += model.k * (model.l0 - q[10])
	dLdq[11] += model.k * (model.l0 - q[11])
	return dLdq
end

function _dLdq̇(model::Racecar12, q, q̇)
	ẋ = q[1]
	ẏ = q[2]
	ż = q[3]
	ω = q̇[4:6]
	αdot = q̇[7]
	l̇1 = q[8]
	l̇2 = q[9]
	l̇3 = q[10]
	l̇4 = q[11]
	ω1 = q[12]
	ω2 = q[13]
	ω3 = q[14]
	ω4 = q[15]

	J1 = jacobian_2(model, q, body = :wheel_hub1)
	J2 = jacobian_2(model, q, body = :wheel_hub2)
	J3 = jacobian_2(model, q, body = :wheel_hub3)
	J4 = jacobian_2(model, q, body = :wheel_hub4)

	dLdq̇ = zeros(eltype(q), model.dim.q)
	dLdq̇[1] += model.mb * ẋ
	dLdq̇[2] += model.mb * ẏ
	dLdq̇[3] += model.mb * ż
	dLdq̇[4:6] = Diagonal(model.Jb*ones(3)) * ω
	dLdq̇[7] += 2model.Jw * αdot
	dLdq̇[12] += model.Jw * ω1
	dLdq̇[13] += model.Jw * ω2
	dLdq̇[14] += model.Jw * ω3
	dLdq̇[15] += model.Jw * ω4

	dLdq̇ += (J1'*J1 + J2'*J2 + J3'*J3 + J4'*J4) * q̇
	return dLdq̇
end


# Methods
function M_func(model::Racecar12, q) #false
	Diagonal(@SVector [model.mb + 4model.mw,
					   model.mb + 4model.mw,
					   model.mb + 4model.mw,
					   model.Jb + ((model.lb/2)^2 + (model.wb/2)^2)* 4model.mw,
					   model.Jb + ((model.lb/2)^2 + (model.wb/2)^2)* 4model.mw,
					   model.Jb + ((model.lb/2)^2 + (model.wb/2)^2)* 4model.mw,
					   2model.Jw,
					   model.mw,
					   model.mw,
					   model.mw,
					   model.mw,
					   model.Jw,
					   model.Jw,
					   model.Jw,
					   model.Jw,
					   ])
 end

function C_func(model::Racecar12, q, q̇) #false
	nq = model.dim.q
	dLdq̇q = zeros(nq, nq)
	return dLdq̇q * q̇ - _dLdq(model, q, q̇)
end

function ϕ_func(model::Racecar12, env::Environment, q)
	# assumes point of contact is on the suspension axis (flat ground)
	p1 = kinematics_3(model, q, body = :wheel_contact1)
	p2 = kinematics_3(model, q, body = :wheel_contact2)
	p3 = kinematics_3(model, q, body = :wheel_contact3)
	p4 = kinematics_3(model, q, body = :wheel_contact4)
	SVector{4}([
		p1[3] - env.surf(p1[1:2]),
		p2[3] - env.surf(p2[1:2]),
		p3[3] - env.surf(p3[1:2]),
		p4[3] - env.surf(p4[1:2]),
		])
end

function J_func(model::Racecar12, env::Environment, q)
	J1 = jacobian_3(model, q, body = :wheel_contact1)
	J2 = jacobian_3(model, q, body = :wheel_contact2)
	J3 = jacobian_3(model, q, body = :wheel_contact3)
	J4 = jacobian_3(model, q, body = :wheel_contact4)
	return [J1; J2; J3; J4;]
end

function B_func(model::Racecar12, q)
	@SMatrix [ 0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
	 		   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  .01  0.0  0.0  0.0;
			   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  .01  0.0  0.0;
			   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  .01  0.0;
			   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  .01] # /100 scaling
end

function A_func(::Racecar12, q)
	@SMatrix [ 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
	           0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0;
	           0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0]
end

# Parameters
g = 9.81 # gravity
μ_world = 0.20 # coefficient of friction
μ_joint = 10.0

mb = 4.0 # body mass
m_payload = 2.0 # payload mass
Jb = 0.1 # body inertia
lb = 1.0 # body length
wb = 0.5 # body width

mw = 1.0 # wheel mass
Jw = 0.1 # wheel inertia
rw = 0.1 # wheel radius

l0 = 0.15 # suspension length
k = 300.0 # suspension stiffness

# Dimensions
nq = 15
nu = 5
nw = 3
nc = 4
nquat = 0

racecar = Racecar12(Dimensions(nq, nu, nw, nc, nquat),
			   mb, Jb, lb, wb, mw, Jw, rw, l0, k,
			   μ_world, μ_joint, g,
			   BaseMethods(), DynamicsMethods(),
			   SVector{nq}(μ_joint * [zeros(7); ones(4); zeros(4)]))

racecar_payload = Racecar12(Dimensions(nq, nu, nw, nc, nquat),
			   mb + m_payload, Jb, lb, wb, mw, Jw, rw, l0, k,
			   μ_world/2, μ_joint, g,
			   BaseMethods(), DynamicsMethods(),
			   SVector{nq}(μ_joint * [zeros(7); ones(4); zeros(4)]))


# m = zeros(3)
# rt = MRP(m...)
#
# @variables m_[1:3]
# @variables v_[1:3]
# rt_ = MRP(m_...)
# expr = Symbolics.jacobian(rt_*v_, m_)
# Symbolics.jacobian(rt_*[1,0,0], m_)
# Symbolics.jacobian(rt_*[0,1,0], m_)
# Symbolics.jacobian(rt_*[0,0,1], m_)
#
# cod = build_function(expr, m_, v_)[1]
# fct = eval(cod)
#
# m = rand(3)
# v = rand(3)
#
# ∇s = fct(m, v)
# ∇f = ForwardDiff.jacobian(m -> MRP(m...)*v, m)
# ∇a = Rotations.∇rotate(MRP(m...), v)
