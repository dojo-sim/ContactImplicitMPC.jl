"""
    Bicycle
"""
mutable struct Bicycle{T} <: ContactModel
    dim::Dimensions

	mb::T # body mass
	Jb::T # body inertia
    lb::T # body length

	mw::T # wheel mass
	Jw::T # wheel inertia
	rw::T # wheel radius

	l0::T # suspension length
	k::T # suspension stiffness

    μ_world::T  # coefficient of friction
    μ_joint::T  # gravity
	g::T

	base::BaseMethods
	dyn::DynamicsMethods

	joint_friction::SVector
end


# Kinematics
function kinematics_1(model::Bicycle, q; body::Symbol=:front)
	x = q[1]
	z = q[2]
	θ = q[3]
	p = [x, z]
	lb = model.lb
	if body == :front
		return p + lb/2 * [cos(θ), sin(θ)]
	elseif body == :rear
		return p - lb/2 * [cos(θ), sin(θ)]
	else
		@error "unknown body."
	end
end

function kinematics_2(model::Bicycle, q; body::Symbol=:front_hub)
	θ = q[3]
	lf = q[4]
	lr = q[5]
	if body == :front_hub
		p_front = kinematics_1(model, q; body = :front)
		return p_front + lf * [sin(θ), -cos(θ)]
	elseif body == :rear_hub
		p_rear  = kinematics_1(model, q; body = :rear)
		return p_rear  + lr * [sin(θ), -cos(θ)]
	else
		@error "unknown body."
	end
end

function kinematics_3(model::Bicycle, q; body::Symbol=:front_contact)
	rw = model.rw
	if body == :front_contact
		p_front_hub = kinematics_2(model, q; body = :front_hub)
		return p_front_hub + [0.0, -rw]
	elseif body == :rear_contact
		p_rear_hub = kinematics_2(model, q; body = :rear_hub)
		return p_rear_hub + [0.0, -rw]
	else
		@error "unknown body."
	end
end

function kinematics(model::Bicycle, q)
	p_front = kinematics_3(model, q, body = :front_contact)
	p_rear  = kinematics_3(model, q, body = :rear_contact)

	SVector{4}([p_front; p_rear])
end


function jacobian_1(model::Bicycle, q; body = :frame, mode = :com)
	jac = zeros(eltype(q), 2, model.dim.q)
	jac[1, 1] = 1.0
	jac[2, 2] = 1.0
	if mode == :rear
		r = - model.lb / 2
		θ = q[3]
		jac[1, 3] = - r * sin(θ)
		jac[2, 3] = + r * cos(θ)
	elseif mode == :front
		r = + model.lb / 2
		θ = q[3]
		jac[1, 3] = - r * sin(θ)
		jac[2, 3] = + r * cos(θ)
	else
		@error "incorrect mode specification"
	end

	return jac
end

function jacobian_2(model::Bicycle, q; body = :front_hub)
	θ = q[3]

	if body == :front_hub
		jac = jacobian_1(model, q, body = :frame, mode = :front)

		lf = q[4]
		jac[1, 3] += lf * cos(θ)
		jac[2, 3] += lf * sin(θ)
		jac[1, 4] += + sin(θ)
		jac[2, 4] += - cos(θ)
	elseif body == :rear_hub
		jac = jacobian_1(model, q, body = :frame, mode = :rear)

		lr = q[5]
		jac[1, 3] += lr * cos(θ)
		jac[2, 3] += lr * sin(θ)
		jac[1, 5] += + sin(θ)
		jac[2, 5] += - cos(θ)
	else
		@error "incorrect body specification"
	end
	return jac
end

function jacobian_3(model::Bicycle, q; body = :front_contact)
	rw = model.rw

	if body == :front_contact
		jac = jacobian_2(model, q, body = :front_hub)
		jac[1, 6] += rw
		jac[2, 6] += 0.0
	elseif body == :rear_contact
		jac = jacobian_2(model, q, body = :rear_hub)
		jac[1, 7] += rw
		jac[2, 7] += 0.0
	else
		@error "incorrect body specification"
	end
	return jac
end



function lagrangian(model::Bicycle, q, q̇)
	L = 0.0

	x = q[1]
	z = q[2]
	θ = q[3]
	lf = q[4]
	lr = q[5]
	θf = q[6]
	θr = q[7]

	ẋ = q̇[1]
	ż = q̇[2]
	ω = q̇[3]
	l̇f = q̇[4]
	l̇r = q̇[5]
	ωf = q̇[6]
	ωr = q̇[7]

	L += 0.5 * model.mb * (ẋ^2 + ż^2)
	L += 0.5 * model.Jb * ω^2

	L += 0.5 * model.mw * l̇f^2
	L += 0.5 * model.Jw * ωf^2

	L += 0.5 * model.mw * l̇r^2
	L += 0.5 * model.Jw * ωr^2

	z_front = kinematics(model, q, body = :front)[2]
	z_rear  = kinematics(model, q, body = :rear)[2]
	L -= model.mb * model.g * z
	L -= model.mw * model.g * z_front
	L -= model.mw * model.g * z_rear
	L -= 0.5 * model.k  * (lf - model.l0)^2
	L -= 0.5 * model.k  * (lr - model.l0)^2
	return L
end

# Methods
function M_func(model::Bicycle, q) #false
	Diagonal(@SVector [model.mb + 2model.mw,
					   model.mb + 2model.mw,
					   model.Jb + (model.lb/2)^2 * 2model.mw,
					   model.mw,
					   model.mw,
					   model.Jw,
					   model.Jw,
					   ])
 end

function C_func(model::Bicycle, q, q̇) #false
	@SVector [0.0,
			  0.0,
			  0.0,
			  0.0,
			  0.0,
			  0.0,
			  0.0,
			  ]
end

function ϕ_func(model::Bicycle, env::Environment, q)
	# assumes point of contact is on the suspension axis (flat ground)
	p_front = kinematics_3(model, q, body = :front_contact)
	p_rear  = kinematics_3(model, q, body = :rear_contact)
	SVector{2}([
		p_front[2] - env.surf(p_front[1]),
		p_rear[2]  - env.surf(p_rear[1]),
		])
end

function J_func(model::Bicycle, env::Environment, q)
	J_front = jacobian_3(model, q, body = :front_contact)
	J_rear  = jacobian_3(model, q, body = :rear_contact)

	return [J_front;
			J_rear;]
end


function B_func(model::Bicycle, q)
	@SMatrix [ 0.0  0.0  0.0  1.0  0.0  0.0      0.0;
	 		   0.0  0.0  0.0  0.0  1.0  0.0      0.0;
			   0.0  0.0  0.0  0.0  0.0  1.0/100  0.0;
			   0.0  0.0  0.0  0.0  0.0  0.0      1.0/100] # /100 scaling
end

function A_func(::Bicycle, q)
	@SMatrix [ 1.0  0.0  0.0  0.0  0.0  0.0  0.0;
			   0.0  1.0  0.0  0.0  0.0  0.0  0.0]
end

# Parameters
g = 9.81 # gravity
μ_world = 0.1 # coefficient of friction
μ_joint = 10.0

mb = 1.0 # body mass
Jb = 0.1 # body inertia
lb = 1.0 # body length

mw = 1.0 # wheel mass
Jw = 0.1 # wheel inertia
rw = 0.1 # wheel radius

l0 = 0.5 # suspension length
k = 300.0 # suspension stiffness

# Dimensions
nq = 7
nu = 4
nw = 2
nc = 2
nquat = 0

bicycle = Bicycle(Dimensions(nq, nu, nw, nc, nquat),
			   mb, Jb, lb, mw, Jw, rw, l0, k,
			   μ_world, μ_joint, g,
			   BaseMethods(), DynamicsMethods(),
			   SVector{nq}(μ_joint * [0.0; 0.0; 0.0; 1.0; 1.0; 0.0; 0.0]))
