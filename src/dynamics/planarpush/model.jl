"""
    PlanarPush
        s = (x, y, t, xp, yp)
            x  - x position of the object
            y  - y position of the object
            z  - z position of the object
            t  - object orientation along the Z axis
            xp - x position of the pusher
            yp - y position of the pusher
            zp - z position of the pusher
"""
#                                                  o  o
#           XXX    XXX                          o        o
#      XXX              XXX                    o     .(xp,yp)
#                                              o          o
#  XXX                      XXX          -      o        o
#                                -   \             o  o
# XXX                   -     XXX     |α
#                .(x,y)---------------------------------------------
# XXX                         XXX
#
# XXX                        XXX
#
#    XXX                  XXX
#
#          XXX     XXX

mutable struct PlanarPush13{T} <: ContactDynamicsModel
    dim::Dimensions

	μ_world::T  # coefficient of friction
	μ_joint::T  # gravity
	g::T

    m::T  # mass of object
	J::T  # inertia of object
    mp::T # mass of pusher

	r::T  # radius of object
	rp::T # radius of pusher

	alt

	base::BaseMethods
	dyn::DynamicsMethods
	con::ContactMethods
	res::ResidualMethods
	linearized::ResidualMethods

	spa::SparseStructure

	joint_friction::SVector

	env::Environment
end

# function kinematics_1(model::PlanarPush13, q; body = :torso, mode = :ee)
# 	x = q[1]
# 	z = q[2]
#
# 	if body == :torso
# 		l = model.l_torso
# 		d = model.d_torso
# 		θ = q[3]
# 		if mode == :ee
# 			return [x - l * sin(θ); z + l * cos(θ)]
# 		elseif mode == :com
# 			return [x - d * sin(θ); z + d * cos(θ)]
# 		end
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
# function jacobian_1(model::PlanarPush13, q; body = :torso, mode = :ee)
# 	jac = zeros(eltype(q), 2, model.dim.q)
# 	jac[1, 1] = 1.0
# 	jac[2, 2] = 1.0
# 	if body == :torso
# 		r = mode == :ee ? model.l_torso : model.d_torso
# 		θ = q[3]
# 		jac[1, 3] = -r * cos(θ)
# 		jac[2, 3] = -r * sin(θ)
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
# function kinematics_2(model::PlanarPush13, q; body = :calf_1, mode = :ee)
#
# 	if body == :calf_1
# 		p = kinematics_1(model, q, body = :thigh_1, mode = :ee)
#
# 		θb = q[5]
#
# 		lb = model.l_calf1
# 		db = model.d_calf1
# 	elseif body == :calf_2
# 		p = kinematics_1(model, q, body = :thigh_2, mode = :ee)
#
# 		θb = q[7]
#
# 		lb = model.l_calf2
# 		db = model.d_calf2
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
# function jacobian_2(model::PlanarPush13, q; body = :calf_1, mode = :ee)
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
# 	else
# 		@error "incorrect body specification"
# 	end
#
# 	return jac
# end
#
# function kinematics_3(model::PlanarPush13, q; body = :foot_1, mode = :ee)
#
# 	if body == :foot_1
# 		p = kinematics_2(model, q, body = :calf_1, mode = :ee)
#
# 		θb = q[8]
#
# 		lb = model.l_foot1
# 		db = model.d_foot1
# 		cb = 0.5 * (model.l_foot1 - model.d_foot1)
# 	elseif body == :foot_2
# 		p = kinematics_2(model, q, body = :calf_2, mode = :ee)
#
# 		θb = q[9]
#
# 		lb = model.l_foot2
# 		db = model.d_foot2
# 		cb = 0.5 * (model.l_foot2 - model.d_foot2)
# 	else
# 		@error "incorrect body specification"
# 	end
#
# 	if mode == :toe
# 		return p + [lb * sin(θb); -1.0 * lb * cos(θb)]
# 	elseif mode == :heel
# 		return p + [-db * sin(θb); 1.0 * db * cos(θb)]
# 	elseif mode == :com
# 		return p + [cb * sin(θb); -1.0 * cb * cos(θb)]
# 	else
# 		@error "incorrect mode specification"
# 	end
# end
#
# function jacobian_3(model::PlanarPush13, q; body = :foot_1, mode = :ee)
#
# 	if body == :foot_1
# 		jac = jacobian_2(model, q, body = :calf_1, mode = :ee)
#
# 		θb = q[8]
#
# 		if mode == :toe
# 			r = model.l_foot1
# 		elseif mode == :heel
# 			r = -1.0 * model.d_foot1
# 		elseif mode == :com
# 			r = 0.5 * (model.l_foot1 - model.d_foot1)
# 		else
# 			@error "incorrect mode specification"
# 		end
#
# 		jac[1, 8] += r * cos(θb)
# 		jac[2, 8] += r * sin(θb)
#
# 	elseif body == :foot_2
# 		jac = jacobian_2(model, q, body = :calf_2, mode = :ee)
#
# 		θb = q[9]
#
# 		if mode == :toe
# 			r = model.l_foot2
# 		elseif mode == :heel
# 			r = -1.0 * model.d_foot2
# 		elseif mode == :com
# 			r = 0.5 * (model.l_foot2 - model.d_foot2)
# 		else
# 			@error "incorrect mode specification"
# 		end
# 		jac[1, 9] += r * cos(θb)
# 		jac[2, 9] += r * sin(θb)
#
# 	else
# 		@error "incorrect body specification"
# 	end
#
# 	return jac
# end

# Lagrangian

function lagrangian(model::PlanarPush13, q, q̇)
	v_object = [q̇[1], q̇[2], q̇[3]]
	v_pusher = [q̇[5], q̇[6], q̇[7]]

	L = 0.0
	# kinetic energy
	L += 0.5 * model.m * transpose(v_object) * v_object
	L += 0.5 * model.mp * transpose(v_pusher) * v_pusher
	L += 0.5 * model.J * q̇[4]^2
	# potential energy
	L -= model.m * model.g * q[3]
	L -= model.mp * model.g * q[7]
	return L
end

function _dLdq(model::PlanarPush13, q, q̇)
	Lq(x) = lagrangian(model, x, q̇)
	ForwardDiff.gradient(Lq, q)
end

function _dLdq̇(model::PlanarPush13, q, q̇)
	Lq̇(x) = lagrangian(model, q, x)
	ForwardDiff.gradient(Lq̇, q̇)
end

function _C_func(model::PlanarPush13, q, q̇)
	tmp_q(z) = _dLdq̇(model, z, q̇)
	tmp_q̇(z) = _dLdq̇(model, q, z)

	ForwardDiff.jacobian(tmp_q, q) * q̇ - _dLdq(model, q, q̇)
end

function kinematics_1(model::PlanarPush13, q; body = :floor)
	# Contact 1 is between object and floor: :floor
	# Contact 2 is between object and pusher: :object, :pusher
	x = [q[1], q[2], 0.0]
	xp = [q[5], q[6], 0.0]
	Δ = x - xp
	Δ ./= norm(Δ)
	if body == :floor
		return x
	elseif body == :object
		return x - Δ*model.r
	elseif body == :pusher
		return xp + Δ*model.rp
	else
		@error "incorrect body specification"
	end
end

function kinematics(model::PlanarPush13, q)
	p_floor = kinematics_1(model, q; body = :floor)
	p_object = kinematics_1(model, q; body = :object)
	p_pusher= kinematics_1(model, q; body = :pusher)
	SVector{9}([p_floor; p_object; p_pusher])
end

# Methods
function M_func(model::PlanarPush13, q)
	Diagonal(@SVector [model.m, model.m, model.m,
					   model.J,
					   model.mp, model.mp, model.mp,])
end

function ϕ_func(model::PlanarPush13, q)
	x = [q[1], q[2], 0.0]
	xp = [q[5], q[6], 0.0]
	Δ = x - xp
	SVector{3}([q[3] - model.env.surf(x[1:2]),
				norm(Δ) - (model.r + model.rp),
				norm(Δ) - (model.r + model.rp)])
end

function B_func(model::PlanarPush13, q)
	@SMatrix [0.0  0.0  0.0  0.0  1.0  0.0  0.0;
			  0.0  0.0  0.0  0.0  0.0  1.0  0.0]
end

function A_func(model::PlanarPush13, q)
	@SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 1.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 1.0 0.0 0.0 0.0 0.0]
end

function contact_angle(model::PlanarPush13, q)
	# Angle of the contact normal
	x = q[1]
	y = q[2]
	xp = q[5]
	yp = q[6]
	α = atan(yp - y, xp - x)
	return α
end

function J_func(model::PlanarPush13, q)
	# Jacobian velocity in world frame of contact points ATTACHED to their respective links.
	J_floor = [1.0  0.0  0.0  0.0  0.0  0.0  0.0;
			   0.0  1.0  0.0  0.0  0.0  0.0  0.0;
			   0.0  0.0  0.0  0.0  0.0  0.0  0.0]
	# p_object = [x + r cos(α) , y + sin(α)]
	# v_object = [ẋ - r sin(α) θd, ẏ + r cos(α) θd]
	r = model.r
	α = contact_angle(model, q)
    J_object = [1.0  0.0  0.0 -r*sin(α) 0.0  0.0  0.0;
			    0.0  1.0  0.0  r*cos(α) 0.0  0.0  0.0;
			    0.0  0.0  0.0  0.0      0.0  0.0  0.0]
	J_pusher = [0.0  0.0  0.0  0.0  1.0  0.0  0.0;
			    0.0  0.0  0.0  0.0  0.0  1.0  0.0;
			    0.0  0.0  0.0  0.0  0.0  0.0  0.0]
	return [J_floor;
			J_object;
			J_pusher] # (nc*np) x nq  = 9x7
end

function contact_forces(model::PlanarPush13, γ1, b1, q2, k)
	# Express the contact forces in the world frame.
	# returns λ1
	# which will be used in dynamics: transpose(J_fast(model, q2)) * λ1
	m = friction_mapping(model.env)
	# get rotation
	α = contact_angle(model, q)
	cα = cos(α)
	sα = sin(α)
	R = [-sα  0.0  cα; # rotation matrix for contact frame to world frame
		  cα  0.0  sα;
		 0.0  1.0  0.0]
	# Express in world frame
	λ_floor = [m * b1[1:4]; γ1[1]]
	λ_object_c = [m * b1[5:8]; γ1[2]]
	λ_pusher_c = [m * b1[9:12]; γ1[3]]
	λ_object =  R * λ_object_c # rotate in world frame
	λ_pusher = -R * λ_pusher_c # action reaction principle
	SVector{9}([λ_floor; λ_object; λ_pusher])
end


# function contact_forces(model::Hopper3D, γ1, b1, q2, k)
# 	# k = kinematics(model, q2)
# 	m = friction_mapping(model.env)
#
# 	SVector{3}(transpose(rotation(model.env, k)) * [m * b1; γ1])
# end


function velocity_stack(model::PlanarPush13, q1, q2, k, h)
	# k = kinematics(model, q2)
	v = J_func(model, q2) * (q2 - q1) / h[1]

	# In the world frame
	v_floor = rotation(model.env, k) * v[1:3]
	v_object = rotation(model.env, k) * v[4:6]
	v_pusher = rotation(model.env, k) * v[7:9]


	SVector{12}([v_floor[1:2];  -v_floor[1:2];
				v_object[1:2]; -v_object[1:2];
				v_pusher[1:2]; -v_pusher[1:2];])
end

# Dimensions
nq = 3 + 1 + 3            # configuration dimension
nu = 2                    # control dimension
nc = 3                    # number of contact points
nf = 4                    # number of parameters for friction cone
nb = nc * nf              # number of friction parameters
nw = 3                    # disturbance dimension

# World parameters
μ_world = 0.9      # coefficient of friction
μ_joint = 0.0
g = 9.81     # gravity

# Model parameters
m = 1.0
J = 0.0025
mp = 0.1
r = 0.1
rp = 0.01


planarpush = PlanarPush13(Dimensions(nq, nu, nw, nc, nb),
			  g, μ_world, μ_joint,
			  m, J, mp, r, rp,
			  zeros(nc),
			  BaseMethods(), DynamicsMethods(), ContactMethods(),
			  ResidualMethods(), ResidualMethods(),
			  SparseStructure(spzeros(0, 0), spzeros(0, 0)),
			  SVector{nq}([zeros(3); 0.0 * μ_joint * ones(nq - 3)]),
			  environment_3D_flat())

@variables q[1:7]
@variables q̇[1:7]

L = lagrangian(planarpush, q, q̇)
ddLdq̇q̇ = Symbolics.sparsehessian(L, q̇, simplify=true)
