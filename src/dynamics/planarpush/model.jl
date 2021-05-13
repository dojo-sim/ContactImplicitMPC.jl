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
	# L -= model.mp * model.g * q[7] # the pusher floats in the air, not subject to gravity
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
	p_pusher = kinematics_1(model, q; body = :pusher)
	# SVector{9}([p_floor; p_object; p_pusher])

	SVector{6}([p_floor; p_object])

	# SVector{9}([p_floor; p_object; p_pusher])
	# p_floor = kinematics_1(model, q; body = :floor)
	# SVector{3}([p_floor; ])
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
	# SVector{3}([q[3] - model.env.surf(x[1:2]),
	# 			norm(Δ) - (model.r + model.rp),
	# 			norm(Δ) - (model.r + model.rp),
	# 			])
	# SVector{3}([q[3] - model.env.surf(x[1:2]),
	# 			x[1],
	# 			xp[1],
	# 			])
	# SVector{2}([q[3] - model.env.surf(x[1:2]),
	# 			x[1],
	# 			])

	SVector{2}([q[3] - model.env.surf(x[1:2]),
				norm(Δ) - (model.r + model.rp),
				# norm(Δ) - (model.r + model.rp),
				])

	# SVector{3}([q[3] - model.env.surf(x[1:2]),
	# 			norm(Δ) - (model.r + model.rp),
	# 			norm(Δ) - (model.r + model.rp),
	# 			])


	# SVector{1}([q[3] - model.env.surf(x[1:2]),
	# 			# norm(Δ) - (model.r + model.rp),
	# 			# norm(Δ) - (model.r + model.rp),
	# 			])
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
			   # 0.0  0.0  0.0  0.0  0.0  0.0  0.0]
			   0.0  0.0  1.0  0.0  0.0  0.0  0.0]
	# p_object = [x + r cos(α) , y + sin(α)]
	# v_object = [ẋ - r sin(α) θd, ẏ + r cos(α) θd]
	r = model.r
	α = contact_angle(model, q)
    J_object = [1.0  0.0  0.0 -r*sin(α) 0.0  0.0  0.0;
			    0.0  1.0  0.0  r*cos(α) 0.0  0.0  0.0;
				# 0.0  0.0  0.0  0.0      0.0  0.0  0.0]
			    0.0  0.0  1.0  0.0  0.0  0.0  0.0]

	J_object = [1.0  0.0  0.0 -r*sin(α)  -1.0  0.0  0.0;
			    0.0  1.0  0.0  r*cos(α)  0.0  -1.0  0.0;
				# 0.0  0.0  0.0  0.0      0.0  0.0  0.0]
			    0.0  0.0  1.0  0.0       0.0  0.0  -1.0] # relative speed in world frame
				# of the two contact points when mutiplied by qdot
				# one contact point attached to a object and the other one attached to pusher.


	# J_pusher = [0.0  0.0  0.0  0.0  1.0  0.0  0.0;
	# 		    0.0  0.0  0.0  0.0  0.0  1.0  0.0;
	# 			# 0.0  0.0  0.0  0.0  0.0  0.0  0.0]
	# 		    0.0  0.0  0.0  0.0  0.0  0.0  1.0]
	# return [J_floor;
	# 		J_object;
	# 		J_pusher] # (nc*np) x nq  = 9x7

	return [J_floor;
			J_object;
			] # (nc*np) x nq  = 9x7

	# return [J_floor;
	# 		J_object;
	# 		J_pusher] # (nc*np) x nq  = 9x7
	# return [J_floor;
	# 		] # (nc*np) x nq  = 9x7
end

function contact_forces(model::PlanarPush13, γ1, b1, q2, k)
	# Express the contact forces in the world frame.
	# returns λ1
	# which will be used in dynamics: transpose(J_fast(model, q2)) * λ1
	m = friction_mapping(model.env)
	# get rotation
	# α = contact_angle(model, q2)
	# cα = cos(α)
	# sα = sin(α)
	# R = [-sα  0.0  cα;
	# 	  cα  0.0  sα;
	# 	 0.0  1.0  0.0]# W -> S
	# # Express in world frame
	# λ_floor = [m * b1[1:4]; γ1[1]]
	# λ_object_c = [m * b1[5:8]; γ1[2]]
	# λ_pusher_c = [m * b1[9:12]; γ1[3]]
	# λ_object =  R * λ_object_c # rotate in world frame
	# # λ_pusher = -R * λ_pusher_c # action reaction principle
	# λ_pusher =  R * λ_pusher_c # action reaction principle
	# SVector{9}([λ_floor; λ_object; λ_pusher])

	# λ_floor = [m * b1[1:4]; γ1[1]]
	# λ_object_c = [γ1[2]; m * b1[5:8]]
	# λ_pusher_c = [γ1[3]; m * b1[9:12]]
	# λ_object =  λ_object_c # rotate in world frame
	# # λ_pusher = -R * λ_pusher_c # action reaction principle
	# λ_pusher =  λ_pusher_c # action reaction principle
	# SVector{9}([λ_floor; λ_object; λ_pusher])


	λ_floor = [m * b1[1:4]; γ1[1]]

	# β = α-π # angle between pusher and object
	# b_object_c = m * b1[5:8]
	# λ_object_c = [-sin(β)*b_object_c[1], cos(β)*b_object_c[1], b_object_c[2]]
	# λ_object_c += [cos(β) * γ1[2], sin(β) * γ1[2], 0.0]
	# λ_object =  λ_object_c # rotate in world frame
	λ_object_c = [m * b1[5:8]; γ1[2]]

	α = contact_angle(model, q2) + π
	cα = cos(α)
	sα = sin(α)
	R = [-sα  0.0  cα;
		  cα  0.0  sα;
		 0.0  1.0  0.0] # S -> W

	λ_object =  R*λ_object_c # rotate in world frame

	# λ_pusher = -R * λ_pusher_c # action reaction principle
	# λ_pusher =  λ_pusher_c # action reaction principle

	SVector{6}([λ_floor; λ_object])

	# γ = α # angle between object and pusher
	# # γ = α-π # angle between object and pusher
	# b_pusher_c = m * b1[9:12]
	# λ_pusher_c = [-sin(γ)*b_pusher_c[1], cos(γ)*b_pusher_c[1], b_pusher_c[2]]
	# λ_pusher_c += [cos(γ) * γ1[3], sin(γ) * γ1[3], 0.0]
	# # λ_pusher_c += [-cos(γ) * γ1[3], -sin(γ) * γ1[3], 0.0]
	# # λ_pusher =  [λ_pusher_c[1], λ_pusher_c[2], λ_pusher_c[3]] # rotate in world frame
	# λ_pusher =  λ_pusher_c # rotate in world frame
	# SVector{9}([λ_floor; λ_object; λ_pusher])
	# # SVector{3}([λ_floor;])
end

function velocity_stack(model::PlanarPush13, q1, q2, k, h)
	# k = kinematics(model, q2)
	v = J_func(model, q2) * (q2 - q1) / h[1]

	α = contact_angle(model, q2) + π
	cα = cos(α)
	sα = sin(α)
	R = [-sα  0.0  cα;
		  cα  0.0  sα;
		 0.0  1.0  0.0] # S -> W
	# In the world frame
	v_floor = rotation(model.env, k) * v[1:3]
	v_object_w = rotation(model.env, k) * v[4:6]
	v_object = R' * v_object_w # W -> S
	# v_object = [cα*v[4] + sα*v[5], -sα*v[4] + cα*v[5], v[6]]

	# v_pusher = rotation(model.env, k) * v[7:9]

	# SVector{12}([v_floor[1:2];  -v_floor[1:2];
	# 			v_object[1:2]; -v_object[1:2];
	# 			v_pusher[1:2]; -v_pusher[1:2];])

	# SVector{8}([v_floor[1:2];  -v_floor[1:2];
	# 			v_object[1:2]; -v_object[1:2];])

	SVector{8}([v_floor[1:2];  -v_floor[1:2];
				# v_object[2:3]; -v_object[2:3];
				# v_object[1:2]; -v_object[1:2];
				friction_mapping(model.env)' * v_object[1:2];
				# v_object[[3,2]]; -v_object[[3,2]];
				])

	# SVector{12}([v_floor[1:2];  -v_floor[1:2];
	# 			v_object[[3,2]]; -v_object[[3,2]];
	# 			# v_pusher[[2,3]]; -v_pusher[[2,3]];
	# 			# -v_pusher[[3,2]];  v_pusher[[3,2]];
	# 			v_pusher[[3,2]]; -v_pusher[[3,2]];
	# 			])

	# SVector{4}([v_floor[1:2];  -v_floor[1:2];
	# 			# v_object[1:2]; -v_object[1:2];
	# 			# v_pusher[1:2]; -v_pusher[1:2];
	# 			])
end
#
# v1_surf = rotation(model.env, k) * v
#
# SVector{4}(friction_mapping(model.env)' * v1_surf[1:2])
#

# Dimensions
nq = 3 + 1 + 3            # configuration dimension
nu = 2                    # control dimension
# nc = 3                    # number of contact points
nc = 2                    # number of contact points
# nc = 1                    # number of contact points
nf = 4                    # number of parameters for friction cone
nb = nc * nf              # number of friction parameters
nw = 3                    # disturbance dimension

# World parameters
μ_world = 0.9      # coefficient of friction
# μ_world = 0.0      # coefficient of friction
μ_joint = 0.0
g = 9.81     # gravity
# g = 0.00     # gravity

# Model parameters
m = 0.2
J = 0.0025
mp = 2.0
r = 0.1
rp = 0.01


planarpush = PlanarPush13(Dimensions(nq, nu, nw, nc, nb),
			  μ_world, μ_joint, g,
			  m, J, mp, r, rp,
			  zeros(nc),
			  BaseMethods(), DynamicsMethods(), ContactMethods(),
			  ResidualMethods(), ResidualMethods(),
			  SparseStructure(spzeros(0, 0), spzeros(0, 0)),
			  SVector{nq}([zeros(3); 0.0 * μ_joint * ones(nq - 3)]),
			  environment_3D_flat())
#
# @variables q[1:7]
# @variables q̇[1:7]
#
# L = lagrangian(planarpush, q, q̇)
# ddLdq̇q̇ = Symbolics.sparsehessian(L, q̇, simplify=true)
