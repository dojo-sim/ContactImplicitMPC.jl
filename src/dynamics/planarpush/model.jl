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

mutable struct PlanarPush{T} <: ContactModel
    dim::Dimensions

	μ_world::T  # coefficient of friction
	μ_joint::T  # gravity
	g::T

    m::T  # mass of object
	J::T  # inertia of object
    mp::T # mass of pusher

	r::T  # radius of object
	rp::T # radius of pusher

	base::BaseMethods
	dyn::DynamicsMethods

	joint_friction::SVector
end

# Lagrangian
function lagrangian(model::PlanarPush, q, q̇)
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

function _dLdq(model::PlanarPush, q, q̇)
	Lq(x) = lagrangian(model, x, q̇)
	ForwardDiff.gradient(Lq, q)
end

function _dLdq̇(model::PlanarPush, q, q̇)
	Lq̇(x) = lagrangian(model, q, x)
	ForwardDiff.gradient(Lq̇, q̇)
end

function _C_func(model::PlanarPush, q, q̇)
	tmp_q(z) = _dLdq̇(model, z, q̇)
	tmp_q̇(z) = _dLdq̇(model, q, z)

	ForwardDiff.jacobian(tmp_q, q) * q̇ - _dLdq(model, q, q̇)
end

function kinematics_1(model::PlanarPush, q; body = :floor)
	# Contact 1 is between object and floor: :floor
	# Contact 2 is between object and pusher: :object, :pusher
	x = [q[1], q[2], 0.0]
	xp = [q[5], q[6], 0.0]
	Δ = x - xp
	Δ ./= norm(Δ)
	rc = model.r/sqrt(2)
	θ = q[4]
	if body == :floor
		return x
	elseif body == :floor1
		# we use sqrt(2) to simulate a uniform pressure field under the disk,
		# it will result in the same frictional torque
		return x + [cos(θ)*rc, sin(θ)*rc, 0.0]
	elseif body == :floor2
		return x - [cos(θ)*rc, sin(θ)*rc, 0.0]
	elseif body == :object
		return x - Δ*model.r
	elseif body == :pusher
		return xp + Δ*model.rp
	else
		@error "incorrect body specification"
	end
end

function kinematics(model::PlanarPush, q)
	p_floor = kinematics_1(model, q; body = :floor)
	p_floor1 = kinematics_1(model, q; body = :floor1)
	p_floor2 = kinematics_1(model, q; body = :floor2)
	p_object = kinematics_1(model, q; body = :object)
	p_pusher = kinematics_1(model, q; body = :pusher)
	SVector{9}([p_floor1; p_floor2; p_object])
end

# Methods
function M_func(model::PlanarPush, q)
	Diagonal(@SVector [model.m, model.m, model.m,
					   model.J,
					   model.mp, model.mp, model.mp,])
end

function ϕ_func(model::PlanarPush, env::Environment, q)
	x = kinematics_1(model, q, body=:floor)
	x1 = kinematics_1(model, q, body=:floor1)
	x2 = kinematics_1(model, q, body=:floor2)
	xp = [q[5], q[6], 0.0]
	Δ = x - xp
	SVector{3}([q[3] - env.surf(x1[1:2]),
				q[3] - env.surf(x2[1:2]),
				norm(Δ) - (model.r + model.rp),
				])
end

function B_func(model::PlanarPush, q)
	@SMatrix [0.0  0.0  0.0  0.0  1.0  0.0  0.0;
			  0.0  0.0  0.0  0.0  0.0  1.0  0.0]
end

function A_func(model::PlanarPush, q)
	@SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 1.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 1.0 0.0 0.0 0.0 0.0]
end

function contact_angle(model::PlanarPush, q)
	# Angle of the contact normal
	x = q[1]
	y = q[2]
	xp = q[5]
	yp = q[6]
	α = atan(yp - y, xp - x)
	return α
end

function rotation_s_to_w(model::PlanarPush, q2)
	# rotation matrix from frame attached to the contact point on the object
	# (b1, and b2 are orthogonal to the object, b2 points in the z direction, γ points inside the object)
	# to the world frame (x,y,z)
	β = contact_angle(model, q2) + π
	cβ = cos(β)
	sβ = sin(β)
	R = [-sβ  0.0  cβ;
		  cβ  0.0  sβ;
		 0.0  1.0  0.0] # S -> W
	return R
end

function J_func(model::PlanarPush, env::Environment{<:World,LinearizedCone}, q)
	r = model.r
	rc = model.r/sqrt(2)
	α = contact_angle(model, q)
	θ = q[4]
	# Jacobian velocity in world frame of contact points ATTACHED to their respective links.
	J_floor1 = [1.0  0.0  0.0 -rc*sin(θ)  0.0  0.0  0.0;
			    0.0  1.0  0.0  rc*cos(θ)  0.0  0.0  0.0;
			    0.0  0.0  1.0  0.0        0.0  0.0  0.0]
	J_floor2 = [1.0  0.0  0.0  rc*sin(θ)  0.0  0.0  0.0;
				0.0  1.0  0.0 -rc*cos(θ)  0.0  0.0  0.0;
				0.0  0.0  1.0  0.0        0.0  0.0  0.0]
	# p_object = [x + r cos(α) , y + sin(α)]
	# v_object = [ẋ - r sin(α) θd, ẏ + r cos(α) θd]
	J_object = [1.0  0.0  0.0 -r*sin(α)  -1.0  0.0  0.0;
			    0.0  1.0  0.0  r*cos(α)   0.0 -1.0  0.0;
			    0.0  0.0  1.0  0.0        0.0  0.0 -1.0] # relative speed in world frame
				# of the two contact points when mutiplied by qdot
				# one contact point attached to a object and the other one attached to pusher.
	return [J_floor1; J_floor2; J_object;] # (nc*np) x nq  = 6x7
end

function contact_forces(model::PlanarPush, env::Environment{<:World,LinearizedCone}, γ1, b1, q2, k)
	# Express the contact forces in the world frame.
	# returns λ1
	# which will be used in dynamics: transpose(J_fast(model, q2)) * λ1
	m = friction_mapping(env)
	# In the surface frame
	λ_floor1 = [m * b1[1:4]; γ1[1]]
	λ_floor2 = [m * b1[5:8]; γ1[2]]
	λ_object_c = [m * b1[9:12]; γ1[3]]

	# In the world frame
	R = rotation_s_to_w(model, q2)
	λ_object =  R*λ_object_c # rotate in world frame
	SVector{9}([λ_floor1; λ_floor2; λ_object])
end

function velocity_stack(model::PlanarPush, env::Environment{<:World,LinearizedCone}, q1, q2, k, h)
	# In the world frame
	v = J_func(model, q2) * (q2 - q1) / h[1]

	R = rotation_s_to_w(model, q2)
	v_floor1 = rotation(env, k) * v[1:3]
	v_floor2 = rotation(env, k) * v[4:6]
	v_object_w = v[7:9]
	# We express in the surface frame
	v_object = R' * v_object_w # W -> S
	SVector{12}([v_floor1[1:2];  -v_floor1[1:2];
				 v_floor2[1:2];  -v_floor2[1:2];
				friction_mapping(env)' * v_object[1:2];
				])
end



# Dimensions
nq = 3 + 1 + 3            # configuration dimension
nu = 2                    # control dimension
nc = 3                    # number of contact points
nf = 4                    # number of parameters for friction cone
nb = nc * nf              # number of friction parameters
nw = 3                    # disturbance dimension

# World parameters
μ_world = 0.50      # coefficient of friction
μ_joint = 0.0
g = 9.81     # gravity

# Model parameters
m = 10.0
J = m * 0.05^2
mp = 100.0
r = 0.2
rp = 0.04


planarpush = PlanarPush(Dimensions(nq, nu, nw, nc),
			  μ_world, μ_joint, g,
			  m, J, mp, r, rp,
			  BaseMethods(), DynamicsMethods(),
			  SVector{nq}([zeros(3); 0.0 * μ_joint * ones(nq - 3)]))
