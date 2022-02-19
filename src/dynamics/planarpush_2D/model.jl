"""
    PlanarPush2D
        s = (x, y, xp, yp)
            x  - x position of the object
            y  - y position of the object
            xp - x position of the pusher
            yp - y position of the pusher
"""
mutable struct PlanarPush2D{T} <: ContactModel
    dim::Dimensions

	μ_world::T  # coefficient of friction
	μ_joint::T  # gravity
	g::T

    m::T  # mass of object
    mp::T # mass of pusher

	r::T  # radius of object
	rp::T # radius of pusher

	base::BaseMethods
	dyn::DynamicsMethods

	joint_friction::SVector
end

# Lagrangian
function lagrangian(model::PlanarPush2D, q, q̇)
	v_object = [q̇[1], q̇[2]]
	v_pusher = [q̇[3], q̇[4]]

	L = 0.0
	# kinetic energy
	L += 0.5 * model.m * transpose(v_object) * v_object
	L += 0.5 * model.mp * transpose(v_pusher) * v_pusher
	# potential energy
	L -= model.m * model.g * q[2]
	# L -= model.mp * model.g * q[4] # the pusher floats in the air, not subject to gravity
	return L
end

function _dLdq(model::PlanarPush2D, q, q̇)
	Lq(x) = lagrangian(model, x, q̇)
	ForwardDiff.gradient(Lq, q)
end

function _dLdq̇(model::PlanarPush2D, q, q̇)
	Lq̇(x) = lagrangian(model, q, x)
	ForwardDiff.gradient(Lq̇, q̇)
end

function _C_func(model::PlanarPush2D, q, q̇)
	tmp_q(z) = _dLdq̇(model, z, q̇)
	tmp_q̇(z) = _dLdq̇(model, q, z)

	ForwardDiff.jacobian(tmp_q, q) * q̇ - _dLdq(model, q, q̇)
end

function kinematics_1(model::PlanarPush2D, q; body = :floor)
	# Contact 1 is between object and floor: :floor
	# Contact 2 is between object and pusher: :object, :pusher
	x = [q[1], q[2]]
	xp = [q[3], q[4]]
	if body == :floor
		return x
	elseif body == :object
		return x - [model.r, 0.0]
	elseif body == :pusher
		return xp + [model.rp, 0.0]
	else
		@error "incorrect body specification"
	end
end

function kinematics(model::PlanarPush2D, q)
	p_floor = kinematics_1(model, q; body = :floor)
	p_object = kinematics_1(model, q; body = :object)
	p_pusher = kinematics_1(model, q; body = :pusher)
	SVector{4}([p_floor; p_object])
end

# Methods
function M_func(model::PlanarPush2D, q)
	Diagonal(@SVector [model.m, model.m, model.mp, model.mp,])
end

function ϕ_func(model::PlanarPush2D, env::Environment, q)
	x = kinematics_1(model, q, body=:floor)
	Δ = q[1] - q[3]
	SVector{2}([q[2] - env.surf(x[1]),
				Δ - (model.r + model.rp),
				])
end

function B_func(model::PlanarPush2D, q)
	@SMatrix [0.0  0.0  1.0  0.0]
end

function A_func(model::PlanarPush2D, q)
	@SMatrix [1.0  0.0  0.0  0.0;
			  0.0  1.0  0.0  0.0]
end

# function contact_angle(model::PlanarPush2D, q)
# 	# Angle of the contact normal
# 	x = q[1]
# 	y = q[2]
# 	xp = q[5]
# 	yp = q[6]
# 	α = atan(yp - y, xp - x)
# 	return α
# end
#
# function rotation_s_to_w(model::PlanarPush2D, q2)
# 	# rotation matrix from frame attached to the contact point on the object
# 	# (b1, and b2 are orthogonal to the object, b2 points in the z direction, γ points inside the object)
# 	# to the world frame (x,y,z)
# 	β = contact_angle(model, q2) + π
# 	cβ = cos(β)
# 	sβ = sin(β)
# 	R = [-sβ  0.0  cβ;
# 		  cβ  0.0  sβ;
# 		 0.0  1.0  0.0] # S -> W
# 	return R
# end

function J_func(model::PlanarPush2D, env::Environment{<:World,LinearizedCone}, q)
	# Jacobian velocity in world frame of contact points ATTACHED to their respective links.
	J_floor  = [1.0  0.0  0.0  0.0;
			    0.0  1.0  0.0  0.0]
	# p_object = [x + r cos(α) , y + sin(α)]
	# v_object = [ẋ - r sin(α) θd, ẏ + r cos(α) θd]
	J_object = [1.0  0.0 -1.0  0.0;
			    0.0  1.0  0.0 -1.0] # relative speed in world frame
				# of the two contact points when mutiplied by qdot
				# one contact point attached to a object and the other one attached to pusher.
	return [J_floor; J_object;] # (nc*np) x nq  = 4x4
end

function contact_forces(model::PlanarPush2D, env::Environment{<:World,LinearizedCone}, γ1, b1, q2, k)
	# Express the contact forces in the world frame.
	# returns λ1
	# which will be used in dynamics: transpose(J_fast(model, q2)) * λ1
	m = friction_mapping(env)
	# In the surface frame
	λ_floor_s = [m * b1[1:2]; γ1[1]] # friction along xs, impact along ys
	λ_object_s = [m * b1[3:4]; γ1[2]] # friction along xs, impact along ys

	# In the world frame
	λ_floor_w = [λ_floor_s[1], λ_floor_s[2]] # friction along xw, impact along yw
	λ_object_w = [λ_object_s[2], -λ_object_s[1]] # friction along -yw, impact along xw

	SVector{4}([λ_floor_w; λ_object_w])
end

function velocity_stack(model::PlanarPush2D, env::Environment{<:World,LinearizedCone}, q1, q2, k, h)
	# In the world frame
	v = J_func(model, env, q2) * (q2 - q1) / h[1]

	v_floor_w = v[1:2]
	v_object_w = v[3:4]
	# We express in the surface frame
	v_floor_s = rotation(env, k[1:2]) * v_floor_w
	R = RotZ(- π / 2) # Rs->w for object contact
	v_object_s = R'[1:2,1:2] * v_object_w # W -> S
	SVector{4}([
				friction_mapping(env)' * v_floor_s[1];
				friction_mapping(env)' * v_object_s[1];
				])
end


# Dimensions
nq = 2 + 2                # configuration dimension
nu = 1                    # control dimension
nc = 2                    # number of contact points
nw = 2                    # disturbance dimension
nquat = 0

# World parameters
# μ_world = 0.50      # coefficient of friction
μ_world = 0.10      # coefficient of friction
μ_joint = 0.0
g = 9.81     # gravity

# Model parameters
m = 10.0
mp = 100.0
r = 0.2
rp = 0.04


planarpush_2D = PlanarPush2D(Dimensions(nq, nu, nw, nc, nquat),
			    μ_world, μ_joint, g,
			    m, mp, r, rp,
			    BaseMethods(), DynamicsMethods(),
			    SVector{nq}([zeros(2); 0.0 * μ_joint * ones(nq - 2)]))
