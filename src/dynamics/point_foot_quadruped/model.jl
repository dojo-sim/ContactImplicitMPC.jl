# """
#     point foot quadruped with linear springs and dampers on each foot
#     q = (p, r, f1, f2, f3, f4)
#         p - body position
#         r - body orientation (modified Rodriques parameters)
#         f1 - foot 1 position
#         f2 - foot 2 position
#         f3 - foot 3 position
#         f4 - foot 4 position
# """
# mutable struct PointFootQuadruped117{T} <: Model{T}
# 	# dimensions
# 	nq::Int # generalized coordinates
# 	nu::Int # controls
# 	nw::Int # parameters
# 	nc::Int # contact points
#
# 	# parameters
# 	body_height::T
# 	foot_x::T
# 	foot_y::T
#
# 	# parameters
# 	mass_body::T
# 	inertia_body::Matrix{T}
# 	mass_foot::T
#
# 	# fast methods
# 	base
# 	dyn
#
# 	# parameters
# 	joint_friction::Vector{T}
# 	damper_foot::Vector{T}
# 	spring_stiffness_joint::Vector{T}
#
# 	# environment
#     friction_foot_world::Vector{T}
# 	g::T
# end
#
# function L_mult(x)
#     [x[1] -transpose(x[2:4]);
#      x[2:4] x[1] * I(3) + skew(x[2:4])]
# end
#
# # right quaternion multiply as matrix
# function R_mult(x)
#     [x[1] -transpose(x[2:4]); x[2:4] x[1] * I(3) - skew(x[2:4])]
# end
#
# # rotation matrix
# function quat_rot_mat(q)
#     H = [zeros(1, 3); I(3)]
#     transpose(H) * L_mult(q) * transpose(R_mult(q)) * H
# end
#
# function quat_from_mrp(p)
#     """Quaternion (scalar first) from MRP"""
#     return (1.0 / (1.0 + dot(p, p))) * [(1 - dot(p, p)); 2.0 * p]
# end
#
# function mrp_rot_mat(x)
#     quat_rot_mat(quat_from_mrp(x))
# end
#
# # Kinematics
# function kinematics(model::PointFootQuadruped117, q)
# 	[
# 	q[6 .+ (1:3)];
# 	q[6 .+ (4:6)];
# 	q[6 .+ (7:9)];
# 	q[6 .+ (10:12)]
# 	]
# end
#
# lagrangian(model::PointFootQuadruped117, q, q̇) = 0.0
#
# function M_func(model::PointFootQuadruped117, q)
#     cat(
#         model.mass_body * Diagonal(ones(3)),     # body position
#         model.inertia_body,                      # body orienation
#         model.mass_foot * Diagonal(ones(3 * 4)), # feet position
#         dims=(1, 2)
#         )
# end
#
# function C_func(model::PointFootQuadruped117, q, q̇)
# 	offsets = [
# 		[+model.foot_x, +model.foot_y, -model.body_height],
# 		[+model.foot_x, -model.foot_y, -model.body_height],
# 		[-model.foot_x, +model.foot_y, -model.body_height],
# 		[-model.foot_x, -model.foot_y, -model.body_height],
# 		]
#
# 	spring = model.spring_stiffness_joint
#     [
#         model.mass_body * [0,0,model.g] + spring .* (+sum(offsets) + 4 * q[1:3] - q[7:9] - q[10:12] - q[13:15] - q[16:18]);     # body position
#         skew(q̇[4:6]) * model.inertia_body * q̇[4:6]; # body orientation
#         model.mass_foot * [0,0,model.g];# + spring .* (-offsets[1] + q[7:9] - q[1:3]);
#         model.mass_foot * [0,0,model.g];# + spring .* (-offsets[2] + q[10:12] - q[1:3]);
#         model.mass_foot * [0,0,model.g];# + spring .* (-offsets[3] + q[13:15] - q[1:3]);
#         model.mass_foot * [0,0,model.g];# + spring .* (-offsets[4] + q[16:18] - q[1:3]);
#     ]
# end
#
# function ϕ_func(model::PointFootQuadruped117, env::Environment, q)
#
#     position_foot1 = q[6 .+ (1:3)]
#     position_foot2 = q[9 .+ (1:3)]
#     position_foot3 = q[12 .+ (1:3)]
# 	position_foot4 = q[15 .+ (1:3)]
#
# 	return [position_foot1[3]; position_foot2[3]; position_foot3[3]; position_foot4[3]]
# end
#
# function B_func(model::PointFootQuadruped117, q)
#     position_body = q[1:3]
#     orientation_body = q[3 .+ (1:3)]
# 	R = mrp_rot_mat(orientation_body)
# 	# R = euler_rotation_matrix(orientation_body)
#
# 	# kinematics in world frame
# 	r1 = q[6 .+ (1:3)] - position_body
# 	r2 = q[9 .+ (1:3)] - position_body
# 	r3 = q[12 .+ (1:3)] - position_body
# 	r4 = q[15 .+ (1:3)] - position_body
#
# 	z3 = zeros(3, 3)
#
# 	transpose([
#         I(3) z3   I(3) I(3) I(3) I(3);
# 		z3   I(3) transpose(R) * skew(r1) transpose(R) * skew(r2) transpose(R) * skew(r3) transpose(R) * skew(r4);
#         z3   z3   -I(3)    z3    z3   z3;
#         z3   z3   z3    -I(3)    z3   z3;
#         z3   z3   z3       z3 -I(3)   z3;
#         z3   z3   z3       z3    z3 -I(3)
#     ])
# end
#
# function A_func(model::PointFootQuadruped117, q)
#     @SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
#               0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
# 			  0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
# end
#
# function J_func(model::PointFootQuadruped117, env::Environment, q)
#     z3 = zeros(3, 3)
#
#     [
#         z3   z3 I(3)   z3   z3   z3;
# 		z3   z3   z3 I(3)   z3   z3;
#         z3   z3   z3   z3 I(3)   z3;
#         z3   z3   z3   z3   z3 I(3);
#     ]
# end
#
# function nominal_configuration(model::PointFootQuadruped117)
#     [
#         0.0; 0.0; model.body_height;
#         0.0; 0.0; 0.0;
#         model.foot_x; model.foot_y; 0.0;
#         model.foot_x;-model.foot_y; 0.0;
#        -model.foot_x; model.foot_y; 0.0;
#        -model.foot_x;-model.foot_y; 0.0;
#     ]
# end
#
# function nominal_state(model::PointFootQuadruped117)
# 	[nominal_configuration(model); zeros(model.nq)]
# end
#
# # friction coefficients
# function friction_coefficients(model::PointFootQuadruped117)
#     model.friction_foot_world
# end
#
# # dimensions
# nq = 3 + 3 + 3 * 4       # generalized coordinates
# nu = 6 + 3 * 4           # controls
# nw = 3                   # parameters
# nc = 4                   # contact points
#
# # inertial properties
# mass_body = 13.5
# i_xx = 0.0178533
# i_xy = 0.0
# i_xz = 0.0
# i_yz = 0.0
# i_yy = 0.0377999
# i_zz = 0.0456542
# inertia_body = Array(Diagonal([i_xx, i_yy, i_zz]))
# mass_foot = 0.2
#
# # parameters
# body_height = 0.3
# foot_x = 0.17
# foot_y = 0.15
#
# joint_friction = 0.0 * ones(nq)
# damper_foot = 0.15 * ones(3)
# spring_stiffness_joint = [5,5,15.0]
# friction_foot_world = 0.3*ones(1)     # coefficient of friction
# gravity = 9.81                   # gravity
#
#
# point_foot_quadruped = PointFootQuadruped117(nq, nu, nw, nc,
# 				body_height,
# 				foot_x,
# 				foot_y,
# 				mass_body,
# 				inertia_body,
# 				mass_foot,
# 				BaseMethods(),
# 				DynamicsMethods(),
# 				joint_friction,
# 				damper_foot,
# 				spring_stiffness_joint,
# 				friction_foot_world,
# 				gravity,
#                 )
#
# # point_foot_quadruped_contact_kinematics = [
# #     q -> q[6  .+ (1:3)],
# #     q -> q[9  .+ (1:3)],
# #     q -> q[12 .+ (1:3)],
# #     q -> q[15 .+ (1:3)],
# # ]
# #
# # point_foot_quadruped_contact_kinematics_jacobians = [
# #     q -> [zeros(3, 6) I(3) zeros(3, 9)],
# #     q -> [zeros(3, 9) I(3) zeros(3, 6)],
# #     q -> [zeros(3, 12) I(3) zeros(3, 3)],
# #     q -> [zeros(3, 15) I(3)],
# # ]
#
# floating_base_dim(::PointFootQuadruped117) = 6
# name(::PointFootQuadruped117) = :point_foot_quadruped
#
# function dynamics(model::PointFootQuadruped117, h, q0, q1, u1, w1, Λ1, q2)
# 	# evalutate at midpoint
# 	qm1 = 0.5 * (q0 + q1)
#     vm1 = (q1 - q0) / h[1]
#     qm2 = 0.5 * (q1 + q2)
#     vm2 = (q2 - q1) / h[1]
#
# 	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
# 	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)
#
# 	# return (0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
# 	# 	+ transpose(B_fast(model, qm2)) * u1
# 	# 	+ transpose(A_fast(model, qm2)) * w1
# 	# 	# + transpose(J_fast(model, q2)) * λ1
# 	# 	+ Λ1
# 	# 	- h[1] * model.joint_friction .* vm2)
# 	d = 0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2 # variational integrator (midpoint)
# 	d .+= transpose(B_fast(model, qm2)) * u1        # control inputs
# 	d .+= transpose(A_fast(model, qm2)) * w1        # control inputs
# 	d .+= Λ1                                        # contact impulses
#
# 	# d[1:3] .-= model.damper_foot .* (vm2[1:3] - vm2[7:9]) # joint friction
# 	# d[1:3] .-= model.damper_foot .* (vm2[1:3] - vm2[10:12]) # joint friction
# 	# d[1:3] .-= model.damper_foot .* (vm2[1:3] - vm2[13:15]) # joint friction
# 	# d[1:3] .-= model.damper_foot .* (vm2[1:3] - vm2[16:18]) # joint friction
# 	#
# 	# d[7:9]   .-= model.damper_foot .* (vm2[7:9] - vm2[1:3]) # joint friction
# 	# d[10:12] .-= model.damper_foot .* (vm2[10:12] - vm2[1:3]) # joint friction
# 	# d[13:15] .-= model.damper_foot .* (vm2[13:15] - vm2[1:3]) # joint friction
# 	# d[16:18] .-= model.damper_foot .* (vm2[16:18] - vm2[1:3]) # joint friction
# 	return d
# end

"""
    point foot quadruped
    q = (p, r, f1, f2, f3, f4)
        p - body position
        r - body orientation (modified Rodriques parameters)
        f1 - foot 1 position
        f2 - foot 2 position
        f3 - foot 3 position
        f4 - foot 4 position
"""
mutable struct PointFootQuadruped117{T} <: Model{T}
    # dimensions
	nq::Int # generalized coordinates
    nu::Int # controls
    nw::Int # parameters
    nc::Int # contact points

	# parameters
	body_height::T
	foot_x::T
	foot_y::T

	# parameters
    mass_body::T
    inertia_body::Matrix{T}
    mass_foot::T

	# environment
	friction_foot_world::T
	joint_friction::Vector{T}
	spring_stiffness_joint::Vector{T}
	g::T

	# fast methods
	base
	dyn
end

function L_mult(x)
    [x[1] -transpose(x[2:4]);
     x[2:4] x[1] * I(3) + skew(x[2:4])]
end

# right quaternion multiply as matrix
function R_mult(x)
    [x[1] -transpose(x[2:4]); x[2:4] x[1] * I(3) - skew(x[2:4])]
end

# rotation matrix
function quat_rot_mat(q)
    H = [zeros(1, 3); I(3)]
    transpose(H) * L_mult(q) * transpose(R_mult(q)) * H
end

function quat_from_mrp(p)
    """Quaternion (scalar first) from MRP"""
    return (1.0 / (1.0 + dot(p, p))) * [(1 - dot(p, p)); 2.0 * p]
end

function mrp_rot_mat(x)
    quat_rot_mat(quat_from_mrp(x))
end

# Kinematics
function kinematics(model::PointFootQuadruped117, q)
	q[6 .+ (1:12)]
end

lagrangian(model::PointFootQuadruped117, q, q̇) = 0.0

function M_func(model::PointFootQuadruped117, q)
    cat(
        model.mass_body * Diagonal(ones(3)),     # body position
        model.inertia_body,                      # body orienation
        model.mass_foot * Diagonal(ones(3 * 4)), # feet position
        dims=(1, 2)
        )
end

function C_func(model::PointFootQuadruped117, q, q̇)
    # [
    #     model.mass_body * [0,0,model.g];            # body position
    #     skew(q̇[4:6]) * model.inertia_body * q̇[4:6]; # body orienation
    #     model.mass_foot * [0,0,model.g];
    #     model.mass_foot * [0,0,model.g];
    #     model.mass_foot * [0,0,model.g];
    #     model.mass_foot * [0,0,model.g];
    # ]

	offsets = [
		[+model.foot_x, +model.foot_y, -model.body_height],
		[+model.foot_x, -model.foot_y, -model.body_height],
		[-model.foot_x, +model.foot_y, -model.body_height],
		[-model.foot_x, -model.foot_y, -model.body_height],
		]

	spring = model.spring_stiffness_joint
	[
	    model.mass_body * [0,0,model.g] + spring .* (+sum(offsets) + 4 * q[1:3] - q[7:9] - q[10:12] - q[13:15] - q[16:18]);     # body position
	    skew(q̇[4:6]) * model.inertia_body * q̇[4:6]; # body orientation
	    model.mass_foot * [0,0,model.g];# + spring .* (-offsets[1] + q[7:9] - q[1:3]);
	    model.mass_foot * [0,0,model.g];# + spring .* (-offsets[2] + q[10:12] - q[1:3]);
	    model.mass_foot * [0,0,model.g];# + spring .* (-offsets[3] + q[13:15] - q[1:3]);
	    model.mass_foot * [0,0,model.g];# + spring .* (-offsets[4] + q[16:18] - q[1:3]);
	]
end

function ϕ_func(model::PointFootQuadruped117, env::Environment, q)

    position_foot1 = q[6 .+ (1:3)]
    position_foot2 = q[9 .+ (1:3)]
    position_foot3 = q[12 .+ (1:3)]
	position_foot4 = q[15 .+ (1:3)]

	return [position_foot1[3]; position_foot2[3]; position_foot3[3]; position_foot4[3]]
end

function B_func(model::PointFootQuadruped117, q)
    position_body = q[1:3]
    orientation_body = q[3 .+ (1:3)]
	# R = mrp_rot_mat(orientation_body)
	R = euler_rotation_matrix(orientation_body)

	# kinematics in world frame
	r1 = q[6 .+ (1:3)] - position_body
	r2 = q[9 .+ (1:3)] - position_body
	r3 = q[12 .+ (1:3)] - position_body
	r4 = q[15 .+ (1:3)] - position_body

	z3 = zeros(3, 3)

	transpose([
        I(3) I(3) I(3) I(3);
		transpose(R) * skew(r1) transpose(R) * skew(r2) transpose(R) * skew(r3) transpose(R) * skew(r4);
        -I(3)    z3    z3   z3;
        z3    -I(3)    z3   z3;
        z3       z3 -I(3)   z3;
        z3       z3    z3 -I(3)
    ])
end

function A_func(model::PointFootQuadruped117, q)
    @SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
              0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
end

function J_func(model::PointFootQuadruped117, env::Environment, q)
    z3 = zeros(3, 3)

    [
        z3   z3 I(3)   z3   z3   z3;
        z3   z3   z3 I(3)   z3   z3;
        z3   z3   z3   z3 I(3)   z3;
        z3   z3   z3   z3   z3 I(3);
    ]
end

function contact_forces(model::PointFootQuadruped117, env::Environment{<:World, LinearizedCone}, γ1, b1, q2, k)
	m = friction_mapping(env)

	SVector{12}([
		m * b1[1:4]; γ1[1];
		m * b1[5:8]; γ1[2];
		m * b1[9:12]; γ1[3];
		m * b1[13:16]; γ1[4];
		])
end

function velocity_stack(model::PointFootQuadruped117, env::Environment{<:World, LinearizedCone}, q1, q2, k, h)
	v = J_func(model, env, q2) * (q2 - q1) / h[1]
	SVector{16}([
		transpose(friction_mapping(env)) * v[1:2];
		transpose(friction_mapping(env)) * v[4:5];
		transpose(friction_mapping(env)) * v[7:8];
		transpose(friction_mapping(env)) * v[10:11];
	])
end

function friction_coefficients(model::PointFootQuadruped117)
	return [model.friction_foot_world]
end

function initialize_z!(z, model::PointFootQuadruped117, idx::RoboDojo.IndicesZ, q)
    z .= 1.0
    z[idx.q] .= q
end

function relative_state_cost(qbody, qorientation, qfoot)
	# cost function on state: 1/2 * qbody'*Qbody*qbody
		# 1/2 * qbody'*Qbody*qbody
		# 1/2 * qorientation'*Qorientation*qorientation
		# 1/2 * (qfoot-qbody)'*Qfoot*(qfoot-qbody)
	Q = zeros(18,18)
	Q[1:3,1:3] = Diagonal(qbody)
	Q[4:6,4:6] = Diagonal(qorientation)
	for i = 1:4
		Q[1:3,1:3] += Diagonal(qfoot)
		Q[3+3i .+ (1:3), 3+3i .+ (1:3)] += Diagonal(qfoot)
		Q[1:3, 3+3i .+ (1:3)] += -Diagonal(qfoot)
		Q[3+3i .+ (1:3), 1:3] += -Diagonal(qfoot)
	end
	return Q
end

function nominal_configuration(model::PointFootQuadruped117)
    [
        0.0; 0.0; model.body_height;
        0.0; 0.0; 0.0;
        model.foot_x; model.foot_y; 0.0;
        model.foot_x;-model.foot_y; 0.0;
       -model.foot_x; model.foot_y; 0.0;
       -model.foot_x;-model.foot_y; 0.0;
    ]
end

function nominal_state(model::PointFootQuadruped117)
	[nominal_configuration(model); zeros(model.nq)]
end

function dynamics(model::PointFootQuadruped117, h, q0, q1, u1, w1, Λ1, q2)
	# evalutate at midpoint
	qm1 = 0.5 * (q0 + q1)
    vm1 = (q1 - q0) / h[1]
    qm2 = 0.5 * (q1 + q2)
    vm2 = (q2 - q1) / h[1]

	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	d = 0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2 # variational integrator (midpoint)
	d .+= transpose(B_fast(model, qm2)) * u1        # control inputs
	d .+= transpose(A_fast(model, qm2)) * w1        # control inputs
	d .+= Λ1                                        # contact impulses

	d[1:3] .-= model.joint_friction .* (vm2[1:3] - vm2[7:9]) # joint friction
	d[1:3] .-= model.joint_friction .* (vm2[1:3] - vm2[10:12]) # joint friction
	d[1:3] .-= model.joint_friction .* (vm2[1:3] - vm2[13:15]) # joint friction
	d[1:3] .-= model.joint_friction .* (vm2[1:3] - vm2[16:18]) # joint friction

	d[7:9]   .-= model.joint_friction .* (vm2[7:9] - vm2[1:3]) # joint friction
	d[10:12] .-= model.joint_friction .* (vm2[10:12] - vm2[1:3]) # joint friction
	d[13:15] .-= model.joint_friction .* (vm2[13:15] - vm2[1:3]) # joint friction
	d[16:18] .-= model.joint_friction .* (vm2[16:18] - vm2[1:3]) # joint friction
	return d
end

# dimensions
nq = 3 + 3 + 3 * 4       # generalized coordinates
nu = 3 * 4               # controls
nw = 3                   # parameters
nc = 4                   # contact points

# parameters
body_height = 0.3
foot_x = 0.17
foot_y = 0.15

# parameters
g = 9.81                 # gravity
friction_foot_world = 0.3
joint_friction = ones(3)
spring_stiffness_joint = 100ones(3)

# inertial properties
mass_body = 13.5
i_xx = 0.0178533
i_xy = 0.0
i_xz = 0.0
i_yz = 0.0
i_yy = 0.0377999
i_zz = 0.0456542
inertia_body = Array(Diagonal([i_xx, i_yy, i_zz]))
mass_foot = 0.2

point_foot_quadruped = PointFootQuadruped117(nq, nu, nw, nc,
				body_height,
				foot_x,
				foot_y,
				mass_body,
                inertia_body,
                mass_foot,
				friction_foot_world,
				joint_friction,
                spring_stiffness_joint,
				g,
				BaseMethods(),
				DynamicsMethods(),
                )
