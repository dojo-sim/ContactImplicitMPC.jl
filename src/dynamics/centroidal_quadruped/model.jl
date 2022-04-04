"""
    centroidal quadruped
    q = (p, r, f1, f2, f3, f4)
        p - body position
        r - body orientation (modified Rodriques parameters)
        f1 - foot 1 position
        f2 - foot 2 position
        f3 - foot 3 position
        f4 - foot 4 position
"""
mutable struct CentroidalQuadruped{T} <: Model{T}
    # dimensions
	nq::Int # generalized coordinates
    nu::Int # controls
    nw::Int # parameters
    nc::Int # contact points

    # environment
    μ_joint::T
    μ_world::T
    g::T

	# parameters
    mass_body::T
    inertia_body::Matrix{T}
    mass_foot::T

	# fast methods
	base
	dyn

	joint_friction
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
function kinematics(model::CentroidalQuadruped, q)
	q[6 .+ (1:12)]
end

lagrangian(model::CentroidalQuadruped, q, q̇) = 0.0

function M_func(model::CentroidalQuadruped, q)
    cat(
        model.mass_body * Diagonal(ones(3)),     # body position
        model.inertia_body,                      # body orienation
        model.mass_foot * Diagonal(ones(3 * 4)), # feet position
        dims=(1, 2)
        )
end

function C_func(model::CentroidalQuadruped, q, q̇)
    [
        model.mass_body * [0,0,model.g];            # body position
        skew(q̇[4:6]) * model.inertia_body * q̇[4:6]; # body orienation
        model.mass_foot * [0,0,model.g];
        model.mass_foot * [0,0,model.g];
        model.mass_foot * [0,0,model.g];
        model.mass_foot * [0,0,model.g];
    ]
end

function ϕ_func(model::CentroidalQuadruped, env::Environment, q)

    position_foot1 = q[6 .+ (1:3)]
    position_foot2 = q[9 .+ (1:3)]
    position_foot3 = q[12 .+ (1:3)]
	position_foot4 = q[15 .+ (1:3)]

	return [position_foot1[3]; position_foot2[3]; position_foot3[3]; position_foot4[3]]
end

function B_func(model::CentroidalQuadruped, q)
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

function A_func(model::CentroidalQuadruped, q)
    @SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
              0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
end

function J_func(model::CentroidalQuadruped, env::Environment, q)
    z3 = zeros(3, 3)

    [
        z3   z3 I(3)   z3   z3   z3;
        z3   z3   z3 I(3)   z3   z3;
        z3   z3   z3   z3 I(3)   z3;
        z3   z3   z3   z3   z3 I(3);
    ]
end

function contact_forces(model::CentroidalQuadruped, env::Environment{<:World, LinearizedCone}, γ1, b1, q2, k)
	m = friction_mapping(env)

	SVector{12}([
		m * b1[1:4]; γ1[1];
		m * b1[5:8]; γ1[2];
		m * b1[9:12]; γ1[3];
		m * b1[13:16]; γ1[4];
		])
end

function velocity_stack(model::CentroidalQuadruped, env::Environment{<:World, LinearizedCone}, q1, q2, k, h)
	v = J_func(model, env, q2) * (q2 - q1) / h[1]
	SVector{16}([
		transpose(friction_mapping(env)) * v[1:2];
		transpose(friction_mapping(env)) * v[4:5];
		transpose(friction_mapping(env)) * v[7:8];
		transpose(friction_mapping(env)) * v[10:11];
	])
end

# dimensions
nq = 3 + 3 + 3 * 4       # generalized coordinates
nu = 3 * 4               # controls
nw = 3                   # parameters
nc = 4                   # contact points

# parameters
g = 9.81                 # gravity
μ_world = 1.0            # coefficient of friction
μ_joint = 1.0            # coefficient of friction

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

centroidal_quadruped = CentroidalQuadruped(nq, nu, nw, nc,
				μ_joint,
				μ_world,
				g,
				mass_body,
                inertia_body,
                mass_foot,
				BaseMethods(), DynamicsMethods(),
                [10ones(3); 30ones(3); 10ones(12)],
                )

function friction_coefficients(model::CentroidalQuadruped)
	return [model.μ_world]
end

function initialize_z!(z, model::CentroidalQuadruped, idx::RoboDojo.IndicesZ, q)
    z .= 1.0
    z[idx.q] .= q
end

# function residual(model::CentroidalQuadruped, env::Environment{<:World,LinearizedCone}, z, θ, κ)
# 	nc = 4
# 	nb = 4

# 	q0 = θ[1:18]
#     q1 = θ[18 .+ (1:18)]
#     u1 = θ[36 .+ (1:12)]
#     w1 = θ[48 .+ (1:3)]
#     μ = θ[51 .+ (1:1)]
#     h = θ[52 .+ (1:1)]

# 	q2 = z[1:18]
#     γ1 = z[18 .+ (1:4)]
#     b1 = z[22 .+ (1:16)]
#     ψ1 = z[38 .+ (1:4)]
#     s1 = z[42 .+ (1:4)]
#     η1 = z[46 .+ (1:16)]
#     s2 = z[62 .+ (1:4)]

# 	ϕ = ϕ_func(model, env, q2)

# 	k = kinematics(model, q2)
# 	λ1 = [
#             friction_mapping(env) * b1[1:4]; γ1[1];
#             friction_mapping(env) * b1[4 .+ (1:4)]; γ1[2];
#             friction_mapping(env) * b1[8 .+ (1:4)]; γ1[3];
#             friction_mapping(env) * b1[12 .+ (1:4)]; γ1[4];
#         ]
# 	Λ1 = transpose(J_func(model, env, q1)) * λ1 #@@@@ maybe need to use J_fast
#     v = J_func(model, env, q2) * (q2 - q1) ./ h[1]

# 	vT_stack = [
#                 v[1:2]; -v[1:2];
#                 v[3:4]; -v[3:4];
#                 v[5:6]; -v[5:6];
#                 v[7:8]; -v[7:8];
#             ]
# 	ψ_stack = [
#         ψ1[1] * ones(4)
#         ψ1[2] * ones(4);
#         ψ1[3] * ones(4);
#         ψ1[4] * ones(4);
#     ]

# 	# @warn "define residual order"
# 	[
#     model.dyn.d(h, q0, q1, u1, w1, Λ1, q2);
# 	 s1 - ϕ;
# 	 η1 - vT_stack - ψ_stack;
# 	 s2 .- (μ[1] * γ1 .- [sum(b1[1:4]); sum(b1[4 .+ (1:4)]); sum(b1[8 .+ (1:4)]); sum(b1[12 .+ (1:4)])]);
# 	 γ1 .* s1 .- κ[1];
# 	 b1 .* η1 .- κ[1];
# 	 ψ1 .* s2 .- κ[1]]
# end

# mutable struct QuadrupedSimple{T} <: Model{T}
# 	nq
#     nu
#     nw
#     nc

# 	g::T
# 	μ_world::T
# 	μ_joint::T

# 	mb::T
# 	mf::T
# 	Ix::T
# 	Iy::T
# 	Iz::T
# 	l_torso::T
# 	w_torso::T

# 	# fast methods
# 	base
# 	dyn

# 	joint_friction
# end

# # Trunk model
# #                             middle   com
# #                                267 mm
# #             o---------------------|---X-----------------o
# #                                   13mm
# # BRhip o___________________________|___________________________o FRhip
# #                    183 mm                    183mm

# # Lagrangian
# function lagrangian(model::QuadrupedSimple, q, q̇)
# 	return 0.0
# end

# function kinematics(model::QuadrupedSimple, q)
# 	pw1 = q[6 .+ (1:3)]
# 	pw2 = q[9 .+ (1:3)]
# 	pw3 = q[12 .+ (1:3)]
# 	pw4 = q[15 .+ (1:3)]
# 	SVector{12}([pw1; pw2; pw3; pw4])
# end

# # Methods
# function M_func(model::QuadrupedSimple, q)
# 	Diagonal(@SVector [model.mb,
# 					   model.mb,
# 					   model.mb,
# 					   model.Ix, model.Iy, model.Iz,
# 					   model.mf, model.mf, model.mf,
# 					   model.mf, model.mf, model.mf,
# 					   model.mf, model.mf, model.mf,
# 					   model.mf, model.mf, model.mf])
# end

# function ϕ_func(model::QuadrupedSimple, env::Environment, q)
# 	k = kinematics(model, q)
# 	@SVector [k[3], k[6], k[9], k[12]]
# end

# function B_func(model::QuadrupedSimple, q)
# 	p_torso = q[1:3]
# 	rot = mrp_rotation_matrix(q[4:6])

# 	# # r in world frame
# 	# r1 = q[6 .+ (1:3)] - p_torso
# 	# r2 = q[9 .+ (1:3)] - p_torso
# 	# r3 = q[12 .+ (1:3)] - p_torso
# 	# r4 = q[15 .+ (1:3)] - p_torso

# 	r1 = p_torso - q[6 .+ (1:3)]
# 	r2 = p_torso - q[9 .+ (1:3)]
# 	r3 = p_torso - q[12 .+ (1:3)]
# 	r4 = p_torso - q[15 .+ (1:3)]
# 	z3 = zeros(3, 3)

# 	transpose(SMatrix{18, 12}([
# 					I I I I;
# 	                transpose(rot) * skew(r1) transpose(rot) * skew(r2) transpose(rot) * skew(r3) transpose(rot) * skew(r4);
# 					-I z3 z3 z3;
# 					z3 -I z3 z3;
# 					z3 z3 -I z3;
# 					z3 z3 z3 -I]))
# end

# function A_func(model::QuadrupedSimple, q)
# 	@SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
# 			  0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
# 			  0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
# end

# function J_func(model::QuadrupedSimple, env::Environment, q)
# 	z3 = zeros(3, 3)
# 	J = transpose([
# 		 zeros(6, 12);
# 	     I z3 z3 z3;
# 	     z3 I z3 z3;
# 	     z3 z3 I z3;
# 	     z3 z3 z3 I])
# end

# function C_func(model::QuadrupedSimple, q, q̇)
# 	SVector{18}([0.0, 0.0, model.g * model.mb,
# 				 cross(q̇[4:6], Diagonal([model.Ix, model.Iy, model.Iz]) * q̇[4:6])...,
# 				 0.0, 0.0, model.g * model.mf,
# 				 0.0, 0.0, model.g * model.mf,
# 				 0.0, 0.0, model.g * model.mf,
# 				 0.0, 0.0, model.g * model.mf])
# end

# function contact_forces(model::QuadrupedSimple, env::Environment{<:World, LinearizedCone}, γ1, b1, q2, k)
# 	m = friction_mapping(env)
# 	SVector{12}([transpose(rotation(env, k[1:2]))   * [m * b1[1:4];   γ1[1]];
# 			     transpose(rotation(env, k[4:5]))   * [m * b1[5:8];   γ1[2]];
# 				 transpose(rotation(env, k[7:8]))   * [m * b1[9:12];  γ1[3]];
# 				 transpose(rotation(env, k[10:11])) * [m * b1[13:16]; γ1[4]]])
# end

# function velocity_stack(model::QuadrupedSimple, env::Environment{<:World, LinearizedCone}, q1, q2, k, h)
# 	v = J_func(model, env, q2) * (q2 - q1) / h[1]

# 	v1_surf = rotation(env, k[1:2]) * v[1:3]
# 	v2_surf = rotation(env, k[4:5]) * v[4:6]
# 	v3_surf = rotation(env, k[7:8]) * v[7:9]
# 	v4_surf = rotation(env, k[10:11]) * v[10:12]

# 	SVector{16}([transpose(friction_mapping(env)) * v1_surf[1:2];
# 				 transpose(friction_mapping(env)) * v2_surf[1:2];
# 				 transpose(friction_mapping(env)) * v3_surf[1:2];
# 				 transpose(friction_mapping(env)) * v4_surf[1:2]])
# end

# function get_stride(model::QuadrupedSimple, traj::ContactTraj)
#     stride = zeros(SizedVector{model.dim.q})
# 	idx = [1,7,10,13,16]
# 	stride[idx] = traj.q[end-1][idx] - traj.q[1][idx]
#     return stride
# end

# ################################################################################
# # Instantiation
# ################################################################################
# # Dimensions
# nq = 3 + 3 + 3 * 4        # configuration dimension
# nu = 4 * 3                # control dimension
# nc = 4                    # number of contact points
# nw = 3
# nquat = 0

# # World parameters
# g = 9.81      # gravity
# μ_world = 0.5 # coefficient of friction
# μ_joint = 0.0 # coefficient of torque friction at the joints

# # # ~Unitree A1
# # mf = 0.4

# #
# # # TRUNK ONLY TODO: parallel axis theorem to add shoulders
# # mb = 4.713
# # Ix = 0.01683993
# # Iy = 0.056579028
# # Iz = 0.064713601
# #
# # l_torso = 0.5 * 0.267 # dimension from com
# # w_torso = 0.5 * 0.194 # dimension from com


# ## Mini Cheetah
# mb = 9.0
# mf = 0.025 * mb
# Ix = 0.07
# Iy = 0.26
# Iz = 0.242

# l_torso = 0.5 * 0.38 # dimension from com
# w_torso = 0.5 * 0.203 # dimension from com

# quadruped_centroidal = QuadrupedSimple(nq, nu, nw, nc,
# 				g, μ_world, μ_joint,
# 				mb, mf, Ix, Iy, Iz, l_torso, w_torso,
# 				BaseMethods(), DynamicsMethods(),
# 				SVector{nq}([zeros(3); μ_joint * ones(nq - 3)]))
