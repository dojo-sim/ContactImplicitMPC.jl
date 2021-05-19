mutable struct QuadrupedSimple{T} <: ContactDynamicsModel
	dim::Dimensions

	g::T
	μ_world::T
	μ_joint::T

	mb::T
	mf::T
	Ix::T
	Iy::T
	Iz::T
	l_torso::T
	w_torso::T

	alt

	# fast methods
	base
	dyn
	con
	res
	linearized

	spa::SparseStructure

	joint_friction

	env::Environment
end

# Trunk model
#                             middle   com
#                                267 mm
#             o---------------------|---X-----------------o
#                                   13mm
# BRhip o___________________________|___________________________o FRhip
#                    183 mm                    183mm

# Lagrangian
function lagrangian(model::QuadrupedSimple, q, q̇)
	return 0.0
end

function kinematics(model::QuadrupedSimple, q)
	pw1 = q[6 .+ (1:3)]
	pw2 = q[9 .+ (1:3)]
	pw3 = q[12 .+ (1:3)]
	pw4 = q[15 .+ (1:3)]
	SVector{12}([pw1; pw2; pw3; pw4])
end

# Methods
function M_func(model::QuadrupedSimple, q)
	Diagonal(@SVector [model.mb,
					   model.mb,
					   model.mb,
					   model.Ix, model.Iy, model.Iz,
					   model.mf, model.mf, model.mf,
					   model.mf, model.mf, model.mf,
					   model.mf, model.mf, model.mf,
					   model.mf, model.mf, model.mf])
end

function ϕ_func(model::QuadrupedSimple, q)
	k = kinematics(model, q)
	@SVector [k[3], k[6], k[9], k[12]]
end

function B_func(model::QuadrupedSimple, q)
	p_torso = q[1:3]
	rot = MRP(q[4:6]...)

	# # r in world frame
	# r1 = q[6 .+ (1:3)] - p_torso
	# r2 = q[9 .+ (1:3)] - p_torso
	# r3 = q[12 .+ (1:3)] - p_torso
	# r4 = q[15 .+ (1:3)] - p_torso

	r1 = p_torso - q[6 .+ (1:3)]
	r2 = p_torso - q[9 .+ (1:3)]
	r3 = p_torso - q[12 .+ (1:3)]
	r4 = p_torso - q[15 .+ (1:3)]
	z3 = zeros(3, 3)

	transpose(SMatrix{18, 12}([
					I I I I;
	                transpose(rot) * skew(r1) transpose(rot) * skew(r2) transpose(rot) * skew(r3) transpose(rot) * skew(r4);
					-I z3 z3 z3;
					z3 -I z3 z3;
					z3 z3 -I z3;
					z3 z3 z3 -I]))
end

function A_func(model::QuadrupedSimple, q)
	@SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
			  0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
end

function J_func(model::QuadrupedSimple, q)
	z3 = zeros(3, 3)
	J = transpose([
		 zeros(6, 12);
	     I z3 z3 z3;
	     z3 I z3 z3;
	     z3 z3 I z3;
	     z3 z3 z3 I])
end

function C_func(model::QuadrupedSimple, q, q̇)
	SVector{18}([0.0, 0.0, model.g * model.mb,
				 0.0 * cross(q̇[4:6], Diagonal([model.Ix, model.Iy, model.Iz]) * q̇[4:6])...,
				 0.0, 0.0, model.g * model.mf,
				 0.0, 0.0, model.g * model.mf,
				 0.0, 0.0, model.g * model.mf,
				 0.0, 0.0, model.g * model.mf])
end

function contact_forces(model::QuadrupedSimple, γ1, b1, q2, k)
	m = friction_mapping(model.env)
	SVector{12}([transpose(rotation(model.env, k[1:2]))   * [m * b1[1:4];   γ1[1]];
			     transpose(rotation(model.env, k[4:5]))   * [m * b1[5:8];   γ1[2]];
				 transpose(rotation(model.env, k[7:8]))   * [m * b1[9:12];  γ1[3]];
				 transpose(rotation(model.env, k[10:11])) * [m * b1[13:16]; γ1[4]]])
end

function velocity_stack(model::QuadrupedSimple, q1, q2, k, h)
	v = J_func(model, q2) * (q2 - q1) / h[1]

	v1_surf = rotation(model.env, k[1:2]) * v[1:3]
	v2_surf = rotation(model.env, k[4:5]) * v[4:6]
	v3_surf = rotation(model.env, k[7:8]) * v[7:9]
	v4_surf = rotation(model.env, k[10:11]) * v[10:12]

	SVector{16}([friction_mapping(model.env)' * v1_surf[1:2];
				 friction_mapping(model.env)' * v2_surf[1:2];
				 friction_mapping(model.env)' * v3_surf[1:2];
				 friction_mapping(model.env)' * v4_surf[1:2]])
end

function get_stride(model::QuadrupedSimple, traj::ContactTraj)
    stride = zeros(SizedVector{model.dim.q})
	idx = [1,7,10,13,16]
	stride[idx] = traj.q[end-1][idx] - traj.q[1][idx]
    return stride
end

################################################################################
# Instantiation
################################################################################
# Dimensions
nq = 3 + 3 + 3 * 4        # configuration dimension
nu = 4 * 3                # control dimension
nc = 4                    # number of contact points
nf = 4                    # number of parameters for friction cone
nb = nc * nf
nw = 3

# World parameters
g = 9.81      # gravity
μ_world = 0.5 # coefficient of friction
μ_joint = 0.0 # coefficient of torque friction at the joints

# # ~Unitree A1
# mf = 0.4

#
# # TRUNK ONLY TODO: parallel axis theorem to add shoulders
# mb = 4.713
# Ix = 0.01683993
# Iy = 0.056579028
# Iz = 0.064713601
#
# l_torso = 0.5 * 0.267 # dimension from com
# w_torso = 0.5 * 0.194 # dimension from com


## Mini Cheetah
mb = 9.0
mf = 0.025 * mb
Ix = 0.07
Iy = 0.26
Iz = 0.242

l_torso = 0.5 * 0.38 # dimension from com
w_torso = 0.5 * 0.203 # dimension from com

quadrupedlinear = QuadrupedSimple(Dimensions(nq, nu, nw, nc, nb),
				g, μ_world, μ_joint,
				mb, mf, Ix, Iy, Iz, l_torso, w_torso,
				zeros(nc),
				BaseMethods(), DynamicsMethods(), ContactMethods(),
				ResidualMethods(), ResidualMethods(),
				SparseStructure(spzeros(0, 0), spzeros(0, 0)),
				SVector{nq}([zeros(3); μ_joint * ones(nq - 3)]),
				environment_3D_flat())
