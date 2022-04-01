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

function skew(x)
    return [0.0  -x[3]  x[2];
            x[3]   0.0 -x[1];
           -x[2]  x[1]   0.0]
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
function quaternion_rotation_matrix(q)
    H = [zeros(1, 3); I(3)]
    transpose(H) * L_mult(q) * transpose(R_mult(q)) * H
end

function quaternion_from_mrp(p)
    """Quaternion (scalar first) from MRP"""
    return (1.0 / (1.0 + dot(p, p))) * [(1 - dot(p, p)); 2.0 * p]
end

function mrp_rotation_matrix(x)
    quaternion_rotation_matrix(quaternion_from_mrp(x))
end

# Kinematics
function kinematics(model::CentroidalQuadruped, q)
	# q[6 .+ (1:12)]
	q[6 .+ (1:3)]
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

	# return [position_foot1[3]; position_foot2[3]; position_foot3[3]; position_foot4[3]]
	return [position_foot1[3]; position_foot2[3]]
    # return [position_foot1[3]]
end

function B_func(model::CentroidalQuadruped, q)
    position_body = q[1:3]
    orientation_body = q[3 .+ (1:3)]
    R = mrp_rotation_matrix(orientation_body)

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
        # z3   z3   z3   z3 I(3)   z3;
        # z3   z3   z3   z3   z3 I(3);
    ]
end

function contact_forces(model::CentroidalQuadruped, env::Environment{<:World, LinearizedCone}, γ1, b1, q2, k)
	m = friction_mapping(env)

	# SVector{12}([
	SVector{6}([
		m * b1[1:4]; γ1[1];
		m * b1[5:8]; γ1[2];
		# m * b1[9:12]; γ1[3];
		# m * b1[13:16]; γ1[4];
		])
end

function velocity_stack(model::CentroidalQuadruped, env::Environment{<:World, LinearizedCone}, q1, q2, k, h)
	v = J_func(model, env, q2) * (q2 - q1) / h[1]
	# SVector{16}([
	SVector{8}([
		friction_mapping(env)' * v[1:2];
		friction_mapping(env)' * v[3:4];
		# friction_mapping(env)' * v[5:6];
		# friction_mapping(env)' * v[7:8];
	])
end

function get_stride(model::CentroidalQuadruped, traj)
    stride = zeros(model.nq)
    stride[1:2] = traj.q[end-1][1:2] - traj.q[1][1:2]
    return stride
end

# dimensions
nq = 3 + 3 + 3 * 4       # generalized coordinates
nu = 3 * 4               # controls
nw = 3                   # parameters
# nc = 4                   # contact points
nc = 2                   # contact points

# parameters
g = 9.81                 # gravity
μ_world = 0.3            # coefficient of friction
μ_joint = 0.0            # coefficient of friction

# inertial properties
mass_body = 1.0
inertia_body = Array(Diagonal(ones(3)))
mass_foot = 0.1

centroidal_quadruped = CentroidalQuadruped(nq, nu, nw, nc,
				μ_joint,
				μ_world,
				g,
				mass_body,
                inertia_body,
                mass_foot,
				BaseMethods(), DynamicsMethods(),
                μ_joint*ones(nq),
                )

function friction_coefficients(model::CentroidalQuadruped)
	# return model.μ_world * ones(4)
	return [model.μ_world]
end

function initialize_z!(z, model::CentroidalQuadruped, idx::RoboDojo.IndicesZ, q)
    z .= 1.0
    z[idx.q] .= q
end
