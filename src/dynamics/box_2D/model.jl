"""
    2D 4-corner box
    - 2D box subject to contact forces

    - configuration: q = (x, z, θ) ∈ R³
    - impact force (magnitude): γ ∈ R₊
    - friction force: β ∈ R²₊
        - friction coefficient: μ ∈ R₊
"""
mutable struct Box2D{T} <: ContactModel
    dim::Dimensions
    m::T # mass
    g::T # gravity
	J::T # inertia about y axis
    μ_world::T # friction coefficient
	μ_joint::T

	n_corners::Int
	corner_offset

	base::BaseMethods
	dyn::DynamicsMethods

	joint_friction::SVector
end

function lagrangian(model::Box2D, q, q̇)
	L = 0.0

	L += 0.5 * model.m * transpose(q̇[1:2]) * q̇[1:2]
	L += 0.5 * model.J * q̇[3]^2
	L -= model.m * model.g * q[2]

	return L
end

function kinematics(model::Box2D, q)
	θ = q[3]
	R = [ cos(θ) sin(θ);
		 -sin(θ) cos(θ)]
	p = q[1:2]
	p1 = p + R * model.corner_offset[1]
	p2 = p + R * model.corner_offset[2]
	p3 = p + R * model.corner_offset[3]
	p4 = p + R * model.corner_offset[4]
    SVector{2 * model.n_corners}([p1;
	                              p2;
								  p3;
								  p4])
end

# mass matrix
function M_func(model::Box2D, q)
	m = model.m
    J = model.J

    Diagonal(@SVector [m, m, J])
end

# gravity
function C_func(model::Box2D, q, q̇)
    m = model.m
    g = model.g

    @SVector [0.0, m * g, 0.0]
end

# signed distance function
function ϕ_func(model::Box2D, env::Environment, q)
	k = kinematics(model, q)
	SVector{model.n_corners}([k[2] - env.surf(k[1:1]),
							  k[4] - env.surf(k[3:3]),
							  k[6] - env.surf(k[5:5]),
							  k[8] - env.surf(k[7:7])])
end

# control Jacobian
function B_func(model::Box2D, q)
    SMatrix{3, 3}([1.0 0.0 0.0;
				   0.0 1.0 0.0;
                   0.0 0.0 1.0])
end

# disturbance Jacobian
function A_func(model::Box2D, q)
	SMatrix{3, 3}([1.0 0.0 0.0;
				   0.0 1.0 0.0;
                   0.0 0.0 1.0])
end

# contact Jacobian
function J_func(model::Box2D, env::Environment, q)
	k(z) = kinematics(model, z)
	ForwardDiff.jacobian(k, q)
end

function contact_forces(model::Box2D, env::Environment{<:World, LinearizedCone}, γ1, b1, q2, k)
	m = friction_mapping(env)

	SVector{8}([transpose(rotation(env, k[1:2])) * [m * b1[1:2]; γ1[1]];
				transpose(rotation(env, k[3:4])) * [m * b1[3:4]; γ1[2]];
				transpose(rotation(env, k[5:6])) * [m * b1[5:6]; γ1[3]];
				transpose(rotation(env, k[7:8])) * [m * b1[7:8]; γ1[4]]])
end

function contact_forces(model::Box2D, env::Environment{<:World, NonlinearCone}, γ1, b1, q2, k)
	SVector{8}([transpose(rotation(env, k[1:2])) * [b1[1:1]; γ1[1]];
				transpose(rotation(env, k[3:4])) * [b1[2:2]; γ1[2]];
				transpose(rotation(env, k[5:6])) * [b1[3:3]; γ1[3]];
				transpose(rotation(env, k[7:8])) * [b1[4:5]; γ1[4]]])
end

function velocity_stack(model::Box2D, env::Environment{<:World, LinearizedCone}, q1, q2, k, h)
	v = J_func(model, env, q2) * (q2 - q1) / h[1]

	v1_surf = rotation(env, k[1:2]) * v[1:2]
	v2_surf = rotation(env, k[3:4]) * v[3:4]
	v3_surf = rotation(env, k[5:6]) * v[5:6]
	v4_surf = rotation(env, k[7:8]) * v[7:8]
	SVector{8}([transpose(friction_mapping(env)) * v1_surf[1:1];
	            transpose(friction_mapping(env)) * v2_surf[1:1];
				transpose(friction_mapping(env)) * v3_surf[1:1];
				transpose(friction_mapping(env)) * v4_surf[1:1]])
end

function velocity_stack(model::Box2D, env::Environment{<:World, NonlinearCone}, q1, q2, k, h)
	nc = model.dim.c
	ne = dim(env)
	v = J_func(model, env, q2) * (q2 - q1) / h[1]

	v_surf = [rotation(env, k[(i-1) * ne .+ (1:ne)]) * v[(i-1) * ne .+ (1:ne)] for i = 1:nc]

	vT_stack = vcat([v_surf[i][1:ne-1] for i = 1:nc]...)
end



# Kinematics
r = 0.5
c1 = @SVector [r, r]
c2 = @SVector [r, -r]
c3 = @SVector [-r, r]
c4 = @SVector [-r, -r]

corner_offset = @SVector [c1, c2, c3, c4]

# Model (flat surface)
box_2D = Box2D(Dimensions(3, 3, 3, 4, 0), 1.0, 9.81, 1/6, 0.5, 0.0, 4, corner_offset,
	BaseMethods(), DynamicsMethods(),
	SVector{3}(0.05*ones(3)))
