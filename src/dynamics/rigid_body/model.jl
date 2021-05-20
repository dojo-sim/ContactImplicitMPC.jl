include(joinpath(pwd(), "src/dynamics/quaternions.jl"))

# rigid body with quaternion representation (sphere w/ radius r)
mutable struct RigidBody{T} <: ContactModel
    dim::Dimensions
    m::T # mass
	J::Vector{T} # inertia
    g::T # gravity
    μ_world::T # friction coefficient
	r::T # radius

	base::BaseMethods
	dyn::DynamicsMethods

	joint_friction::SVector
end

function kinematics(model::RigidBody, q)
	return q[1:3] - [0.0; 0.0; model.r]
end

function lagrangian(model::RigidBody, q, q̇)
	m = 1.0
	J = Diagonal([1.0, 1.0, 1.0])

	p = q[1:3]
	quat = q[4:7]
	v = q̇[1:3]
	ω = q̇[4:6]

	L = 0.0

	# linear
	L += 0.5 * m * transpose(v) * v
	L -= m * model.g * p[3]

	# angular
	L += 0.5 * transpose(ω) * J * ω

	return L
end

# mass matrix
function M_func(model::RigidBody, q)
    m = model.m
	J = model.J
    Diagonal(@SVector [m, m, m, J...])
end

# gravity
function C_func(model::RigidBody, q, q̇)
    m = model.m
    g = model.g

    @SVector [0.0, 0.0, m * g, 0.0, 0.0, 0.0]
end

# signed distance function
function ϕ_func(model::RigidBody, env::Environment, q)
	SVector{1}(q[3] - model.r - env.surf(q[1:2]))
end


# control Jacobian
function B_func(model::RigidBody, q)
	quat = q[4:7]
	r = [0.0; 0.0; -model.r]

    SMatrix{6, 6}([1.0 0.0 0.0 0.0 0.0 0.0;
                   0.0 1.0 0.0 0.0 0.0 0.0;
                   0.0 0.0 1.0 0.0 0.0 0.0;
				   0.0 0.0 0.0 1.0 0.0 0.0;
                   0.0 0.0 0.0 0.0 1.0 0.0;
                   0.0 0.0 0.0 0.0 0.0 1.0])
end

# disturbance Jacobian
function A_func(model::RigidBody, q)
	SMatrix{3, 6}([1.0 0.0 0.0 0.0 0.0 0.0
				   0.0 1.0 0.0 0.0 0.0 0.0
				   0.0 0.0 1.0 0.0 0.0 0.0])
end

# contact Jacobian
function J_func(model::RigidBody, q)
	quat = q[4:7]
	r = [0.0; 0.0; -model.r]
    SMatrix{3, 6}([Diagonal(ones(3)) transpose(skew(r))])
end

function contact_forces(model::RigidBody, γ1, b1, q2, k)
	m = friction_mapping(model.env)

	SVector{3}(transpose(rotation(model.env, k)) * [m * b1; γ1])
end

function velocity_stack(model::RigidBody, q1, q2, k, h)
	p1 = q1[1:3]
	quat1 = q1[4:7]

	p2 = q2[1:3]
	quat2 = q2[4:7]

	v = J_func(model, q2) * [(p2 - p1) / h[1]; ω_finite_difference(quat1, quat2, h[1])]

	v1_surf = rotation(model.env, k) * v

	SVector{4}(friction_mapping(model.env)' * v1_surf[1:2])
end

# quaternion midpoint http://web.cse.ohio-state.edu/~parent.1/classes/682/Lectures/Lectures07/L05_Interpolation/Quaternions.pdf
function dynamics(model::RigidBody, h, q0, q1, u1, w1, λ1, q2)

	p0 = q0[1:3]
	quat0 = q0[4:7]

	p1 = q1[1:3]
	quat1 = q1[4:7]

	p2 = q2[1:3]
	quat2 = q2[4:7]

	# evalutate at midpoint
	qm1 = [0.5 * (p0 + p1); (quat0 + quat1) ./ norm(quat0 + quat1)]
    vm1 = [(p1 - p0) / h[1]; ω_finite_difference(quat0, quat1, h[1])]
    qm2 = [0.5 * (p1 + p2); (quat1 + quat2) ./ norm(quat1 + quat2)]
    vm2 = [(p2 - p1) / h[1]; ω_finite_difference(quat1, quat2, h[1])]

	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	return [(0.5 * h[1] * D1L1 + D2L1 + 0.5 * h[1] * D1L2 - D2L2
		+ transpose(B_fast(model, qm2)) * u1
		+ transpose(A_fast(model, qm2)) * w1
		+ transpose(J_fast(model, q2)) * λ1
		- h[1] * model.joint_friction .* vm2);
		 sqrt(quat2[1]^2.0 + quat2[2]^2.0 + quat2[3]^2.0 + quat2[4]^2.0) - 1.0]
end

function G_func(x)
	q = x[4:7]
	[1.0 0.0 0.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 0.0 0.0;
	 0.0 0.0 1.0 0.0 0.0 0.0;
     0.0 0.0 0.0 -transpose(q[2:4]);
	 zeros(3, 3) Diagonal(q[1] * ones(3)) + skew(q[2:4])]
end

# Model (flat surface)
model = RigidBody(Dimensions(7, 6, 3, 1),
	1.0, [1.0, 1.0, 1.0], 0.0 * 9.81, 0.5, 0.25,
	BaseMethods(), DynamicsMethods(),
	SVector{6}(zeros(6)))

# # Symbolics
# @variables q[1:nq]
# @variables q̇[1:nv]

# nq = model.dim.q
# nv = model.dim.q - 1
#
# lagrangian(model, rand(nq), rand(nv))

# # Lagrangian
# L = lagrangian(model, q, q̇)
# L = Symbolics.simplify.(L)
#
# _dLq = Symbolics.gradient(L, q, simplify=true)
# dLq = G_func(q)' * _dLq
#
# dLq̇ = Symbolics.gradient(L, q̇, simplify=true)
# ddL = Symbolics.hessian(L, [q; q̇], simplify=true)
# ddLq̇q = ddL[nq .+ (1:nv), 1:nq] * G_func(q)
# M = ddL[nq .+ (1:nv), nq .+ (1:nv)]
#
# # Coriolis and Centrifugal forces Jacobians
# C = ddLq̇q * q̇ - dLq
# C = Symbolics.simplify.(C)

dir = joinpath(pwd(), "src/dynamics/rigid_body/")

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "flat/residual.jld2")
path_jac = joinpath(dir, "flat/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "flat/linearized.jld2")

expr_base = generate_base_expressions(model,
	M_analytical = false,
	mapping = G_func, nv = model.dim.q-1)

save_expressions(expr_base, path_base, overwrite=true)
instantiate_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model, path_dyn)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
@load path_jac rz_sp rθ_sp
instantiate_residual!(model, path_res)

model.spa.rz_sp = rz_sp
model.spa.rθ_sp = rθ_sp

expr_linearized = generate_linearized_expressions(model)
save_expressions(expr_linearized, path_linearized, overwrite=true)
instantiate_linearized!(model, path_linearized)

# time
h = 0.01
T = 100

# initial conditions
q1 = @SVector [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0]
q0 = @SVector [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0]

# simulator
sim = ContactControl.simulator(model, q0, q1, h, T,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-6, κ_tol = 1.0e-6,
		diff_sol = false),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)
@test status

include(joinpath(pwd(), "src/dynamics/rigid_body/visuals.jl"))
vis = Visualizer()
render(vis)

visualize!(vis, model, sim.traj.q, Δt = h)

t = 1
qq1 = UnitQuaternion(RotZ(0.0 * π))
qq2 = UnitQuaternion(RotZ(0.5 * π))

qqm = [qq1.w + qq2.w; qq1.x + qq2.x; qq1.y + qq2.y; qq1.z + qq2.z] ./ norm([qq1.w + qq2.w; qq1.x + qq2.x; qq1.y + qq2.y; qq1.z + qq2.z])

settransform!(vis["satellite"],
	  compose(Translation((sim.traj.q[t][1:3] + [-0.25; -0.25; -0.25])...),
			LinearMap(UnitQuaternion(qqm...))))
