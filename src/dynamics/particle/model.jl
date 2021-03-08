"""
    particle dynamics
    - 3D particle subject to contact forces

    - configuration: q = (x, y, z) ∈ R³
    - impact force (magnitude): γ ∈ R₊
    - friction force: β ∈ R⁴₊
        - friction coefficient: μ ∈ R₊

    Discrete Mechanics and Variational Integrators
        pg. 363
"""
struct Particle{T} <: ContactDynamicsModel
    dim::Dimensions
    m::T # mass
    g::T # gravity
    μ_world::T # friction coefficient
	μ_joint::T

	base::BaseMethods
	dyn::DynamicsMethods
	res::ResidualMethods

	joint_friction::SVector
end

function lagrangian(model::Particle, q, q̇)
	L = 0.0

	L += 0.5 * model.m * transpose(q̇) * q̇
	L -= model.m * model.g * q[3]

	return L
end

# mass matrix
function M_func(model, q)
    m = model.m

    Diagonal(@SVector [m, m, m])
end

# gravity
function C_func(model, q, q̇)
    m = model.m
    g = model.g

    @SVector [0.0, 0.0, m * g]
end

# signed distance function
function ϕ_func(model, q)
    q[3:3]
end

# control Jacobian
function B_func(model, q)
    SMatrix{3, 3}([1.0 0.0 0.0;
                   0.0 1.0 0.0;
                   0.0 0.0 1.0])
end

# normal Jacobian
function N_func(model, q)
    SMatrix{1, 3}([0.0, 0.0, 1.0])
end

# tangent Jacobian
function P_func(model, q)
    SMatrix{4, 3}([1.0 0.0 0.0;
                   0.0 1.0 0.0;
                   -1.0 0.0 0.0;
                   0.0 -1.0 0.0])
end

# Model
particle = Particle(Dimensions(3, 3, 3, 1, 4), 1.0, 9.81, 1.0, 0.0,
	BaseMethods(), DynamicsMethods(), ResidualMethods(),
	@SVector zeros(3))

# path_base = "base.jld2"
# path_dyn = "dynamics.jld2"
# path_res = "residual.jld2"
# path_jac = "sparse_jacobians.jld2"
#
# expr_base = generate_base_expressions(model)
# save_expressions(expr_base, path_base, overwrite=true)
# instantiate_base!(model, path_base)
#
# expr_dyn = generate_dynamics_expressions(model)
# save_expressions(expr_dyn, path_dyn, overwrite=true)
# instantiate_dynamics!(model, path_dyn)
#
# expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
# save_expressions(expr_res, path_res, overwrite=true)
# @save path_jac rz_sp rθ_sp
# instantiate_residual!(model, path_res)
