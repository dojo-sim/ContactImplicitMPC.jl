# Base functions
function lagrangian_fast(model::ContactDynamicsModel, q, q̇)
    return model.bas.L(q, q̇)
end

function M_fast(model::ContactDynamicsModel, q)
    return model.bas.M(q)
end

function B_fast(model::ContactDynamicsModel, q)
    return model.bas.B(q)
end

function N_fast(model::ContactDynamicsModel, q)
    return model.bas.N(q)
end

function P_fast(model::ContactDynamicsModel, q)
    return model.bas.P(q)
end

function C_fast(model::ContactDynamicsModel, q, q̇)
    return model.bas.C(q, q̇)
end


# Dynamics functions
function dynamics_fast(model::ContactDynamicsModel, q0, q1, u1, γ1, b1, q2)
    return model.dyn.d(model.dt, q0, q1, u1, γ1, b1, q2)
end

function ∇y_dynamics_fast!(model::ContactDynamicsModel, ∇, q0, q1, u1, γ1, b1, q2)
    return model.dyn.dy(∇, model.dt, q0, q1, u1, γ1, b1, q2)
end

function ∇q0_dynamics_fast!(model::ContactDynamicsModel, ∇, q0, q1, u1, γ1, b1, q2)
    return model.dyn.dq0(∇, model.dt, q0, q1, u1, γ1, b1, q2)
end

function ∇q1_dynamics_fast!(model::ContactDynamicsModel, ∇, q0, q1, u1, γ1, b1, q2)
    return model.dyn.dq1(∇, model.dt, q0, q1, u1, γ1, b1, q2)
end

function ∇u1_dynamics_fast!(model::ContactDynamicsModel, ∇, q0, q1, u1, γ1, b1, q2)
    return model.dyn.du1(∇, model.dt, q0, q1, u1, γ1, b1, q2)
end

function ∇γ1_dynamics_fast!(model::ContactDynamicsModel, ∇, q0, q1, u1, γ1, b1, q2)
    return model.dyn.dγ1(∇, model.dt, q0, q1, u1, γ1, b1, q2)
end

function ∇b1_dynamics_fast!(model::ContactDynamicsModel, ∇, q0, q1, u1, γ1, b1, q2)
    return model.dyn.db1(∇, model.dt, q0, q1, u1, γ1, b1, q2)
end

function ∇q2_dynamics_fast!(model::ContactDynamicsModel, ∇, q0, q1, u1, γ1, b1, q2)
    return model.dyn.dq2(∇, model.dt, q0, q1, u1, γ1, b1, q2)
end


# Residual methods
function r_fast!(model::ContactDynamicsModel, z, θ)
    return model.res.r(z, θ)
end

function rz_fast!(model::ContactDynamicsModel, z, θ)
    return model.res.rz(z, θ)
end

function rθ_fast!(model::ContactDynamicsModel, z, θ)
    return model.res.rθ(z, θ)
end
