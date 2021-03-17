function ϕ_fast(model::ContactDynamicsModel, q)
    return model.base.ϕ(Array(q))
end

function M_fast(model::ContactDynamicsModel, q)
    return model.base.M(q)
end

function B_fast(model::ContactDynamicsModel, q)
    return model.base.B(q)
end

function A_fast(model::ContactDynamicsModel, q)
    return model.base.A(q)
end

function J_fast(model::ContactDynamicsModel, q)
    return model.base.J(q)
end

function C_fast(model::ContactDynamicsModel, q, q̇)
    return model.base.C(q, q̇)
end

# Dynamics functions
function d_fast(model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)
    return model.dyn.d(h, q0, q1, u1, w1, γ1, b1, q2)
end

function dy_fast!(v, model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)
    return model.dyn.dy(v, h, q0, q1, u1, w1, γ1, b1, q2)
end

function dq0_fast!(v, model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)
    return model.dyn.dq0(v, h, q0, q1, u1, w1, γ1, b1, q2)
end

function dq1_fast!(v, model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)
    return model.dyn.dq1(v, h, q0, q1, u1, w1, γ1, b1, q2)
end

function du1_fast!(v, model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)
    return model.dyn.du1(v, h, q0, q1, u1, w1, γ1, b1, q2)
end

function dw1_fast!(v, model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)
    return model.dyn.dw1(v, h, q0, q1, u1, w1, γ1, b1, q2)
end

function dγ1_fast!(v, model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)
    return model.dyn.dγ1(v, h, q0, q1, u1, w1, γ1, b1, q2)
end

function db1_fast!(v, model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)
    return model.dyn.db1(v, h, q0, q1, u1, w1, γ1, b1, q2)
end

function dq2_fast!(v, model::ContactDynamicsModel, h, q0, q1, u1, w1, γ1, b1, q2)
    return model.dyn.dq2(v, h, q0, q1, u1, w1, γ1, b1, q2)
end

# Residual methods
function r_fast!(v, model::ContactDynamicsModel, z, θ, κ)
    return model.res.r(v, z, θ, κ)
end

function rz_fast!(v, model::ContactDynamicsModel, z, θ)
    return model.res.rz(v, z, θ)
end

function rθ_fast!(v, model::ContactDynamicsModel, z, θ)
    return model.res.rθ(v, z, θ)
end
