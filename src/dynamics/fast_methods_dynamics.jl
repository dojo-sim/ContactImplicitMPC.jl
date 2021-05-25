function M_fast(model::ContactModel, q)
    return model.base.M(q)
end

function B_fast(model::ContactModel, q)
    return model.base.B(q)
end

function A_fast(model::ContactModel, q)
    return model.base.A(q)
end

function J_fast(model::ContactModel, q)
    return model.base.J(q)
end

function C_fast(model::ContactModel, q, q̇)
    return model.base.C(q, q̇)
end

# Dynamics functions
function d_fast(model::ContactModel, h, q0, q1, u1, w1, λ1, q2)
    return model.dyn.d(h, q0, q1, u1, w1, λ1, q2)
end
