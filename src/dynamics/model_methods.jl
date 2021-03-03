function lagrangian_fast(model::ContactDynamicsModel, q, q̇)
    try
        return model.fct.L(q, q̇)
    catch e
        if isa(e, MethodError)
            println("You need to instantiate the fast dynamics methods,
                using instantiate_dynamics!(model, path)")
        end
    end
end

function M_fast(model::ContactDynamicsModel, q)
    try
        return model.fct.M(q)
    catch e
        if isa(e, MethodError)
            println("You need to instantiate the fast dynamics methods,
                using instantiate_dynamics!(model, path)")
        end
    end
end

function B_fast(model::ContactDynamicsModel, q)
    try
        return model.fct.B(q)
    catch e
        if isa(e, MethodError)
            println("You need to instantiate the fast dynamics methods,
                using instantiate_dynamics!(model, path)")
        end
    end
end

function N_fast(model::ContactDynamicsModel, q)
    try
        return model.fct.N(q)
    catch e
        if isa(e, MethodError)
            println("You need to instantiate the fast dynamics methods,
                using instantiate_dynamics!(model, path)")
        end
    end
end

function P_fast(model::ContactDynamicsModel, q)
    try
        return model.fct.P(q)
    catch e
        if isa(e, MethodError)
            println("You need to instantiate the fast dynamics methods,
                using instantiate_dynamics!(model, path)")
        end
    end
end

function C_fast(model::ContactDynamicsModel, q, q̇)
    try
        return model.fct.C(q, q̇)
    catch e
        if isa(e, MethodError)
            println("You need to instantiate the fast dynamics methods,
                using instantiate_dynamics!(model, path)")
        end
    end
end

function dynamics_fast(model::ContactDynamicsModel, q_1, q, u, γ, b, q1)
    try
        return model.fct.d(model.dt, q_1, q, u, γ, b, q1)
    catch e
        if isa(e, MethodError)
            println("You need to instantiate the fast dynamics methods,
                using instantiate_dynamics!(model, path)")
        end
    end
end

function ∇z_dynamics_fast!(model::ContactDynamicsModel, ∇, q_1, q, u, γ, b, q1)
    try
        return model.fct.dz(∇, model.dt, q_1, q, u, γ, b, q1)
    catch e
        if isa(e, MethodError)
            println("You need to instantiate the fast dynamics methods,
                using instantiate_dynamics!(model, path)")
        end
    end
end

function ∇q_1_dynamics_fast!(model::ContactDynamicsModel, ∇, q_1, q, u, γ, b, q1)
    try
        return model.fct.dq_1(∇, model.dt, q_1, q, u, γ, b, q1)
    catch e
        if isa(e, MethodError)
            println("You need to instantiate the fast dynamics methods,
                using instantiate_dynamics!(model, path)")
        end
    end
end

function ∇q_dynamics_fast!(model::ContactDynamicsModel, ∇, q_1, q, u, γ, b, q1)
    try
        return model.fct.dq(∇, model.dt, q_1, q, u, γ, b, q1)
    catch e
        if isa(e, MethodError)
            println("You need to instantiate the fast dynamics methods,
                using instantiate_dynamics!(model, path)")
        end
    end
end

function ∇u_dynamics_fast!(model::ContactDynamicsModel, ∇, q_1, q, u, γ, b, q1)
    try
        return model.fct.du(∇, model.dt, q_1, q, u, γ, b, q1)
    catch e
        if isa(e, MethodError)
            println("You need to instantiate the fast dynamics methods,
                using instantiate_dynamics!(model, path)")
        end
    end
end

function ∇γ_dynamics_fast!(model::ContactDynamicsModel, ∇, q_1, q, u, γ, b, q1)
    try
        return model.fct.dγ(∇, model.dt, q_1, q, u, γ, b, q1)
    catch e
        if isa(e, MethodError)
            println("You need to instantiate the fast dynamics methods,
                using instantiate_dynamics!(model, path)")
        end
    end
end

function ∇b_dynamics_fast!(model::ContactDynamicsModel, ∇, q_1, q, u, γ, b, q1)
    # try
        return model.fct.db(∇, model.dt, q_1, q, u, γ, b, q1)
    # catch e
    #     if isa(e, MethodError)
    #         println("You need to instantiate the fast dynamics methods,
    #             using instantiate_dynamics!(model, path)")
    #     end
    # end
end

function ∇q1_dynamics_fast!(model::ContactDynamicsModel, ∇, q_1, q, u, γ, b, q1)
    # try
        return model.fct.dq1(∇, model.dt, q_1, q, u, γ, b, q1)
    # catch e
    #     if isa(e, MethodError)
    #         println("You need to instantiate the fast dynamics methods,
    #             using instantiate_dynamics!(model, path)")
    #     end
    # end
end
