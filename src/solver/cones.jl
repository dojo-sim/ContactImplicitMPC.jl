"""
    second-order cone
"""
function second_order_cone_projection(z)
    n = length(z)

    z0 = z[1]
    z1 = view(z, 2:n)

    if norm(z1) <= z0
        return z, true
    elseif norm(z1) <= -z0
        return zero(z), false
    else
        a = 0.5 * (1.0 + z0 / norm(z1))
        z_proj = zero(z)
        z_proj[1:end-1] = a * z1
        z_proj[end] = a * norm(z1)
        return z_proj, false
    end
end

function second_order_cone_product(z, s)
    n = length(z)
    SVector{n}([z' * s; z[1] * view(s, 2:n) + s[1] * view(z, 2:n)])
end

# check that second-order cone constraints are satisfied
function second_order_cone_check(x, idx_soc)
    for idx in idx_soc
        !second_order_cone_projection(view(x, idx))[2] && return false
    end
    return true
end

# check that inequality (non-negative orthant) constraints are satisfied
function inequality_check(x, idx_ineq)
    for i in idx_ineq
        if x[i] <= 0.0
            return false
        end
    end
    return true
end

function cone_check(x, idx_ineq, idx_soc)
    !inequality_check(x, idx_ineq) && return true
    !second_order_cone_check(x, idx_soc) && return true
end
