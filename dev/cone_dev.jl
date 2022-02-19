using ForwardDiff
using LinearAlgebra

function cone_prod(u, v)
    w = [u' * v; u[1] * v[2:end] + v[1] * u[2:end]]
    return w
end

function lift(u::AbstractVector{T}) where {T}
    n = length(u)
    U = Array(Diagonal(u[1] * ones(n)))
    U[2:end,1] = u[2:end]
    U[1,2:end] = u[2:end]
    return U
end

function mult(u::AbstractVector{T}, v::AbstractVector{T}) where {T}
    # complexity O(n)
    n = length(u)
    @assert length(v) == n
    w = u[1] * v # n
    w[1] += u[2:end]' * v[2:end] # n-1
    w[2:end] += v[1] * u[2:end] # n-1
    return w
end

function divi(u::AbstractVector{T}, v::AbstractVector{T}) where {T}
    # adapted from https://www.sciencedirect.com/science/article/pii/S0893965914000469
    # complexity O(n)
    n = length(u)
    @assert length(v) == n
    s = u[2:end] ./ u[1] # n-1
    α = - s' * s # n-1
    t = v - 1 / (1 + α) * (v[1] - s' * v[2:end]) * [0.0; s] # 2n-1
    w = t - [s' * t[2:end]; zeros(n-1)] # n
    w = w ./ u[1] # n
    return w
end

# lift
u = [2,3,4]
U = lift(u)
@test U == [2 3 4;
            3 2 0;
            4 0 2;]

# mult
n = 4
u = rand(n)
v = rand(n)
@test norm(mult(u,v) - lift(u) * v) < 1e-10

# div
n = 4
u = rand(n)
v = rand(n)
@test norm(divi(u,v) - inv(lift(u)) * v) < 1e-10


plt = plot()
U = lift(rand(10))
plot!(plt, Gray.(1e10 * abs.(U)))
display(plt)

################################################################################
# Grapical Example of IFT + IP Power
################################################################################

"""
    outer problem
        A small car has to reach a goal point, the outer optimization problem tries to reach this goal point.
        It provides fuel (costly) to the car to reach the destination (reward).

    inner problem
        The inner problem formulates the car dynamics. Fuel is used to move forward.
        Ideally the car has just enough fuel to reach the goal.
        In practice, there is a road block midway preventing the car from getting to the goal.

    solution
        In the optimal solution to this problem, the outer agent chose just the
        right amount of fuel to reach the road block. Providing more fuel is useless
        because the car won't be able to move further and it will cost more.


                start location                      road block      goal location
                *                                   |               *
                                                                         -------
                                                      -------------------
                                    ----------------X-
    fuel cost : --------------------
    goal cost : -------------

                             -----------------------X---------             -----

                                                              -------------
                                                    X Optimal point

"""

# Problem data
x_start = 0.0
x_goal = 2.0
x_block = 1.0

function goal_cost(x)
    return (x - x_goal)^2
end

function fuel_cost(f)
    return f
end

function dynamics(f)
    x = min(x_start + f, x_block)
    return x
end

function outer_cost(x, f)
    return goal_cost(x) + fuel_cost(f)
end

function outer_cost(f::Real)
    return goal_cost(dynamics(f)) + fuel_cost(f)
end

# Finite difference method:
# We use finite difference to ompute the sensitivity of the inner problem.
# This works however it is computationnally intensive, especially in high dimensions.
function grad_outer(f)
    grad = (outer_cost(f + ϵ) - outer_cost(f - ϵ)) / (2ϵ)
    return grad
end

function stepping(f0)
    grad = grad_outer(f0)
    f1 = f0 - α * grad
    return f1
end

function optimize(N)
    plt = plot(layout = (3,1), legend = false)
    f = f0
    for i = 1:N
        f = stepping(f)
        plot!(plt[1,1], F, C)
        plot!(plt[2,1], F, G)
        c = outer_cost(f)
        scatter!(plt[1,1], [f], [c], xlabel = "cost")
        scatter!(plt[3,1], [i], [log(abs(f - x_block))],
            xlims = (0,N),
            ylims = (-8,1),
            xlabel = "dist to optimal")
        display(plt)
    end
    return f
end

f0 = 0.0
ϵ = 0.10
α = 0.05

plt = plot()
F = Vector(range(-0.1, stop = 2.0, length = 1000))
C = [outer_cost(f) for f in F]
G = [grad_outer(f) for f in F]
plot!(plt, F, C)
plot!(plt, F, G)

optimize(25)

# Exact subgradient method:
# This methods relies on the exact subgradients of the inner problem. It should
# stall to a non optimal solution.
function grad_dyn(f) # dx/df
    x = dynamics(f)
    if x < x_block
        return 1.0
    elseif x == x_block
        return 0.0
    elseif x > x_block
        return 0.0
    end
end

function grad_goal(x) # dG/dx
    return 2 * (x - x_goal)
end

function grad_fuel(f) # dF\df
    return 1.0
end

function grad_outer(f)
    x = dynamics(f)
    grad = grad_fuel(f) + grad_goal(x) * grad_dyn(f)
    return grad
end

function stepping(f0)
    grad = grad_outer(f0)
    f1 = f0 - α * grad
    return f1
end

function optimize(N)
    plt = plot(layout = (3,1), legend = false)
    f = f0
    for i = 1:N
        f = stepping(f)
        plot!(plt[1,1], F, C)
        plot!(plt[2,1], F, G)
        c = outer_cost(f)
        scatter!(plt[1,1], [f], [c], xlabel = "cost")
        scatter!(plt[3,1], [i], [log(abs(f - x_block))],
            xlims = (0,N),
            ylims = (-8,1),
            xlabel = "dist to optimal")
        display(plt)
    end
    return f
end

f0 = 0.0
ϵ = 0.10
α = 0.05

plt = plot()
F = Vector(range(-0.1, stop = 2.0, length = 1000))
C = [outer_cost(f) for f in F]
G = [grad_outer(f) for f in F]
plot!(plt, F, C)
plot!(plt, F, G)

optimize(25)


# Smoothed approximation method:
# This method uses a smooth approximation of the subgradient of the inner problem.
# It uses an interior point relaxation of the inner problem. This should work and be more
# computationally efficient than the finite difference method.

function grad_dyn(f) # dx/df
    x = dynamics(f)
    if x < x_block
        return 1.0
    elseif x == x_block
        return 0.0
    elseif x > x_block
        return 0.0
    end
end

function grad_goal(x) # dG/dx
    return 2 * (x - x_goal)
end

function grad_fuel(f) # dF\df
    return 1.0
end

function grad_outer(f)
    x = dynamics(f)
    grad = grad_fuel(f) + grad_goal(x) * grad_dyn(f)
    return grad
end

function stepping(f0)
    grad = grad_outer(f0)
    f1 = f0 - α * grad
    return f1
end

function optimize(N)
    plt = plot(layout = (3,1), legend = false)
    f = f0
    for i = 1:N
        f = stepping(f)
        plot!(plt[1,1], F, C)
        plot!(plt[2,1], F, G)
        c = outer_cost(f)
        scatter!(plt[1,1], [f], [c], xlabel = "cost")
        scatter!(plt[3,1], [i], [log(abs(f - x_block))],
            xlims = (0,N),
            ylims = (-8,1),
            xlabel = "dist to optimal")
        display(plt)
    end
    return f
end

function inner_stepping(y0, z0, f)
    Y0 = [y0, z0]
    res = [1 * (y0 - x_block + x_start + f) - z0, y0 * z0 - ρ]
    jac = [1   -1;
           z0  y0]
    jac_data = [1;
                0]

    α = 1.0
    Δ = - inv(jac) * res
    for i = 1:20
        all(Y0 + α * Δ .> 0) && break
        α = 0.5 * α
    end
    Y = Y0 + α * Δ
    y1, z1 = Y
    grad_Y = - inv(jac) * jac_data
    return y1, z1, norm(res, Inf), grad_Y
end

function inner_optimize(f)
    y = 1.0
    z = 1.0
    for i = 1:100
        y, z, err, grad_Y = inner_stepping(y, z, f)
        grad = - grad_Y[1]
        if err < tol
            x = - y + x_block
            return x, grad
        end
    end
    return nothing
end

function grad_dyn(f)
    return inner_optimize(f)[2]
end

function inner_dynamics(f)
    return inner_optimize(f)[1]
end

function outer_cost(f::Real)
    return goal_cost(inner_dynamics(f)) + fuel_cost(f)
end

function grad_outer(f)
    x = inner_dynamics(f)
    grad = grad_fuel(f) + grad_goal(x) * grad_dyn(f)
    return grad
end


inner_stepping(1.0, 1.0, 0.5)
y = 1.0
z = 1.0
y, z, = inner_stepping(y, z, 0.75)
inner_optimize(0.50)
inner_optimize(0.95)
inner_optimize(0.99)
inner_optimize(1.99)

f0 = 0.0
ϵ = 0.01
α = 0.05
ρ = 1e-3
tol = 1e-8

plt = plot()
F = Vector(range(-0.1, stop = 2.0, length = 1000))
C = [outer_cost(f) for f in F]
G = [grad_outer(f) for f in F]
plot!(plt, F, C)
plot!(plt, F, G)

optimize(25)
