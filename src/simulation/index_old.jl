using ForwardDiff
using LinearAlgebra

function lift(u)
    n = length(u)
    U = Array(Diagonal(u[1] * ones(n)))
    U[2:end,1] = u[2:end]
    U[1,2:end] = u[2:end]
    return U
end

function cone_prod(u, v)
    w = [u' * v; u[1] * v[2:end] + v[1] * u[2:end]]
    return w
end

function lift_diff(u)
    n = length(u)
    v = rand(n)
    cp(v) = cone_prod(u, v)
    U = ForwardDiff.jacobian(cp, u)
    return U
end

n = 5
v = [rand(); 1e-1 * rand(n-1)]

lift(v) - lift_diff(v)
lift(v)
qr(lift(v))
inv(lift(rand(10)))



function inv_lift(u)
    plt = plot()

    n = length(u)
    # permutation
    P = zeros(n, n)
    for i = 1:n
        P[end-i+1, i] = 1.0
    end
    plot!(plt, Gray.(1e10 * abs.(P)))
    display(plt)
    # Lift and perm
    U = lift(u)
    plot!(plt, Gray.(1e10 * abs.(U)))
    display(plt)
    Up = P * U * P
    plot!(plt, Gray.(1e10 * abs.(Up)))
    display(plt)
    # Lower triangular
    S1 = zeros(n, n)
    S1[end, 1:end-1] = Up[end, 1:end-1]
    plot!(plt, Gray.(1e10 * abs.(S1)))
    display(plt)

    # Upper triangular
    S2 = zeros(n, n)
    S2[1:end-1, end] = Up[1:end-1, end]
    plot!(plt, Gray.(1e10 * abs.(S2)))
    display(plt)

    α = sum([-Up[i, end] * Up[end, i] for i=1:n-1])
    @show α
    Upi = (I - S1) * (I - 1 / (1 + α) * (S2 * (I - S1)))
    @show norm(Upi - inv(Up))
    Ui = P * Upi  * P
    return Ui
end
inv_lift(v)
inv_lift(v) - inv(lift(v))

lift(rand(5)) * lift(rand(5))



function inv_dot(u, a)
    # b = inv(U) * a
    n = length(u)

    # permutation
    P = zeros(n, n)
    for i = 1:n
        P[end-i+1, i] = 1.0
    end

    # ap = P * a
    # up = P * u
    # α = - up[1:end-1]' * up[1:end-1]
    #
    # tp = ap - 1 / (1 + α) * (ap[end] * [up[1:end-1]; 0.0] - up[1:end-1]' * ap[1:end-1] * [up[1:end-1]; 0.0])
    # bp = tp - [zeros(n-1); up[1:end-1]' * tp[1:end-1]]
    # b = P * b

    α = - u[2:end]' * u[2:end]

    t = a - 1 / (1 + α) * (a[1] * [0.0; u[2:end]] - u[2:end]' * a[2:end] * [0.0; u[2:end]])
    b = t - [u[2:end]' * t[2:end]; zeros(n-1)]

    # testing
    U = lift(u)
    Ui = inv(U)
    bt = Ui * a
    @show norm(bt - b)

    return b
end

inv_dot([1.0; rand(4)], rand(5))
