"""
    simulate
    - solves 1-step feasibility problem for T time steps
    - initial configurations: q0, q1
    - time step: h
"""
function simulate(q0, q1, T, h;
        u = [],
        w = [],
        r_tol = 1.0e-5,
        κ_tol = 1.0e-5,
        z_init = 1.0,
        κ_init = 1.0)

    println("simulation")

    # initialize histories
    q = [q0, q1]
    γ = []
    b = []
    Δq0 = []
    Δq1 = []
    Δu = []

    # step
    for t = 1:T
        println("   t = $t")
        q2, γ1, b1, _Δq0, _Δq1, _Δu1, status = step(q[end-1], q[end], h,
            u[t], w[t],
            r_tol = r_tol,
            κ_tol = κ_tol,
            z_init = z_init,
            κ_init = κ_init)

        if !status
            @error "failed step (t = $t)"
            return q, γ, b, Δq0, Δq1, Δu
        else
            push!(q, q2)
            push!(γ, γ1)
            push!(b, b1)
            push!(Δq0, _Δq0)
            push!(Δq1, _Δq1)
            push!(Δu,  _Δu1)
        end
    end

    return q, γ, b, Δq0, Δq1, Δu
end
