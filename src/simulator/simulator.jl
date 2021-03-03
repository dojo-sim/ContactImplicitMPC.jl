struct SimulatorData
    q
    u
    γ
    b
    w

    dq2q0
    dq2dq1
    dq2du
    dγdq0
    dγdq1
    dγdu
    dbdq0
    dbdq1
    dbdu

    ip_data
    ip_opts
end

function simulator_data(model, q0, q1, h, T;
    u = [@SVector zeros(model.nu) for t = 1:T-1],
    w = [@SVector zeros(model.nw) for t = 1:T-1],
    diff_sol = false,
    ip_opts = InteriorPointOptions(diff_sol = diff_sol))

    q = [q0, q1, [@SVector zeros(model.nq) for t = 1:T-1]]
    γ = [@SVector zeros(model.nc) for t = 1:T-1]
    b = [@SVector zeros(model.nb) for t = 1:T-1]

    dq2q0 = [nothing for t = 1:T-1]
    dq2dq1 = [nothing for t = 1:T-1]
    dq2du = [nothing for t = 1:T-1]
    dγdq0 = [nothing for t = 1:T-1]
    dγdq1 = [nothing for t = 1:T-1]
    dγdu = [nothing for t = 1:T-1]
    dbdq0 = [nothing for t = 1:T-1]
    dbdq1 = [nothing for t = 1:T-1]
    dbdu = [nothing for t = 1:T-1]

    ip_data = interior_point_data(num_var(model),
        num_data(model),
        inequality_indices(model))
    ip_data.data.info[:model] = model

    SimulatorData(
        q,
        u,
        γ,
        b,
        w,
        dq2q0,
        dq2dq1,
        dq2du,
        dγdq0,
        dγdq1,
        dγdu,
        dbdq0,
        dbdq1,
        dbdu,
        ip_data,
        ip_opts)
end

function step!(sim, t)
    # unpack
    q = sim.q
    u = sim.u
    w = sim.w
    ip_data = sim.ip_data
    z = ip_data.z
    θ = ip_data.data.θ
    model = ip_data.data.info[:model]

    # initialize
    z_initialize!(z, model, q[t+1])
    θ_initialize!(θ, model, q[t], q[t+1], u[t], w[t])

    # solve
    interior_point!(sim.ip_data, opts = sim.ip_opts)
end


"""
    simulate
    - solves 1-step feasibility problem for T time steps
    - initial configurations: q0, q1
    - time step: h
"""
function simulate!(sim::SimulatorData)

    println("simulation")

    for t = 1:T-1
        println("   t = $t")
        status = step!()

        if !status
            @error "failed step (t = $t)"
        else
            push!(q, q2)
            push!(γ, γ1)
            push!(b, b1)
            push!(Δq0, _Δq0)
            push!(Δq1, _Δq1)
            push!(Δu,  _Δu1)
        end
    end
end
