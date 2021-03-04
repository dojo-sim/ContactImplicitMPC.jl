struct Simulator
    T

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

function simulator(model, q0, q1, h, T;
    u = [@SVector zeros(model.nu) for t = 1:T],
    w = [@SVector zeros(model.nw) for t = 1:T],
    diff_sol = false,
    ip_opts = InteriorPointOptions(diff_sol = diff_sol))

    nq = model.nq
    nu = model.nu
    nw = model.nw
    nc = model.nc
    nb = model.nb

    q = [q0, q1, [@SVector zeros(nq) for t = 1:T]...]
    γ = [@SVector zeros(nc) for t = 1:T]
    b = [@SVector zeros(nb) for t = 1:T]

    dq2q0 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:T]
    dq2dq1 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:T]
    dq2du = [SizedMatrix{nq,nu}(zeros(nq, nu)) for t = 1:T]
    dγdq0 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:T]
    dγdq1 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:T]
    dγdu = [SizedMatrix{nc, nu}(zeros(nc, nu)) for t = 1:T]
    dbdq0 = [SizedMatrix{nb, nq}(zeros(nb, nq)) for t = 1:T]
    dbdq1 = [SizedMatrix{nb, nq}(zeros(nb, nq)) for t = 1:T]
    dbdu = [SizedMatrix{nb, nu}(zeros(nb, nu)) for t = 1:T]

    ip_data = interior_point_data(num_var(model),
        num_data(model),
        inequality_indices(model))
    ip_data.data.info[:model] = model

    Simulator(
        T,
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
    fill!(z, 0.0)
    fill!(θ, 0.0)
    z_initialize!(z, model, q[t+1])
    θ_initialize!(θ, model, q[t], q[t+1], u[t], w[t])

    # solve
    status = interior_point!(ip_data, opts = sim.ip_opts)

    if status
        # parse result
        q2, γ, b, _ = unpack_z(z)
        sim.q[t+2] = copy(q2)
        sim.γ[t] = γ
        sim.b[t] = b

        if sim.ip_opts.diff_sol
            nq = model.nq
            nu = model.nu
            nc = model.nc
            nb = model.nb

            sim.dq2q0[t] = view(ip_data.δz, 1:nq, 1:nq)
            sim.dq2dq1[t] = view(ip_data.δz, 1:nq, nq .+ (1:nq))
            sim.dq2du[t] = view(ip_data.δz, 1:nq, 2 * nq .+ (1:nu))
            sim.dγdq0[t] = view(ip_data.δz, nq .+ (1:nc), 1:nq)
            sim.dγdq1[t] = view(ip_data.δz, nq .+ (1:nc), nq .+ (1:nq))
            sim.dγdu[t] = view(ip_data.δz, nq .+ (1:nc), 2 * nq .+ (1:nu))
            sim.dbdq0[t] = view(ip_data.δz, nq + nc .+ (1:nb), 1:nq)
            sim.dbdq1[t] = view(ip_data.δz, nq + nc .+ (1:nb), nq .+ (1:nq))
            sim.dbdu[t] = view(ip_data.δz, nq + nc .+ (1:nb), 2 * nq .+ (1:nu))
        end
    end
    return status
end

"""
    simulate
    - solves 1-step feasibility problem for T time steps
    - initial configurations: q0, q1
    - time step: h
"""
function simulate!(sim::Simulator; verbose = false)

    println("\nSimulation")

    for t = 1:sim.T
        verbose && println("   t = $t")
        status = step!(sim, t)
        !status && (@error "failed step (t = $t)")
    end
end
