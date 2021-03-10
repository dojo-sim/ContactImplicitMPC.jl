@with_kw struct SimulatorOptions{T}
    warmstart::Bool = true
    z_warmstart::T = 0.001
    κ_warmstart::T = 0.001
end

struct Simulator{S,nq,nu,nc,nb,nw}
    model

    H::Int
    h::S

    q::Vector{SArray{Tuple{nq},S,1,nq}}
    u::Vector{SArray{Tuple{nu},S,1,nu}}
    γ::Vector{SArray{Tuple{nc},S,1,nc}}
    b::Vector{SArray{Tuple{nb},S,1,nb}}
    w::Vector{SArray{Tuple{nw},S,1,nw}}

    dq2dq0::Vector{SizedArray{Tuple{nq,nq},S,2,2}}
    dq2dq1::Vector{SizedArray{Tuple{nq,nq},S,2,2}}
    dq2du::Vector{SizedArray{Tuple{nq,nu},S,2,2}}
    dγdq0::Vector{SizedArray{Tuple{nc,nq},S,2,2}}
    dγdq1::Vector{SizedArray{Tuple{nc,nq},S,2,2}}
    dγdu::Vector{SizedArray{Tuple{nc,nu},S,2,2}}
    dbdq0::Vector{SizedArray{Tuple{nb,nq},S,2,2}}
    dbdq1::Vector{SizedArray{Tuple{nb,nq},S,2,2}}
    dbdu::Vector{SizedArray{Tuple{nb,nu},S,2,2}}

    ip::InteriorPoint{S}
    ip_opts::InteriorPointOptions{S}

    sim_opts::SimulatorOptions{S}
end


function simulator_base(model, q0::SVector, q1::SVector, h::S, H::Int;
    u = [@SVector zeros(model.dim.u) for t = 1:H],
    w = [@SVector zeros(model.dim.w) for t = 1:H],
    ip_opts = InteriorPointOptions{S}(),
    sim_opts = SimulatorOptions{S}()) where S

    simulator(model, q0, q1, h, H;
        u = u,
        w = w,
        ip_opts = ip_opts,
        r! = model.res.r,
        rz! = model.res.rz,
        rθ! = model.res.rθ,
        rz = model.spa.rz_sp,
        rθ = model.spa.rθ_sp,
        sim_opts = sim_opts)
end


function simulator(model, q0::SVector, q1::SVector, h::S, H::Int;
    u = [@SVector zeros(model.dim.u) for t = 1:H],
    w = [@SVector zeros(model.dim.w) for t = 1:H],
    ip_opts = InteriorPointOptions{S}(),
    r! = r!, rz! = rz!, rθ! = rθ!,
    rz = spzeros(num_var(model), num_var(model)),
    rθ = spzeros(num_var(model), num_data(model)),
    sim_opts = SimulatorOptions{S}()) where S

    nq = model.dim.q
    nu = model.dim.u
    nw = model.dim.w
    nc = model.dim.c
    nb = model.dim.b

    q = [q0, q1, [@SVector zeros(nq) for t = 1:H]...]
    γ = [@SVector zeros(nc) for t = 1:H]
    b = [@SVector zeros(nb) for t = 1:H]

    dq2dq0 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:H]
    dq2dq1 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:H]
    dq2du = [SizedMatrix{nq,nu}(zeros(nq, nu)) for t = 1:H]
    dγdq0 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:H]
    dγdq1 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:H]
    dγdu = [SizedMatrix{nc,nu}(zeros(nc, nu)) for t = 1:H]
    dbdq0 = [SizedMatrix{nb,nq}(zeros(nb, nq)) for t = 1:H]
    dbdq1 = [SizedMatrix{nb,nq}(zeros(nb, nq)) for t = 1:H]
    dbdu = [SizedMatrix{nb,nu}(zeros(nb, nu)) for t = 1:H]

    ip = interior_point(
        num_var(model),
        num_data(model),
        inequality_indices(model),
        r! = r!, rz! = rz!, rθ! = rθ!,
        rz = rz,
        rθ = rθ)

    Simulator(
        model,
        H, h,
        q,
        u,
        γ,
        b,
        w,
        dq2dq0,
        dq2dq1,
        dq2du,
        dγdq0,
        dγdq1,
        dγdu,
        dbdq0,
        dbdq1,
        dbdu,
        ip,
        ip_opts,
        sim_opts)
end


function step!(sim::Simulator, t)
    # unpack
    model = sim.model
    q = sim.q
    u = sim.u
    w = sim.w
    h = sim.h
    ip = sim.ip
    z = ip.z
    θ = ip.θ

    # initialize
    if sim.sim_opts.warmstart
        z .+= sim.sim_opts.z_warmstart * rand(ip.num_var)
        sim.ip_opts.κ_init = sim.sim_opts.κ_warmstart
    else
        z_initialize!(z, model, q[t+1])
    end
    θ_initialize!(θ, model, q[t], q[t+1], u[t], w[t], h)

    # solve
    status = interior_point!(ip, opts = sim.ip_opts)

    if status
        # parse result
        q2, γ, b, _ = unpack_z(model, z)
        sim.q[t+2] = copy(q2)
        sim.γ[t] = γ
        sim.b[t] = b

        if sim.ip_opts.diff_sol
            nq = model.dim.q
            nu = model.dim.u
            nc = model.dim.c
            nb = model.dim.b

            sim.dq2dq0[t] = view(ip_data.δz, 1:nq, 1:nq)
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
    - solves 1-step feasibility problem for H time steps
    - initial configurations: q0, q1
    - time step: h
"""
function simulate!(sim::Simulator; verbose = false)

    verbose && println("\nSimulation")

    # initialize configurations for first step
    z_initialize!(sim.ip.z, sim.model, sim.q[2])

    status = true

    # simulate
    for t = 1:sim.H
        verbose && println("t = $t")
        status = step!(sim, t)
        !status && (@error "failed step (t = $t)")
    end

    return status
end
