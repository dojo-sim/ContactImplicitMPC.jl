struct Simulator{S,nq,nu,nc,nb,nw}
    model

    T::Int

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
end

function simulator(model, q0::SVector, q1::SVector, h::S, T::Int;
    u = [@SVector zeros(model.nu) for t = 1:T],
    w = [@SVector zeros(model.nw) for t = 1:T],
    ip_opts = InteriorPointOptions(),
    rz = spzeros(num_var(model), num_var(model)),
    rθ = spzeros(num_var(model), num_data(model))) where S

    nq = model.nq
    nu = model.nu
    nw = model.nw
    nc = model.nc
    nb = model.nb

    q = [q0, q1, [@SVector zeros(nq) for t = 1:T]...]
    γ = [@SVector zeros(nc) for t = 1:T]
    b = [@SVector zeros(nb) for t = 1:T]

    dq2dq0 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:T]
    dq2dq1 = [SizedMatrix{nq,nq}(zeros(nq, nq)) for t = 1:T]
    dq2du = [SizedMatrix{nq,nu}(zeros(nq, nu)) for t = 1:T]
    dγdq0 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:T]
    dγdq1 = [SizedMatrix{nc,nq}(zeros(nc, nq)) for t = 1:T]
    dγdu = [SizedMatrix{nc,nu}(zeros(nc, nu)) for t = 1:T]
    dbdq0 = [SizedMatrix{nb,nq}(zeros(nb, nq)) for t = 1:T]
    dbdq1 = [SizedMatrix{nb,nq}(zeros(nb, nq)) for t = 1:T]
    dbdu = [SizedMatrix{nb,nu}(zeros(nb, nu)) for t = 1:T]

    ip = interior_point(
        num_var(model),
        num_data(model),
        inequality_indices(model),
        rz = rz,
        rθ = rθ)

    Simulator(
        model,
        T,
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
        ip_opts)
end

function step!(sim, t)
    # unpack
    model = sim.model
    q = sim.q
    u = sim.u
    w = sim.w
    ip = sim.ip
    z = ip.z
    θ = ip.θ

    # initialize
    fill!(z, 0.0)
    fill!(θ, 0.0)
    z_initialize!(z, model, q[t+1])
    θ_initialize!(θ, model, q[t], q[t+1], u[t], w[t])

    # solve
    status = interior_point!(ip, opts = sim.ip_opts)

    if status
        # parse result
        q2, γ, b, _ = unpack_z(z, model)
        sim.q[t+2] = copy(q2)
        sim.γ[t] = γ
        sim.b[t] = b

        if sim.ip_opts.diff_sol
            nq = model.nq
            nu = model.nu
            nc = model.nc
            nb = model.nb

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
    - solves 1-step feasibility problem for T time steps
    - initial configurations: q0, q1
    - time step: h
"""
function simulate!(sim::Simulator; verbose = false)

    verbose && println("\nSimulation")

    for t = 1:sim.T
        verbose && println("t = $t")
        status = step!(sim, t)
        !status && (@error "failed step (t = $t)")
    end
end
