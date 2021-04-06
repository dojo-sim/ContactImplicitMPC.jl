@with_kw struct SimulatorOptions{T}
    warmstart::Bool = true
    z_warmstart::T = 0.001
    κ_warmstart::T = 0.001
end

struct Simulator{T}
    model::ContactDynamicsModel

    traj::ContactTraj
    deriv_traj::ContactDerivTraj

    p::Policy
    d::Disturbances

    ip::InteriorPoint{T}

    opts::SimulatorOptions{T}
end

function simulator(model, q0::SVector, q1::SVector, h::S, H::Int;
    p = no_policy(model),
    d = no_disturbances(model),
    r! = model.res.r!,
    rz! = model.res.rz!,
    rθ! = model.res.rθ!,
    rz = model.spa.rz_sp,
    rθ = model.spa.rθ_sp,
    ip_opts = InteriorPointOptions{S}(),
    sim_opts = SimulatorOptions{S}()) where S

    # initialize trajectories
    traj = contact_trajectory(H, h, model)
    traj.q[1] = q0
    traj.q[2] = q1
    traj.u[1] = policy(p, traj.q[2], traj, 1)
    traj.w[1] = disturbances(d, traj.q[2], 1)

    traj_deriv = contact_derivative_trajectory(H, model)

    # initialize interior point solver (for pre-factorization)
    z = zeros(num_var(model))
    θ = zeros(num_data(model))
    z_initialize!(z, model, traj.q[2])
    θ_initialize!(θ, model, traj.q[1], traj.q[2], traj.u[1], traj.w[1], h)

    ip = interior_point(z, θ,
        idx_ineq = inequality_indices(model),
        r! = r!,
        rz! = rz!,
        rθ! = rθ!,
        rz = rz,
        rθ = rθ,
        opts = ip_opts)

    Simulator(
        model,
        traj,
        traj_deriv,
        p,
        d,
        ip,
        sim_opts)
end


function step!(sim::Simulator, t)
    # unpack
    model = sim.model
    q = sim.traj.q
    u = sim.traj.u
    w = sim.traj.w
    h = sim.traj.h
    ip = sim.ip
    z = ip.z
    θ = ip.θ

    # policy
    u[t] = policy(sim.p, q[t], sim.traj, t)
    
    # disturbances
    w[t] = disturbances(sim.d, q[t], t)

    # initialize
    if sim.opts.warmstart
        z .+= sim.opts.z_warmstart * rand(ip.num_var)
        sim.ip.opts.κ_init = sim.opts.κ_warmstart
    else
        z_initialize!(z, model, q[t+1])
    end
    θ_initialize!(θ, model, q[t], q[t+1], u[t], w[t], h)

    # solve
    status = interior_point!(ip)

    if status
        # parse result
        q2, γ, b, _ = unpack_z(model, z)
        sim.traj.z[t] = copy(z) # TODO: maybe not use copy
        sim.traj.θ[t] = copy(θ)
        sim.traj.q[t+2] = copy(q2)
        sim.traj.γ[t] = copy(γ)
        sim.traj.b[t] = copy(b)
        sim.traj.κ[1] = ip.κ[1] # the last κ used in the solve.

        if sim.ip.opts.diff_sol
            nq = model.dim.q
            nu = model.dim.u
            nc = model.dim.c
            nb = model.dim.b

            sim.deriv_traj.dq2dq0[t] = view(ip.δz, 1:nq, 1:nq)
            sim.deriv_traj.dq2dq1[t] = view(ip.δz, 1:nq, nq .+ (1:nq))
            sim.deriv_traj.dq2du[t] = view(ip.δz, 1:nq, 2 * nq .+ (1:nu))
            sim.deriv_traj.dγdq0[t] = view(ip.δz, nq .+ (1:nc), 1:nq)
            sim.deriv_traj.dγdq1[t] = view(ip.δz, nq .+ (1:nc), nq .+ (1:nq))
            sim.deriv_traj.dγdu[t] = view(ip.δz, nq .+ (1:nc), 2 * nq .+ (1:nu))
            sim.deriv_traj.dbdq0[t] = view(ip.δz, nq + nc .+ (1:nb), 1:nq)
            sim.deriv_traj.dbdq1[t] = view(ip.δz, nq + nc .+ (1:nb), nq .+ (1:nq))
            sim.deriv_traj.dbdu[t] = view(ip.δz, nq + nc .+ (1:nb), 2 * nq .+ (1:nu))
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
    z_initialize!(sim.ip.z, sim.model, sim.traj.q[2])

    status = true

    # simulate
    for t = 1:sim.traj.H
        verbose && println("t = $t / $(sim.traj.H)")
        status = step!(sim, t)
        !status && (@error "failed step (t = $t)")
    end

    return status
end
