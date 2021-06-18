@with_kw struct SimulatorOptions{T}
    warmstart::Bool = true
    z_warmstart::T = 0.001
    κ_warmstart::T = 0.001
end

struct Simulator{T}
    s::Simulation

    traj::ContactTraj
    deriv_traj::ContactDerivTraj

    p::Policy
    uL::Vector{T}
    uU::Vector{T}

    d::Disturbances

    ip::AbstractIPSolver

    opts::SimulatorOptions{T}
end

function simulator(s::Simulation, q0::SVector, q1::SVector, h::S, H::Int;
    p = no_policy(s.model),
    uL = -Inf * ones(s.model.dim.u),
    uU = Inf * ones(s.model.dim.u),
    d = no_disturbances(s.model),
    r! = s.res.r!,
    rm! = s.res.rm!,
    rz! = s.res.rz!,
    rθ! = s.res.rθ!,
    rz = s.rz,
    rθ = s.rθ,
    space = Euclidean(num_var(s.model, s.env)),
    ip_type::Symbol = :interior_point, # mehrotra
    ip_opts = interior_point_options(ip_type, T=S),
    sim_opts = SimulatorOptions{S}()) where S

    model = s.model
    env = s.env

    # initialize trajectories
    traj = contact_trajectory(model, env, H, h)
    traj.q[1] = q0
    traj.q[2] = q1
    # traj.u[1] = control_saturation(policy(p, traj.q[2], traj, 1), uL, uU) #@@@
    # traj.w[1] = disturbances(d, traj.q[2], 1) #@@@

    # initialize interior point solver (for pre-factorization)
    z = zeros(num_var(model, s.env))
    θ = zeros(num_data(model))
    z_initialize!(z, model, env, traj.q[2])
    θ_initialize!(θ, model, traj.q[1], traj.q[2], traj.u[1], traj.w[1], model.μ_world, h)

    ip = eval(ip_type)(z, θ,
        s = space,
        idx_ineq = inequality_indices(model, env),
        idx_soc = soc_indices(model, env),
        r! = r!,
        rm! = rm!,
        rz! = rz!,
        rθ! = rθ!,
        rz = rz,
        rθ = rθ,
        iy1 = linearization_var_index(model, env)[2],
        iy2 = linearization_var_index(model, env)[3],
        ibil = linearization_term_index(model, env)[3],
        opts = ip_opts)

    # pre-allocate for gradients
    traj_deriv = contact_derivative_trajectory(model, env, ip.δz, H)

    Simulator(
        s,
        traj,
        traj_deriv,
        p, uL, uU,
        d,
        ip,
        sim_opts)
end

function step!(sim::Simulator, t)
    # simulation
    model = sim.s.model
    env = sim.s.env

    # unpack
    q = sim.traj.q
    u = sim.traj.u
    w = sim.traj.w
    h = sim.traj.h
    ip = sim.ip
    z = ip.z
    θ = ip.θ

    # t = 1 2 3
    # u1 = traj.u[t]
    # w1 = traj.w[t]
    # γ1 = traj.γ[t]
    # b1 = traj.b[t]
    # q0 = traj.q[t]
    # q1 = traj.q[t+1]
    # q2 = traj.q[t+2]

    # policy
    u[t] = control_saturation(policy(sim.p, q[t+1], sim.traj, t), sim.uL, sim.uU)

    # disturbances
    w[t] = disturbances(sim.d, q[t+1], t)

    # initialize
    if sim.opts.warmstart
        z_warmstart!(z, model, env, q[t+1], sim.opts.z_warmstart, ip.idx_ineq)
        sim.ip.opts.κ_init = sim.opts.κ_warmstart
    else
        z_initialize!(z, model, env, q[t+1])
    end
    θ_initialize!(θ, model, q[t], q[t+1], u[t], w[t], model.μ_world, h)

    # solve
    status = interior_point_solve!(ip)

    if status
        # parse result
        q2, γ, b, _ = unpack_z(model, env, z)
        sim.traj.z[t] = z # TODO: maybe not use copy
        sim.traj.θ[t] = θ
        sim.traj.q[t+2] = q2
        sim.traj.γ[t] = γ
        sim.traj.b[t] = b
        sim.traj.κ[1] = ip.κ[1] # the last κ used in the solve.

        if sim.ip.opts.diff_sol
            sim.deriv_traj.dq2dq0[t] = sim.deriv_traj.vqq
            sim.deriv_traj.dq2dq1[t] = sim.deriv_traj.vqqq
            sim.deriv_traj.dq2du[t] = sim.deriv_traj.vqu
            sim.deriv_traj.dγdq0[t] = sim.deriv_traj.vγq
            sim.deriv_traj.dγdq1[t] = sim.deriv_traj.vγqq
            sim.deriv_traj.dγdu[t] = sim.deriv_traj.vγu
            sim.deriv_traj.dbdq0[t] = sim.deriv_traj.vbq
            sim.deriv_traj.dbdq1[t] = sim.deriv_traj.vbqq
            sim.deriv_traj.dbdu[t] = sim.deriv_traj.vbu
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
    z_initialize!(sim.ip.z, sim.s.model, sim.s.env, sim.traj.q[2])

    status = true

    # simulate
    for t = 1:sim.traj.H
        verbose && println("t = $t / $(sim.traj.H)")
        status = step!(sim, t)
        !status && (@error "failed step (t = $t)")
    end

    return status
end

T = Float64
vis = Visualizer()
open(vis)
include(joinpath(module_dir(), "src", "dynamics", "quadruped", "visuals.jl"))

const ContactControl = Main
s = get_simulation("quadruped", "flat_2D_lc", "flat")

ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
    load_type = :split_traj_alt))

H = 10
H = ref_traj.H
h = ref_traj.h
q0_sim = SVector{s.model.dim.q}(deepcopy(ref_traj.q[1]))
q1_sim = SVector{s.model.dim.q}(deepcopy(ref_traj.q[2]))

p = open_loop_policy(ref_traj.u; N_sample = 1)

sim_b = simulator(s, q0_sim, q1_sim, h, H, p=p, ip_type=:interior_point,
    sim_opts = SimulatorOptions(warmstart=true))
sim_m = simulator(s, q0_sim, q1_sim, h, H, p=p, ip_type=:mehrotra,
    ip_opts = MehrotraOptions(max_iter_inner=100),
    sim_opts = SimulatorOptions(warmstart=false))

@time simulate!(sim_b)
@time simulate!(sim_m)
# @btime simulate!(deepcopy(sim_b))
# @btime simulate!(deepcopy(sim_m))

plot_surface!(vis, s.env, ylims=[-0.5, 0.5])
anim = visualize_meshrobot!(vis, s.model, sim_b.traj, α=0.3, sample=1, name=:basic)
anim = visualize_meshrobot!(vis, s.model, sim_m.traj, α=1.0, sample=1, name=:mehrotra, anim=anim)


plot(hcat(sim_b.traj.q...)', legend=false)
plot!(hcat(sim_m.traj.q...)', legend=false)
sim_m.traj.q[3]

# filename = "quadruped_mehrotra_vs_basic"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
