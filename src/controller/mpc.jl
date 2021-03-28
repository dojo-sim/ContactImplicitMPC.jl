# MPC11 options
@with_kw mutable struct MPC11Options11{T}
    N_sample::Int=3           # Hypersampling factor for dynamics simulation
    M::Int=3                  # Number of MPC11 loops
    H_mpc::Int=10             # Horizon of the MPC11 solver (H_mpc < H)
    κ::T=1e-3                 # Central path parameter used in the implicit dynamics
    κ_sim::T=1e-8             # Central path parameter used in the simulation dynamics
    r_tol_sim::T=1e-8         # Residual tolerance for the dynamics simulation
    open_loop_mpc::Bool=false # Execute the reference trajectory open-loop instead of doing real MPC11 (~useful to assert the disturbance impact)
    w_amp::Vector{T}=[0.05]   # Amplitude of the disturbance
    live_plotting::Bool=false # Use the live plotting tool to debug
end

mutable struct MPC11{T,nq,nu,nw,nc,nb}
    q0_con::SizedArray{Tuple{nq},T,1,1,Vector{T}} # initial state for the controller
    q1_con::SizedArray{Tuple{nq},T,1,1,Vector{T}} # initial state for the controller
    q_sim::Vector{SizedArray{Tuple{nq},T,1,1,Vector{T}}} # history of simulator configurations
    u_sim::Vector{SizedArray{Tuple{nu},T,1,1,Vector{T}}} # history of simulator controls
    w_sim::Vector{SizedArray{Tuple{nw},T,1,1,Vector{T}}} # history of simulator disturbances
    γ_sim::Vector{SizedArray{Tuple{nc},T,1,1,Vector{T}}} # history of simulator impact forces
    b_sim::Vector{SizedArray{Tuple{nb},T,1,1,Vector{T}}} # history of simulator friction forces
    m_opts::MPC11Options11
    ref_traj::ContactTraj
    impl::ImplicitTraj
    q_stride::SizedArray{Tuple{nq},T,1,1,Vector{T}} # change in q required to loop the ref trajectories
end

function MPC11(model::ContactDynamicsModel, ref_traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ};
        m_opts::MPC11Options11=MPC11Options11()) where {T,nq,nu,nw,nc,nb,nz,nθ}
    q0_con = zeros(SizedArray{Tuple{nq},T,1,1,Vector{T}})
    q1_con = zeros(SizedArray{Tuple{nq},T,1,1,Vector{T}})
    q_sim = Vector{SizedArray{Tuple{nq},T,1,1,Vector{T}}}([])
    u_sim = Vector{SizedArray{Tuple{nu},T,1,1,Vector{T}}}([])
    w_sim = Vector{SizedArray{Tuple{nw},T,1,1,Vector{T}}}([])
    γ_sim = Vector{SizedArray{Tuple{nc},T,1,1,Vector{T}}}([])
    b_sim = Vector{SizedArray{Tuple{nb},T,1,1,Vector{T}}}([])
    H = ref_traj.H
    impl = ImplicitTraj(ref_traj, model, κ=m_opts.κ)
    q_stride = get_stride(model, ref_traj)
    return MPC11{T,nq,nu,nw,nc,nb}(q0_con, q1_con, q_sim, u_sim, w_sim, γ_sim, b_sim, m_opts, ref_traj, impl, q_stride)
end

function rot_n_stride!(traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
    q_stride::SizedArray{Tuple{nq},T,1,1,Vector{T}}) where {T,nq,nu,nw,nc,nb,nz,nθ}
    rotate!(traj)
    mpc_stride!(traj, q_stride)
    return nothing
end

function rotate!(traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}) where {T,nq,nu,nw,nc,nb,nz,nθ}
    q1 = copy(traj.q[1])
    u1 = copy(traj.u[1])
    w1 = copy(traj.w[1])
    γ1 = copy(traj.γ[1])
    b1 = copy(traj.b[1])
    z1 = copy(traj.z[1])
    θ1 = copy(traj.θ[1])
    for t = 1:H+1
        traj.q[t] .= traj.q[t+1]
    end
    for t = 1:H-1
        traj.u[t] .= traj.u[t+1]
        traj.w[t] .= traj.w[t+1]
        traj.γ[t] .= traj.γ[t+1]
        traj.b[t] .= traj.b[t+1]
        traj.z[t] .= traj.z[t+1]
        traj.θ[t] .= traj.θ[t+1]
    end
    traj.q[end] .= q1
    traj.u[end] .= u1
    traj.w[end] .= w1
    traj.γ[end] .= γ1
    traj.b[end] .= b1
    traj.z[end] .= z1
    traj.θ[end] .= θ1
    return nothing
end

"""
    Update the last two cofiguaraton to be equal to the fisrt two up to an constant offset.
"""
function mpc_stride!(traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
    q_stride::SizedArray{Tuple{nq},T,1,1,Vector{T}}) where {T,nq,nu,nw,nc,nb,nz,nθ}
    H = traj.H
    for t = H+1:H+2
        traj.q[t] .= deepcopy(traj.q[t-H]+q_stride)
        update_z!(traj, t-2)
        update_θ!(traj, t-2)
        update_θ!(traj, t-3)
    end
    return nothing
end

function get_stride(model::Hopper2D, traj::ContactTraj)
    q_stride = zeros(SizedVector{nq})
    q_stride[1] = traj.q[end-1][1] - traj.q[1][1]
    return q_stride
end

function get_stride(model::Quadruped, traj::ContactTraj)
    q_stride = zeros(SizedVector{nq})
    q_stride[1] = traj.q[end-1][1] - traj.q[1][1]
    return q_stride
end

function get_stride(model::Biped, traj::ContactTraj)
    q_stride = zeros(SizedVector{nq})
    q_stride[1] = traj.q[end-1][1] - traj.q[1][1]
    return q_stride
end

function dummy_mpc(model::ContactDynamicsModel, core::Newton, mpc::MPC11; verbose::Bool=false)
    elap = 0.0
    # Unpack
    impl = mpc.impl
    ref_traj = mpc.ref_traj
    m_opts = mpc.m_opts
    cost = core.cost
    opts = core.opts
    N_sample = mpc.m_opts.N_sample
    h_sim = mpc.ref_traj.h/N_sample
    nq = model.dim.q

    # Initialization
    q0 = deepcopy(ref_traj.q[1])
    q1 = deepcopy(ref_traj.q[2])
    q0_sim = SVector{nq}(deepcopy(q1 + (q0 - q1)/N_sample))
    q1_sim = SVector{nq}(deepcopy(q1))
    push!(mpc.q_sim, [q0_sim, q1_sim]...)

    for l = 1:m_opts.M
        elap += @elapsed impl = ImplicitTraj(ref_traj, model, κ=m_opts.κ)
        # Get control
        if m_opts.open_loop_mpc
            u_zoh  = SVector{nu}.([deepcopy(ref_traj.u[1])/N_sample for j=1:N_sample])
        else
            warm_start = l > 1
            elap += @elapsed newton_solve!(core, model, impl, ref_traj; verbose=verbose, warm_start=warm_start, q0=q0, q1=q1)
            u_zoh  = SVector{nu}.([deepcopy(core.traj.u[1])/N_sample for j=1:N_sample])
        end
        # Get disturbances
        # Apply control and rollout dynamics
        sim = simulator(model, SVector{nq}(deepcopy(mpc.q_sim[end-1])), SVector{nq}(deepcopy(mpc.q_sim[end])), h_sim, N_sample,
            p = open_loop_policy(u_zoh, h_sim),
            d = random_disturbances(model, m_opts.w_amp/N_sample, N_sample, h_sim),
            r! = model.res.r!,
            rz! = model.res.rz!,
            rθ! = model.res.rθ!,
            rz = model.spa.rz_sp,
            rθ = model.spa.rθ_sp,
            ip_opts = InteriorPointOptions(r_tol = m_opts.r_tol_sim, κ_tol = 2m_opts.κ_sim, κ_init = m_opts.κ_sim),
            sim_opts = SimulatorOptions(warmstart = true)
            )
        simulate!(sim)
        ########################
        if m_opts.live_plotting
            plt = plot(layout=grid(2,1,heights=[0.7, 0.3], figsize=[(1000, 1000),(400,400)]), legend=false, xlims=(0,20))
            plot!(plt[1,1], hcat(Vector.(core.traj.q)...)', color=:blue, linewidth=1.0)
            plot!(plt[1,1], hcat(Vector.(ref_traj.q)...)', linestyle=:dot, color=:red, linewidth=3.0)

            scatter!((2-1/N_sample)ones(model.dim.q), sim.traj.q[1], markersize=8.0, color=:lightgreen)
            scatter!(2*ones(model.dim.q), sim.traj.q[2], markersize=8.0, color=:lightgreen)
            scatter!(plt[1,1], 3*ones(model.dim.q), sim.traj.q[end], markersize=8.0, color=:lightgreen)

            scatter!(plt[1,1], 1*ones(model.dim.q), q0, markersize=6.0, color=:blue)
            scatter!(plt[1,1], 2*ones(model.dim.q), q1, markersize=6.0, color=:blue)

            scatter!(plt[1,1], 1*ones(model.dim.q), ref_traj.q[1], markersize=4.0, color=:red)
            scatter!(plt[1,1], 2*ones(model.dim.q), ref_traj.q[2], markersize=4.0, color=:red)
            scatter!(plt[1,1], 3*ones(model.dim.q), ref_traj.q[3], markersize=4.0, color=:red)

            plot!(plt[2,1], hcat(Vector.(core.traj.u)...)', color=:blue, linewidth=1.0)
            plot!(plt[2,1], hcat(Vector.(ref_traj.u)...)', linestyle=:dot, color=:red, linewidth=3.0)
            display(plt)
        end
        ########################
        push!(mpc.q_sim, deepcopy(sim.traj.q[3:end])...)
        push!(mpc.u_sim, deepcopy(sim.traj.u)...)
        push!(mpc.w_sim, deepcopy(sim.traj.w)...)
        push!(mpc.γ_sim, deepcopy(sim.traj.γ)...)
        push!(mpc.b_sim, deepcopy(sim.traj.b)...)
        ########################
        q0 = deepcopy(q1)
        q1 = deepcopy(sim.traj.q[end])
        rot_n_stride!(ref_traj, mpc.q_stride)
    end
    @show elap
    return nothing
end
