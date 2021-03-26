# MPC options
@with_kw mutable struct MPCOptions{T}
    N_sample::Int=3           # Hypersampling factor for dynamics simulation
    M::Int=3                  # Number of MPC loops
    H_mpc::Int=10             # Horizon of the MPC solver (H_mpc < H)
    κ::T=1e-3                 # Central path parameter used in the dynamics
    open_loop_mpc::Bool=false # Execute the reference trajectory open-loop instead of doing real MPC (~useful to assert the disturbance impact)
    w_amp::T=0.05             # Amplitude of the disturbance
end

mutable struct MPC{T,nq,nu,nw,nc,nb}
    q0_con::SizedArray{Tuple{nq},T,1,1,Vector{T}} # initial state for the controller
    q1_con::SizedArray{Tuple{nq},T,1,1,Vector{T}} # initial state for the controller
    q_sim::Vector{SizedArray{Tuple{nq},T,1,1,Vector{T}}} # history of simulator configurations
    u_sim::Vector{SizedArray{Tuple{nu},T,1,1,Vector{T}}} # history of simulator controls
    w_sim::Vector{SizedArray{Tuple{nw},T,1,1,Vector{T}}} # history of simulator disturbances
    γ_sim::Vector{SizedArray{Tuple{nc},T,1,1,Vector{T}}} # history of simulator impact forces
    b_sim::Vector{SizedArray{Tuple{nb},T,1,1,Vector{T}}} # history of simulator friction forces
    m_opts::MPCOptions
    ref_traj::ContactTraj
    impl::ImplicitTraj
    q_stride::SizedArray{Tuple{nq},T,1,1,Vector{T}} # change in q required to loop the ref trajectories
end

function MPC(model::ContactDynamicsModel, ref_traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ};
        m_opts::MPCOptions=MPCOptions()) where {T,nq,nu,nw,nc,nb,nz,nθ}
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
    return MPC{T,nq,nu,nw,nc,nb}(q0_con, q1_con, q_sim, u_sim, w_sim, γ_sim, b_sim, m_opts, ref_traj, impl, q_stride)
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

function get_stride(model::Quadruped, traj::ContactTraj)
    q_stride = zeros(SizedVector{nq})
    q_stride[1] = traj.q[end-1][1] - traj.q[1][1]
    return q_stride
end

function get_stride(model::Hopper2D, traj::ContactTraj)
    q_stride = zeros(SizedVector{nq})
    q_stride[1] = traj.q[end-1][1] - traj.q[1][1]
    return q_stride
end




# T = Float64
# κ = 1e-4
# # model = get_model("quadruped")
# ref_traj0 = get_trajectory("quadruped", "gait1")
# # time
# h = ref_traj0.h
# H = ref_traj0.H
#
# nq = model.dim.q
# nu = model.dim.u
# nc = model.dim.c
# nb = model.dim.b
#
# # initial conditions
# q0 = SVector{model.dim.q}(ref_traj0.q[1])
# q1 = SVector{model.dim.q}(ref_traj0.q[2])
#
#
# sim0 = simulator(model, q0, q1, h, H,
#     p = open_loop_policy(SVector{model.dim.u}.(ref_traj0.u), h),
#     r! = model.res.r,
#     rz! = model.res.rz,
#     rθ! = model.res.rθ,
#     rz = model.spa.rz_sp,
#     rθ = model.spa.rθ_sp,
#     ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 2κ, κ_init = κ),
#     sim_opts = SimulatorOptions(warmstart = true)
#     )
#
# simulate!(sim0; verbose = false)
# ref_traj0 = deepcopy(sim0.traj)
# traj0 = deepcopy(sim0.traj)

function dummy_mpc(model::ContactDynamicsModel, core::Newton, mpc::MPC)
    impl = mpc.impl
    ref_traj = mpc.ref_traj
    m_opts = mpc.m_opts
    cost = core.cost
    opts = core.opts

    q0 = deepcopy(ref_traj.q[1])
    q1 = deepcopy(ref_traj.q[2])
    for l = 1:m_opts.M
        println("mpc:", l)
        impl = ImplicitTraj(H, model)
        linearization!(model, ref_traj, impl, ref_traj.κ)
        warm_start = l > 1
        newton_solve!(model, core, impl, ref_traj; warm_start=warm_start, q0=q0, q1=q1)

        # Apply control and rollout dynamics
        N_sample = mpc.m_opts.N_sample
        h_sim = mpc.ref_traj.h/N_sample
        q0_sim = SVector{model.dim.q}(deepcopy(core.traj.q[2] + (core.traj.q[1] - core.traj.q[2])/N_sample ))
        q1_sim = SVector{model.dim.q}(deepcopy(core.traj.q[2]))
        sim = simulator(model, q0_sim, q1_sim, h_sim, N_sample,
            p = open_loop_policy([SVector{model.dim.u}(deepcopy(core.traj.u[1])) for i=1:N_sample], h_sim),
            r! = model.res.r,
            rz! = model.res.rz,
            rθ! = model.res.rθ,
            rz = model.spa.rz_sp,
            rθ = model.spa.rθ_sp,
            ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 2κ, κ_init = κ),
            sim_opts = SimulatorOptions(warmstart = true)
            )
        simulate!(sim; verbose = false)
        ########################
        plt = plot(legend=false, xlims=(0,15))
        plot!(hcat(Vector.(core.traj.q)...)', color=:blue, linewidth=1.0)
        plot!(hcat(Vector.(ref_traj.q)...)', linestyle=:dot, color=:red, linewidth=3.0)
        scatter!((2-1/N_sample)ones(model.dim.q), q0_sim, markersize=8.0, color=:lightgreen)
        scatter!(2*ones(model.dim.q), q1_sim, markersize=8.0, color=:lightgreen)
        @show length(sim.traj.q)
        scatter!(3*ones(model.dim.q), sim.traj.q[N_sample+2], markersize=8.0, color=:lightgreen)
        scatter!(1*ones(model.dim.q), q0, markersize=6.0, color=:blue)
        scatter!(2*ones(model.dim.q), q1, markersize=6.0, color=:blue)
        scatter!(1*ones(model.dim.q), ref_traj.q[1], markersize=4.0, color=:red)
        scatter!(2*ones(model.dim.q), ref_traj.q[2], markersize=4.0, color=:red)
        scatter!(3*ones(model.dim.q), ref_traj.q[3], markersize=4.0, color=:red)
        display(plt)
        @show norm(q1 - ref_traj.q[3])
        ########################
        q0 = deepcopy(q1)
        q1 = deepcopy(sim.traj.q[N_sample+2])
        rot_n_stride!(ref_traj, mpc.q_stride)
    end
    return nothing
end
