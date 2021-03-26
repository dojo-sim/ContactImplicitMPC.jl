include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
include("mpc.jl")
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
model = get_model("hopper_2D")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
ref_traj = get_trajectory("hopper_2D", "gait_forward", load_type=:joint_traj)
H = ref_traj.H
h = ref_traj.h
κ = 1.0e-4

function dummy_mpc(model::ContactDynamicsModel, core::Newton, mpc::MPC13)
    # Unpack
    impl = mpc.impl
    ref_traj = mpc.ref_traj
    m_opts = mpc.m_opts
    cost = core.cost
    opts = core.opts
    N_sample = mpc.m_opts.N_sample
    h_sim = mpc.ref_traj.h/N_sample

    # Initialization
    q0 = deepcopy(ref_traj.q[1])
    q1 = deepcopy(ref_traj.q[2])
    q0_sim = SVector{model.dim.q}(deepcopy(q1 + (q0 - q1)/N_sample))
    q1_sim = SVector{model.dim.q}(deepcopy(q1))
    push!(mpc.q_sim, [q0_sim, q1_sim]...)

    for l = 1:m_opts.M
        impl = ImplicitTraj(ref_traj, model, κ=mpc.m_opts.κ)
        warm_start = l > 1
        newton_solve!(core, model, impl, ref_traj; verbose=true, warm_start=warm_start, q0=q0, q1=q1)
        # Apply control and rollout dynamics
        # u_zoh  = SVector{nu}.([deepcopy(ref_traj.u[1])/N_sample for j=1:N_sample])
        u_zoh  = SVector{nu}.([deepcopy(core.traj.u[1])/N_sample for j=1:N_sample])
        w_amp = 0.05
        w_zoh  = SVector{nw}.([w_amp*[-rand(1); zeros(nw-1)]/N_sample for j=1:N_sample])
        sim = simulator(model, SVector{nq}(deepcopy(mpc.q_sim[end-1])), SVector{nq}(deepcopy(mpc.q_sim[end])), h_sim, N_sample,
            p = open_loop_policy(u_zoh, h_sim), #TODO works better need to investigate
            d = open_loop_disturbances(w_zoh, h_sim), #TODO works better need to investigate
            r! = model.res.r!,
            rz! = model.res.rz!,
            rθ! = model.res.rθ!,
            rz = model.spa.rz_sp,
            rθ = model.spa.rθ_sp,
            ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 2e-8, κ_init = 1e-8),
            sim_opts = SimulatorOptions(warmstart = true)
            )
        simulate!(sim; verbose = false)
        println("Δq: ", scn.(norm(ref_traj.q[3] - sim.traj.q[end]), digits=2))
        ########################
        # plt = plot(layout=grid(2,1,heights=[0.7, 0.3], figsize=[(1000, 1000),(400,400)]), legend=false, xlims=(0,20))
        # plot!(plt[1,1], hcat(Vector.(core.traj.q)...)', color=:blue, linewidth=1.0)
        # plot!(plt[1,1], hcat(Vector.(ref_traj.q)...)', linestyle=:dot, color=:red, linewidth=3.0)
        #
        # scatter!((2-1/N_sample)ones(model.dim.q), sim.traj.q[1], markersize=8.0, color=:lightgreen)
        # scatter!(2*ones(model.dim.q), sim.traj.q[2], markersize=8.0, color=:lightgreen)
        # scatter!(plt[1,1], 3*ones(model.dim.q), sim.traj.q[end], markersize=8.0, color=:lightgreen)
        #
        # scatter!(plt[1,1], 1*ones(model.dim.q), q0, markersize=6.0, color=:blue)
        # scatter!(plt[1,1], 2*ones(model.dim.q), q1, markersize=6.0, color=:blue)
        #
        # scatter!(plt[1,1], 1*ones(model.dim.q), ref_traj.q[1], markersize=4.0, color=:red)
        # scatter!(plt[1,1], 2*ones(model.dim.q), ref_traj.q[2], markersize=4.0, color=:red)
        # scatter!(plt[1,1], 3*ones(model.dim.q), ref_traj.q[3], markersize=4.0, color=:red)
        #
        # plot!(plt[2,1], hcat(Vector.(core.traj.u)...)', color=:blue, linewidth=1.0)
        # plot!(plt[2,1], hcat(Vector.(ref_traj.u)...)', linestyle=:dot, color=:red, linewidth=3.0)
        # display(plt)
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
    return nothing
end

ref_traj0 = deepcopy(ref_traj)
n_opts0 = NewtonOptions(r_tol=3e-4, κ_init=κ, κ_tol=2κ, solver_inner_iter=5)
m_opts0 = MPC13Options{T}(N_sample=2, M=270, H_mpc=8, κ=κ)
cost0 = CostFunction(H, model.dim,
    q = [Diagonal(1.0e-2 * [1,1,1,1])   for t = 1:H],
    u = [Diagonal(1.0e-0 * [1e-3, 1e1]) for t = 1:H],
    γ = [Diagonal(1.0e-100 * ones(nc)) for t = 1:H],
    b = [Diagonal(1.0e-100 * ones(nb)) for t = 1:H])
core0 = Newton(H, h, model, cost=cost0, opts=n_opts0)
mpc0 = MPC13(model, ref_traj0, m_opts=m_opts0)
@time dummy_mpc(model, core0, mpc0)

visualize!(vis, model, mpc0.q_sim[1:1:end], Δt=h/m_opts0.N_sample, name=:mpc)
plt = plot(layout=(2,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], m_opts0.N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.(mpc0.q_sim)...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:2], m_opts0.N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:2] for u in mpc0.u_sim]*m_opts0.N_sample)...)', color=:blue, linewidth=1.0)



# filename = "hopper_open_loop"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
