@testset "Test MPC" begin
    T = Float64
    # get hopper model
    # model = ContactControl.get_model("hopper_2D")
    model = ContactControl.get_model("quadruped")

    # get trajectory
    # ref_traj = ContactControl.get_trajectory("hopper_2D", "gait_in_place", load_type=:joint_traj)
    # ref_traj = ContactControl.get_trajectory("hopper_2D", "gait_forward", load_type=:joint_traj)
    ref_traj = ContactControl.get_trajectory("quadruped", "gait1", load_type=:split_traj)
    H = ref_traj.H
    h = ref_traj.h
    κ = 1.0e-4

    n_opts = ContactControl.NewtonOptions(r_tol=3e-4, κ_init=κ, κ_tol=2κ, solver_inner_iter=5)
    m_opts = ContactControl.MPCOptions{T}(
                N_sample=10,
                M=2*H,
                H_mpc=10,
                κ=κ,
                κ_sim=1e-8,
                r_tol_sim=1e-8,
                open_loop_mpc=false,
                w_amp=[0.05, 0.00],
                live_plotting=false)
    # cost = ContactControl.CostFunction(H, model.dim,
    #     q = [Diagonal(1.0e-1 * [1,1,1,1])   for t = 1:m_opts.H_mpc],
    #     u = [Diagonal(1.0e-0 * [1e-3, 1e1]) for t = 1:m_opts.H_mpc],
    #     γ = [Diagonal(1.0e-100 * ones(nc)) for t = 1:m_opts.H_mpc],
    #     b = [Diagonal(1.0e-100 * ones(nb)) for t = 1:m_opts.H_mpc])
    cost = ContactControl.CostFunction(H, model.dim,
        q = [Diagonal(1e-2 * [0.02,0.02,1,.15,.15,.15,.15,.15,.15,.15,.15,]) for t = 1:m_opts.H_mpc],
        u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:m_opts.H_mpc],
        γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:m_opts.H_mpc],
        b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:m_opts.H_mpc])
    core = ContactControl.Newton(m_opts.H_mpc, h, model, cost=cost, opts=n_opts)
    mpc = ContactControl.MPC(model, ref_traj, m_opts=m_opts)
    ContactControl.dummy_mpc(model, core, mpc)

    # Check length
    N_sample = mpc.m_opts.N_sample
    M = mpc.m_opts.M
    H_mpc = mpc.m_opts.H_mpc

    @test length(mpc.q_sim) == N_sample*M+2
    @test length(mpc.u_sim) == N_sample*M
    @test length(mpc.w_sim) == N_sample*M
    @test length(mpc.γ_sim) == N_sample*M
    @test length(mpc.b_sim) == N_sample*M

    @test core.traj.H == H_mpc
    @test length(core.traj.q) == H_mpc+2
    @test length(core.traj.u) == H_mpc
    @test core.traj_cand.H == H_mpc
    @test length(core.traj.q) == H_mpc+2
    @test length(core.traj.u) == H_mpc
    @test length(core.ν) == H_mpc
    @test length(core.ν_cand) == H_mpc

    # Check tracking performance

    function tracking_error(ref_traj::ContactControl.ContactTraj, mpc::ContactControl.MPC)
        N_sample = mpc.m_opts.N_sample
        M = mpc.m_opts.M
        q_error = []
        u_error = []
        γ_error = []
        b_error = []
        q_sim = mpc.q_sim[3:N_sample:end]
        u_sim = mpc.u_sim[1:N_sample:end]
        γ_sim = mpc.γ_sim[1:N_sample:end]
        b_sim = mpc.b_sim[1:N_sample:end]
        # plt = plot(legend=false)
        for t = 1:M
            push!(q_error, norm(ref_traj.q[3+(t-1)%H][2:end] - q_sim[t][2:end]))
            push!(u_error, norm(ref_traj.u[1+(t-1)%H] - u_sim[t]*N_sample))
            push!(γ_error, norm(ref_traj.γ[1+(t-1)%H] - γ_sim[t]*N_sample))
            push!(b_error, norm(ref_traj.b[1+(t-1)%H] - b_sim[t]*N_sample))
            # scatter!([t,t], ref_traj.q[3+(t-1)%H][2:3], color=:red)
            # scatter!([t,t], q_sim[t][2:3], color=:blue)
            # scatter!([t,t], ref_traj.u[1+(t-1)%H][1:2]/N_sample, color=:red)
            # scatter!([t,t], u_sim[t][1:2], color=:blue)
        end
        # display(plt)
        return mean(q_error), mean(u_error), mean(γ_error), mean(b_error)
    end

    # Check the tracking error with disturbances
    q_error, u_error, γ_error, b_error = tracking_error(ref_traj, mpc)
    # @test q_error < 0.05
    # @test u_error < 0.10
    # @test γ_error < 0.05
    # @test b_error < 0.05
    @test q_error < 0.10
    @test u_error < 0.10
    @test γ_error < 0.19
    @test b_error < 0.15


    # get trajectory
    # ref_traj = ContactControl.get_trajectory("hopper_2D", "gait_forward", load_type=:joint_traj)
    ref_traj = ContactControl.get_trajectory("quadruped", "gait1", load_type=:split_traj)
    H = ref_traj.H
    h = ref_traj.h
    κ = 1.0e-4

    n_opts = ContactControl.NewtonOptions(r_tol=3e-4, κ_init=κ, κ_tol=2κ, solver_inner_iter=5)
    m_opts = ContactControl.MPCOptions{T}(
                N_sample=10,
                M=2*H,
                H_mpc=10,
                κ=κ,
                κ_sim=1e-8,
                r_tol_sim=1e-8,
                open_loop_mpc=false,
                w_amp=zeros(model.dim.w),
                live_plotting=false)
    # cost = ContactControl.CostFunction(H, model.dim,
    #     q = [Diagonal(1.0e-1 * [1,1,1,1])   for t = 1:m_opts.H_mpc],
    #     u = [Diagonal(1.0e-0 * [1e-3, 1e1]) for t = 1:m_opts.H_mpc],
    #     γ = [Diagonal(1.0e-100 * ones(nc)) for t = 1:m_opts.H_mpc],
    #     b = [Diagonal(1.0e-100 * ones(nb)) for t = 1:m_opts.H_mpc])
    cost = ContactControl.CostFunction(H, model.dim,
        q = [Diagonal(1e-2 * [0.02,0.02,1,.15,.15,.15,.15,.15,.15,.15,.15,]) for t = 1:m_opts.H_mpc],
        u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:m_opts.H_mpc],
        γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:m_opts.H_mpc],
        b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:m_opts.H_mpc])
    core = ContactControl.Newton(m_opts.H_mpc, h, model, cost=cost, opts=n_opts)
    mpc = ContactControl.MPC(model, ref_traj, m_opts=m_opts)
    ContactControl.dummy_mpc(model, core, mpc)

    # Check the tracking error with no disturbances
    q_error, u_error, γ_error, b_error = tracking_error(ref_traj, mpc)
    # @test q_error < 0.05
    # @test u_error < 0.01
    # @test γ_error < 0.05
    # @test b_error < 0.02
    @test q_error < 0.10
    @test u_error < 0.01
    @test γ_error < 0.15
    @test b_error < 0.10
end

# T = Float64
# # get hopper model
# # model = ContactControl.get_model("hopper_2D")
# model = ContactControl.get_model("quadruped")
#
# # get trajectory
# # ref_traj = ContactControl.get_trajectory("hopper_2D", "gait_in_place", load_type=:joint_traj)
# # ref_traj = ContactControl.get_trajectory("hopper_2D", "gait_forward", load_type=:joint_traj)
# ref_traj = ContactControl.get_trajectory("quadruped", "gait1", load_type=:split_traj)
# H = ref_traj.H
# h = ref_traj.h
# κ = 1.0e-4
#
# n_opts = ContactControl.NewtonOptions(r_tol=3e-4, κ_init=κ, κ_tol=2κ, solver_inner_iter=5)
# m_opts = ContactControl.MPCOptions{T}(
#             N_sample=10,
#             M=2*H,
#             H_mpc=10,
#             κ=κ,
#             κ_sim=1e-8,
#             r_tol_sim=1e-8,
#             open_loop_mpc=false,
#             w_amp=[0.05, 0.00],
#             live_plotting=false)
# # cost = ContactControl.CostFunction(H, model.dim,
# #     q = [Diagonal(1.0e-1 * [1,1,1,1])   for t = 1:m_opts.H_mpc],
# #     u = [Diagonal(1.0e-0 * [1e-3, 1e1]) for t = 1:m_opts.H_mpc],
# #     γ = [Diagonal(1.0e-100 * ones(nc)) for t = 1:m_opts.H_mpc],
# #     b = [Diagonal(1.0e-100 * ones(nb)) for t = 1:m_opts.H_mpc])
# cost = CostFunction(H, model.dim,
#     q = [Diagonal(1e-2 * [0.02,0.02,1,.15,.15,.15,.15,.15,.15,.15,.15,]) for t = 1:m_opts.H_mpc],
#     u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:m_opts.H_mpc],
#     γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:m_opts.H_mpc],
#     b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:m_opts.H_mpc])
# core = ContactControl.Newton(m_opts.H_mpc, h, model, cost=cost, opts=n_opts)
# mpc = ContactControl.MPC(model, ref_traj, m_opts=m_opts)
# get_stride(model, ref_traj)
# ContactControl.dummy_mpc(model, core, mpc)
