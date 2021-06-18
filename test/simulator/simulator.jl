@testset "Simulator" begin

    T = Float64
    # vis = Visualizer()
    # open(vis)
    # include(joinpath(module_dir(), "src", "dynamics", "quadruped", "visuals.jl"))

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
    qerr, uerr, γerr, berr = tracking_error(ref_traj, sim_b.traj, 1, idx_shift = 1)
    @test qerr < 0.0555 * 1.5
    @test uerr < 0.0001 * 1.5
    @test γerr < 0.0491 * 1.5
    @test berr < 0.0871 * 1.5
    qerr > 0.0555 * 1.2 && @warn "simulator regression: q"
    uerr > 0.0001 * 1.2 && @warn "simulator regression: u"
    γerr > 0.0491 * 1.2 && @warn "simulator regression: γ"
    berr > 0.0871 * 1.2 && @warn "simulator regression: b"

    qerr, uerr, γerr, berr = tracking_error(ref_traj, sim_m.traj, 1, idx_shift = 1)
    @test qerr < 0.0170 * 1.5
    @test uerr < 0.0001 * 1.5
    @test γerr < 0.0138 * 1.5
    @test berr < 0.0660 * 1.5
    qerr > 0.0170 * 1.2 && @warn "simulator regression: q"
    uerr > 0.0001 * 1.2 && @warn "simulator regression: u"
    γerr > 0.0138 * 1.2 && @warn "simulator regression: γ"
    berr > 0.0660 * 1.2 && @warn "simulator regression: b"


    # @btime simulate!(deepcopy(sim_b))
    # @btime simulate!(deepcopy(sim_m))

    # plot_surface!(vis, s.env, ylims=[-0.5, 0.5])
    # anim = visualize_meshrobot!(vis, s.model, sim_b.traj, α=0.3, sample=1, name=:basic)
    # anim = visualize_meshrobot!(vis, s.model, sim_m.traj, α=1.0, sample=1, name=:mehrotra, anim=anim)

    # plot(hcat(sim_b.traj.q...)', legend=false)
    # plot!(hcat(sim_m.traj.q...)', legend=false)

end



#
# T = Float64
# # vis = Visualizer()
# # open(vis)
# include(joinpath(module_dir(), "src", "dynamics", "quadruped", "visuals.jl"))
#
# # s = get_simulation("quadruped", "flat_2D_lc", "flat")
#
# ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
#     joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
#     load_type = :split_traj_alt))
#
# im_traj = ImplicitTraj(ref_traj, s, ip_type=:mehrotra)
#
# H = 10
# H = ref_traj.H
# h = ref_traj.h
# q0_sim = SVector{s.model.dim.q}(deepcopy(ref_traj.q[1]))
# q1_sim = SVector{s.model.dim.q}(deepcopy(ref_traj.q[2]))
#
# p = open_loop_policy(ref_traj.u; N_sample = 1)
#
# sim_b = simulator(s, q0_sim, q1_sim, h, H, p=p, ip_type=:interior_point,
#     sim_opts = SimulatorOptions(warmstart=true))
# sim_m = simulator(s, q0_sim, q1_sim, h, H, p=p, ip_type=:mehrotra,
#     ip_opts = MehrotraOptions(max_iter_inner=100),
#     sim_opts = SimulatorOptions(warmstart=false))
