@testset "Trajectory" begin

    s = get_simulation("flamingo", "flat_2D_lc", "flat")
    ref_traj = deepcopy(ContactImplicitMPC.get_trajectory(s.model, s.env,
        joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
        load_type = :split_traj_alt))

    N_repeat = 5
    N_sample = 1
    dupl_traj = repeat_ref_traj(ref_traj, N_repeat, idx_shift=[1])
    @test all(ContactImplicitMPC.tracking_error(ref_traj, dupl_traj, N_sample, idx_shift=[1]) .== 0.0)

end
