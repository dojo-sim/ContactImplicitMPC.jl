@testset "Open-loop disturbance" begin

    nw = 4
    nx = 6
    H = 4
    N_sample = 3
    H_sim = N_sample*H
    w = [rand(nw) for t=1:H]
    d = open_loop_disturbances(w, N_sample)

    idx = [1,1,1,2,2,2,3,3,3,4,4,4]
    cnt = [1,2,3,1,2,3,1,2,3,1,2,3]
    for t = 1:H_sim
        x = rand(nx)
        wt = disturbances(d, x, t)
        @test d.idx == idx[t]
        @test d.cnt == cnt[t]
        @test norm(wt - d.w[d.idx] / d.N_sample) < 1e-10
    end

end

@testset "Open-loop policy" begin

    nu = 4
    nx = 6
    H = 4
    N_sample = 3
    H_sim = N_sample*H
    u = [rand(nu) for t=1:H]
    p = open_loop_policy(u, N_sample=N_sample)

    idx = [1,1,1,2,2,2,3,3,3,4,4,4]
    cnt = [1,2,3,1,2,3,1,2,3,1,2,3]
    for t = 1:H_sim
        x = nothing
        traj = nothing
        ut = policy(p, x, traj, t)
        @test p.idx == idx[t]
        @test p.cnt == cnt[t]
        @test norm(ut - p.u[p.idx] / p.N_sample) < 1e-10
    end

end
