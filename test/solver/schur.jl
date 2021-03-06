@testset "Schur Complement" begin
    T = Float64
    n = 11
    m = 16
    M = rand(n+m,n+m)
    S = ContactImplicitMPC.Schur(M, n=n, m=m)

    u = rand(T,n)
    v = rand(T,m)
    D = rand(m,m)
    ContactImplicitMPC.schur_factorize!(S, D)
    # @benchmark ContactImplicitMPC.schur_factorize!(S, D)
    ContactImplicitMPC.schur_solve!(S, u, v)
    # @benchmark ContactImplicitMPC.schur_solve!(S, u, v)
    M1 = [S.A S.B; S.C D]
    @test norm(M1*[S.x; S.y] - [u; v], Inf) < 1e-9
end

@testset "Schur Complement: ill-conditioning" begin
    T = Float64
    s = get_simulation("flamingo", "flat_2D_lc", "flat")
    ref_traj = deepcopy(get_trajectory(s.model, s.env,
        joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
        load_type = :split_traj_alt))

    ix, iy1, iy2 = linearization_var_index(s.model, s.env)
    idyn, irst, ibil = linearization_term_index(s.model, s.env)
    nz = num_var(s.model, s.env)
    nθ = num_data(s.model)
    nx = length(ix)
    ny = length(iy1)

    r0 = zeros(T, nz)
    rz0 = zeros(T, nz, nz)
    rθ0 = zeros(T, nz, nθ)
    z0 = zeros(T, nz)
    θ0 = zeros(T, nθ)
    z0 .= max.(1e-6, ref_traj.z[10])
    θ0 .= ref_traj.θ[10]
    κ0 = 0.0
    s.res.r!(r0, z0, θ0, κ0)
    s.res.rz!(rz0, z0, θ0)
    s.res.rθ!(rθ0, z0, θ0)

    A = rz0[idyn, ix]
    B = rz0[idyn, iy1]
    C = rz0[irst, ix]
    D = rz0[irst, iy1] - rz0[irst, iy2] * Diagonal(z0[iy2] ./ z0[iy1])
    M = [A B ; C D]
    S = ContactImplicitMPC.Schur(M, n=nx, m=ny)

    u = r0[idyn]
    v = r0[irst] - r0[ibil] ./ z0[iy1]

    # @benchmark schur_factorize!(S, D)
    ContactImplicitMPC.schur_factorize!(S, D)
    # @benchmark ContactImplicitMPC.schur_solve!(S, u, v)
    ContactImplicitMPC.schur_solve!(S, u, v)
    M1 = [S.A S.B; S.C D]
    norm(M1*[S.x; S.y] - [u; v], Inf)
    @test norm(M1*[S.x; S.y] - [u; v], Inf) < 1e-7
end
