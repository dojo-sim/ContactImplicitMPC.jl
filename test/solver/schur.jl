@testset "Schur Complement" begin
    T = Float64
    n = 11
    m = 16
    M = rand(n+m,n+m)
    S = ContactControl.Schur(M, n=n, m=m)

    u = rand(T,n)
    v = rand(T,m)
    D = rand(m,m)
    ContactControl.schur_factorize!(S, D)
    # @benchmark ContactControl.schur_factorize!(S, D)
    ContactControl.schur_solve!(S, u, v)
    # @benchmark ContactControl.schur_solve!(S, u, v)
    M1 = [S.A S.B; S.C D]
    @test norm(M1*[S.x; S.y] - [u; v], Inf) < 1e-10
end


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

r0 = zeros(nz)
rz0 = zeros(nz, nz)
rθ0 = zeros(nz, nθ)
z0 = max.(1e-10, ref_traj.z[10])
θ0 = ref_traj.θ[10]
κ0 = 0.0
s.res.r!(r0, z0, θ0, κ0)
s.res.rz!(rz0, z0, θ0)
s.res.rθ!(rθ0, z0, θ0)

A = rz[idyn, ix]
B = rz[idyn, iy1]
C = rz[irst, ix]
D = rz[irst, iy1] - rz[irst, iy2] * Diagonal(z0[iy2] ./ z0[iy1])
M = [A B ; C D]
S = ContactControl.Schur(M, n=nx, m=ny)


u = r0[idyn]
v = r0[irst] - r0[ibil] ./ z0[iy1]
ContactControl.schur_factorize!(S, D)
ContactControl.schur_solve!(S, u, v)
M1 = [S.A S.B; S.C D]
norm(M1*[S.x; S.y] - [u; v], Inf)
@test norm(M1*[S.x; S.y] - [u; v], Inf) < 1e-10
