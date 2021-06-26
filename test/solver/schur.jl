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

const ContactControl = Main
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
S = Schur11(M, n=nx, m=ny)

u = r0[idyn]
v = r0[irst] - r0[ibil] ./ z0[iy1]

# @benchmark schur_factorize!(S, D)
# @benchmark ContactControl.schur_solve!(S, u, v)
schur_factorize!(S, D)
ContactControl.schur_solve!(S, u, v)
M1 = [S.A S.B; S.C D]
norm(M1*[S.x; S.y] - [u; v], Inf)
@test norm(M1*[S.x; S.y] - [u; v], Inf) < 1e-10




M = [A B ; C D]
S = Schur(M, n=nx, m=ny)

u = r0[idyn]
v = r0[irst] - r0[ibil] ./ z0[iy1]

schur_factorize!(S, D)
ContactControl.schur_solve!(S, u, v)
M1 = [S.A S.B; S.C D]
norm(M1*[S.x; S.y] - [u; v], Inf)
@test norm(M1*[S.x; S.y] - [u; v], Inf) < 1e-10



sca = Vector(ny ./ sum(abs.(S.D - S.CAiB), dims=2)[:,1])
cond(sca .* (S.D - S.CAiB))





# # S.x = Ai*(us + B*(As*(C*(Ai*us) - vs)))
# # S.y = As*(- C*(Ai*us) + vs)
# qr_solve!(gs_data, CAi*us - vs)
# temp = gs_data.xs
# S.x = Ai*(us + B*temp)
# S.y = - temp
α = 1e-15
y00 = (S.D - S.CAiB + α * I) \ (- S.C * (S.Ai * u) + v)
r00 = (- S.C * (S.Ai * u) + v) - (S.D - S.CAiB) * y00
y01 = y00 + (S.D - S.CAiB + α * I) \ r00
r01 = (- S.C * (S.Ai * u) + v) - (S.D - S.CAiB) * y01
y02 = y01 + (S.D - S.CAiB + α * I) \ r01
y0 = y00

x0 = S.Ai * (u + S.B * (-y0))
norm(M1*[x0; y0] - [u; v], Inf)
@test norm(M1*[x0; y0] - [u; v], Inf) < 1e-10

diag(S.D - S.CAiB)

cond(S.D)
cond(S.CAiB)
cond(S.D - S.CAiB + 0e-10*Diagonal([zeros(12); ones(2); zeros(2)]))
cond(S.D - S.CAiB + 1e-1*Diagonal([zeros(12); ones(2); zeros(2)]))
cond(ones(ny) .* (S.D - S.CAiB))
cond(1 ./ abs.(diag(S.D - S.CAiB)) .* (S.D - S.CAiB))
sca = Vector(ny ./ sum(abs.(S.D - S.CAiB), dims=2)[:,1])
sca = Vector(ny ./ sum(abs.(S.D - S.CAiB), dims=1)[1,:])
cond(sca .* (S.D - S.CAiB))


v = [1,2]'
v .* A
A = [1 2;
     3 4]
sum(A, dims = 2)



plot(Gray.(1e10*abs.(Matrix(S.D))))
plot(Gray.(1e10*abs.(Matrix(S.CAiB))))
plot(Gray.(1e10*abs.(Matrix(S.D - S.CAiB))))

plot(Gray.(1e0*abs.(Matrix(S.D))))
plot(Gray.(1e0*abs.(Matrix(S.CAiB))))
plot(Gray.(1e0*abs.(Matrix(S.D - S.CAiB))))
