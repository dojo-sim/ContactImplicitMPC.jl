@testset "Schur Complement" begin
    T = Float64
    n = 11
    m = 16
    M = rand(n+m,n+m)
    S = ContactControl.Schur(M, n=n, m=m)

    u = rand(T,n)
    v = rand(T,m)
    D = rand(m,m)
    ContactControl.schur_update!(S, D)
    # @benchmark ContactControl.schur_update!(S, D)
    ContactControl.schur_solve!(S, u, v)
    # @benchmark ContactControl.schur_solve!(S, u, v)
    M1 = [S.A S.B; S.C D]
    @test norm(M1*[S.x; S.y] - [u; v], Inf) < 1e-10
end
