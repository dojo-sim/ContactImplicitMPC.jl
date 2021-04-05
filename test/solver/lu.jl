@testset "LU solver" begin
    n = 20
    m = 10
    A = rand(n, n)
    X = rand(n, m)
    B = rand(n, m)
    sol = ContactControl.lu_solver(A)
    ContactControl.linear_solve!(sol, X, A, B)
    @test norm(A * X - B, Inf) < 1e-10
end
