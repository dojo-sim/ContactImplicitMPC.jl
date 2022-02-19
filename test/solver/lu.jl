@testset "LU solver" begin
    n = 20
    m = 10
    A = rand(n, n)
    X = rand(n, m)
    B = rand(n, m)
    x = rand(n)
    b = rand(n)

    sol = ContactControl.lu_solver(A)
    ContactControl.linear_solve!(sol, X, A, B)
    @test norm(A * X - B, Inf) < 1e-10

    sol = ContactControl.lu_solver(A)
    ContactControl.linear_solve!(sol, x, A, b)
    @test norm(A * x - b, Inf) < 1e-10
end