@testset "Solver: LDL" begin
    Random.seed!(100)
    n = 50
    d = 0.7

    # Generate diffrent matrices with the same sparsity pattern
    A0_ = sprand(n, n, d)
    A1_ = deepcopy(A0_)
    A1_.nzval[1:10] .+= 1.0
    A0 = A0_ + A0_'
    A1 = A1_ + A1_'
    @test rank(A0) == n
    @test rank(A1) == n
    # Generate random vector
    b = rand(n)

    # QDLDL works
    F0 = qdldl(A0)
    x0 = deepcopy(b)
    QDLDL.solve!(F0, x0)
    @test norm(A0 * x0 - b, Inf) < 1e-8

    # The in-place QDLDL works on the same matrix A0
    G0 = ContactImplicitMPC.LDLSolver(A0, F0)
    ContactImplicitMPC.factorize!(G0, A0)
    x0 = deepcopy(b)
    QDLDL.solve!(G0.F, x0)
    @test norm(A0 * x0 - b, Inf) < 1e-8

    # The in-place QDLDL works on a different marix A1 but with the same sparsity pattern
    G1 = ContactImplicitMPC.LDLSolver(A0, F0)
    ContactImplicitMPC.factorize!(G0, A1)
    x0 = deepcopy(b)
    QDLDL.solve!(G0.F, x0)
    @test !(norm(A0 * x0 - b, Inf) < 1e-8)
    @test norm(A1 * x0 - b, Inf) < 1e-8
end
