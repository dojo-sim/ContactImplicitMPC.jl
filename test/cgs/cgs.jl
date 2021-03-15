@testset "Gram Schmidt" begin
    T = Float64
    n = 20
    Random.seed!(100)
    A = sprand(n,n,0.16)
    while rank(A) < n
        A = sprand(n,n,0.16)
    end
    cgs_data = CGSData13!(A,n)

    cgs!(cgs_data, A)
    @test norm(A - hcat(cgs_data.qs...)*triangularize_atom(cgs_data.rs,n), Inf) < 1e-10

    A1 = deepcopy(A)
    A1.nzval .+= rand(nnz(A))
    cgs!(cgs_data, A1)
    @test norm(A1 - hcat(cgs_data.qs...)*triangularize_atom(cgs_data.rs,n), Inf) < 1e-10

    b = rand(SizedVector{n,T})
    c = zeros(SizedVector{n,T})
    qr_solve!(cgs_data, c, b)
    x = A1\b
    @test norm(c - x, Inf) < 1e-10
end
