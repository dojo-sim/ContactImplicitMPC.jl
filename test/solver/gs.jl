@testset "Solver: Classical Gram Schmidt" begin
    T = Float64
    n = 20
    Random.seed!(100)
    A = sprand(n, n, 0.16)
    while rank(A) < n
        A = sprand(n, n, 0.16)
    end
    gs_data = ContactControl.CGSData(A, n)

    ContactControl.gs!(gs_data, A)
    @test norm(A - hcat(gs_data.qs...) * ContactControl.triangularize(gs_data.rs, n), Inf) < 1e-10

    A1 = deepcopy(A)
    A1.nzval .+= rand(nnz(A))
    ContactControl.cgs!(gs_data, A1)
    @test norm(A1 - hcat(gs_data.qs...) * ContactControl.triangularize(gs_data.rs, n), Inf) < 1e-10

    b = rand(SizedVector{n,T})
    c = zeros(SizedVector{n,T})
    ContactControl.qr_solve!(gs_data, c, b)
    x = A1 \ b
    @test norm(c - x, Inf) < 1e-10
end

@testset "Solver: Modified Gram Schmidt" begin
    T = Float64
    n = 20
    Random.seed!(100)
    A = sprand(n, n, 0.16)
    while rank(A) < n
        A = sprand(n, n, 0.16)
    end
    gs_data = ContactControl.MGSData(A, n)

    ContactControl.mgs!(gs_data, A)
    @test norm(A - hcat(gs_data.qs...) * ContactControl.triangularize(gs_data.rs, n), Inf) < 1e-10

    A1 = deepcopy(A)
    A1.nzval .+= rand(nnz(A))
    ContactControl.gs!(gs_data, A1)
    @test norm(A1 - hcat(gs_data.qs...) * ContactControl.triangularize(gs_data.rs, n), Inf) < 1e-10

    b = rand(SizedVector{n,T})
    c = zeros(SizedVector{n,T})
    ContactControl.qr_solve!(gs_data, c, b)
    x = A1 \ b
    @test norm(c - x, Inf) < 1e-10
end


@testset "Solver: Dense Modified Gram Schmidt" begin
    T = Float64
    n = 20
    Random.seed!(100)
    A = sprand(n, n, 0.16)
    while rank(A) < n
        A = sprand(n, n, 0.16)
    end
    Ad = Matrix(A)
    gs_data = ContactControl.DMGSData!(Ad, n)

    ContactControl.dmgs!(gs_data, Ad)
    @test norm(A - hcat(gs_data.qs...) * ContactControl.triangularize(gs_data.rs, n), Inf) < 1e-10
    b = rand(SizedVector{n,T})
    c = zeros(SizedVector{n,T})
    ContactControl.qr_solve!(gs_data, c, b)
    x = Ad \ b
    @test norm(c - x, Inf) < 1e-10
end
