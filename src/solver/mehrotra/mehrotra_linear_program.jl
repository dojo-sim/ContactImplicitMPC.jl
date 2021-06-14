using Random
using LinearAlgebra

# Dimensions
n = 10
m = 2

Random.seed!(100)
A = rand(m,n)
b = rand(m)
c = rand(n)

x = 1e-3*ones(n)
λ = 1e-3*ones(m)
s = 1e-3*ones(n)
μ = 1e-6


function residual(x, λ, s; A=A, b=b, c=c)
    v1 = maximum(abs.(A'*λ + s - c))
    v2 = maximum(abs.(A*x - b))
    v3 = maximum(abs.(x .* s))
    v4 = maximum(abs.(min.(0, x)))
    v5 = maximum(abs.(min.(0, s)))
    # @show scn.([v1,v2,v3,v4,v5])
    return maximum([v1,v2,v3,v4,v5])
end

function affine_direction(x, λ, s; A=A, b=b, c=c)
    m, n = size(A)
    rb = A*x - b
    rc = A'*λ + s - c
    r = [-rc; -rb; -x .* s]
    H = [zeros(n,n) A' I;
         A zeros(m,m) zeros(m,n);
         Diagonal(s) zeros(n,m) Diagonal(x)]
    Δ = H \ r
    return Δ
end

function corrector_direction(x, λ, s, σ, μ, Δxaff, Δsaff; A=A, b=b, c=c)
    m, n = size(A)
    rb = A*x - b
    rc = A'*λ + s - c
    r = [-rc; -rb; - x .* s - Δxaff .* Δsaff .+ σ * μ]
    H = [zeros(n,n) A' I;
         A zeros(m,m) zeros(m,n);
         Diagonal(s) zeros(n,m) Diagonal(x)]
    Δ = H \ r
    return Δ
end

function unpack(Δ; n=n, m=m)
    off = 0
    x = Δ[off .+ (1:n)]; off += n
    λ = Δ[off .+ (1:m)]; off += m
    s = Δ[off .+ (1:n)]; off += n
    return x, λ, s
end

function pc_solve(; A=A, b=b, c=c, newton_iter::Int=100, μ=1e-1, η=0.9)
    m, n = size(A)
    x, λ, s = initial_state(A=A, b=b, c=c)
    for k = 1:newton_iter
        Δaff = affine_direction(x, λ, s, A=A, b=b, c=c)
        Δxaff, Δλaff, Δsaff = unpack(Δaff, n=n, m=m)
        αaff_p, αaff_d, μaff = step_length(x, s, Δxaff, Δsaff)
        σ = (μaff / μ)^3
        @show scn(σ)
        Δ = corrector_direction(x, λ, s, σ, μ, Δxaff, Δsaff, A=A, b=b, c=c)
        α_p = min(1, η*αaff_p)
        α_d = min(1, η*αaff_d)
        @show α_p
        @show α_d
        Δx, Δλ, Δs = unpack(Δ, n=n, m=m)
        x = x + α_p * Δx
        λ = λ + α_d * Δλ
        s = s + α_d * Δs
        res = residual(x, λ, s, A=A, b=b, c=c)
        println("res :     ", scn(res))
        if res < 1e-8
            break
        end
    end
    return x, λ, s
end

function step_length(x, s, Δxaff, Δsaff)
    n = length(x)

    αaff_p = 1.0
    αaff_d = 1.0
    for i = 1:n
        if Δxaff[i] < 0.0
            αaff_p = min(αaff_p, - x[i] / Δxaff[i])
        end
        if Δsaff[i] < 0.0
            αaff_d = min(αaff_d, - s[i] / Δsaff[i])
        end
    end
    μaff = (x + αaff_p * Δxaff)' * (s + αaff_d * Δsaff) / n
    return αaff_p, αaff_d, μaff
end

function initial_state(;A=A, b=b, c=c)
    xt = A' * ((A * A') \ b)
    λt = (A * A') \ (A * c)
    st = c - A'*λt

    δx = max(-1.5 * minimum(xt), 0)
    δs = max(-1.5 * minimum(st), 0)

    xh = xt .+ δx
    sh = st .+ δs

    δhx = 0.5 * xh'*sh / sum(sh)
    δhs = 0.5 * xh'*sh / sum(xh)

    x0 = xh .+ δhx
    λ0 = λt
    s0 = sh .+ δhs
    # @show δx
    # @show δs
    # @show norm(xt - A' * inv(A * A') * b)
    # @show norm(λt - inv(A * A') * (A * c))
    #
    # @show x0
    # @show λ0
    # @show s0
    #
    # @show δhx
    # @show δhs
    return x0, λ0, s0
end

residual(x, λ, s)
initial_state(A=A, b=b, c=c)
pc_solve(newton_iter=20)
