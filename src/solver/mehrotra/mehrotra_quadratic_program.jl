using Random
using LinearAlgebra

# Dimensions
n = 10
m = 5

Random.seed!(100)
G = rand(n,n)
G = G'*G
A = rand(m,n)
c = rand(n)
b = rand(m)

x = 1e-3*ones(n)
λ = 1e-3*ones(m)
y = 1e-3*ones(m)
μ = 1e-6


function residual(x, λ, y; G=G, A=A, b=b, c=c)
    v1 = maximum(abs.(G*x - A'*λ + c))
    v2 = maximum(abs.(A*x - y - b))
    v3 = maximum(abs.(y .* λ))
    v4 = maximum(abs.(min.(0, y)))
    v5 = maximum(abs.(min.(0, λ)))
    # @show scn.([v1,v2,v3,v4,v5])
    return maximum([v1,v2,v3,v4,v5])
end

function affine_direction(x, λ, y; G=G, A=A, b=b, c=c)
    m, n = size(A)
    rd = G*x - A'*λ + c
    rp = A*x - y - b
    r = [-rd; -rp; -λ .* y]
    H = [G zeros(n,m) -A';
         A -I zeros(m,m);
         zeros(m,n) Diagonal(λ) Diagonal(y)]
    Δ = H \ r
    return Δ
end

function corrector_direction(x, λ, y, σ, μ, Δyaff, Δλaff; G=G, A=A, b=b, c=c)
    m, n = size(A)
    rd = G*x - A'*λ + c
    rp = A*x - y - b
    r = [-rd; -rp; -λ .* y - Δλaff .* Δyaff .+ σ * μ]
    H = [G zeros(n,m) -A';
         A -I zeros(m,m);
         zeros(m,n) Diagonal(λ) Diagonal(y)]
    Δ = H \ r
    return Δ
end

function unpack(Δ; n=n, m=m)
    off = 0
    x = Δ[off .+ (1:n)]; off += n
    y = Δ[off .+ (1:m)]; off += m
    λ = Δ[off .+ (1:m)]; off += m
    return x, λ, y
end

function pc_solve(; G=G, A=A, b=b, c=c, newton_iter::Int=100, μ=1e-1, η=0.9, τ=0.9995)
    m, n = size(A)
    x, λ, y = initial_state(G=G, A=A, b=b, c=c)
    for k = 1:newton_iter
        Δaff = affine_direction(x, λ, y; G=G, A=A, b=b, c=c)
        Δxaff, Δλaff, Δyaff = unpack(Δaff, n=n, m=m)
        μ = y'*λ / m
        αhaff, μaff = step_length(y, λ, Δyaff, Δλaff)
        σ = (μaff / μ)^3
        Δ = corrector_direction(x, λ, y, σ, μ, Δyaff, Δλaff; G=G, A=A, b=b, c=c)
        Δx, Δλ, Δy = unpack(Δ, n=n, m=m)
        αh = corrector_step_length(y, λ, Δy, Δλ; τ=τ)
        x = x + αh * Δx
        λ = λ + αh * Δλ
        y = y + αh * Δy
        res = residual(x, λ, y, G=G, A=A, b=b, c=c)
        println("res :     ", scn(res))
        if res < 1e-8
            break
        end
    end
    return x, λ, y
end

function corrector_step_length(y, λ, Δy, Δλ; τ=0.9995)
    m = length(y)

    ατ_p = 1.0
    ατ_d = 1.0
    for i = 1:m
        if Δy[i] < 0.0
            ατ_p = min(ατ_p, - τ * y[i] / Δy[i])
        end
        if Δλ[i] < 0.0
            ατ_d = min(ατ_d, - τ * λ[i] / Δλ[i])
        end
    end
    # @show ατ_p
    # @show ατ_d
    αh = min(ατ_p, ατ_d)
    return αh
end

function step_length(y, λ, Δyaff, Δλaff)
    m = length(y)

    αhaff = 1.0
    for i = 1:m
        if Δyaff[i] < 0.0
            αhaff = min(αhaff, - y[i] / Δyaff[i])
        end
        if Δλaff[i] < 0.0
            αhaff = min(αhaff, - λ[i] / Δλaff[i])
        end
    end
    # @show αhaff
    μaff = (y + αhaff * Δyaff)' * (λ + αhaff * Δλaff) / m
    return αhaff, μaff
end

function initial_state(; G=G, A=A, b=b, c=c)
    xt = (A' * A + I) \ (A' * b)
    yt = - b + A * xt
    ct = c + G * xt
    λt = (A *A') \ (A * ct)

    δλ = max(-1.5 * minimum(λt), 0)
    δy = max(-1.5 * minimum(yt), 0)

    λh = λt .+ δλ
    yh = yt .+ δy

    δhλ = 0.5 * λh'*yh / sum(yh)
    δhy = 0.5 * λh'*yh / sum(λh)

    x0 = xt
    λ0 = λh .+ δhλ
    y0 = yh .+ δhy
    return x0, λ0, y0
end

residual(x, λ, y)
initial_state(G=G, A=A, b=b, c=c)
pc_solve(newton_iter=10)
