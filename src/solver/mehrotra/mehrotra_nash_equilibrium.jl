using Random
using LinearAlgebra

# Dimensions
n = 11
m = 16

Random.seed!(100)
M = rand(n+m, n+m)
M = M'*M
E = M[1:n, 1:n]
F = M[1:n, n .+ (1:m)]
G = M[n .+ (1:m), 1:n]
H = -M[n .+ (1:m), n .+ (1:m)]
b = rand(n)
c = rand(m)
A = [E F zeros(n,m);
     G H Diagonal(ones(m))]

w1 = 1e-3*ones(n)
w2 = 1e-3*ones(m)
w3 = 1e-3*ones(m)
A * [w1; w2; w3]

function residual(w1, w2, w3; E=E, F=F, G=G, H=H, b=b, c=c)
    v1 = maximum(abs.(E*w1 + F*w2 -b))
    v2 = maximum(abs.(G*w1 + H*w2 + w3 - c))
    v3 = maximum(abs.(w2 .* w3))
    v4 = maximum(abs.(min.(0, w2)))
    v5 = maximum(abs.(min.(0, w3)))
    # @show scn.([v1,v2,v3,v4,v5])
    return maximum([v1,v2,v3,v4,v5])
end

function affine_direction(w1, w2, w3; E=E, F=F, G=G, H=H, b=b, c=c)
    m, n = size(G)

    r = [- (E*w1 + F*w2 - b);
         - (G*w1 + H*w2 + w3 - c);
         - w2 .* w3 .+ 0.0*1e-9]
    H = [E F zeros(n,m);
         G H Diagonal(ones(m));
         zeros(m,n) Diagonal(w3) Diagonal(w2)]
    Δ = H \ r
    return Δ
end

function corrector_direction(w1, w2, w3, σ, μ, Δw2aff, Δw3aff; E=E, F=F, G=G, H=H, b=b, c=c)
    m, n = size(G)

    r = [- (E*w1 + F*w2 -b);
         - (G*w1 + H*w2 + w3 - c);
         - w2 .* w3 - Δw2aff .* Δw3aff .+ σ * μ]
    H = [E F zeros(n,m);
         G H Diagonal(ones(m));
         zeros(m,n) Diagonal(w3) Diagonal(w2)]
    Δ = H \ r
    return Δ
end

function unpack(Δ; n=n, m=m)
    off = 0
    w1 = Δ[off .+ (1:n)]; off += n
    w2 = Δ[off .+ (1:m)]; off += m
    w3 = Δ[off .+ (1:m)]; off += m
    return w1, w2, w3
end

function pc_solve(; E=E, F=F, G=G, H=H, b=b, c=c, newton_iter::Int=100, μ=1e-1, η=0.9, τ=0.99995)
    m, n = size(G)
    w1, w2, w3 = initial_state(E=E, F=F, G=G, H=H, b=b, c=c)
    iter = 0
    for k = 1:newton_iter
        iter += 1
        Δaff = affine_direction(w1, w2, w3; E=E, F=F, G=G, H=H, b=b, c=c)
        Δw1aff, Δw2aff, Δw3aff = unpack(Δaff, n=n, m=m)
        μ = w2'*w3 / m
        αhaff, μaff = step_length(w2, w3, Δw2aff, Δw3aff)
        σ = (μaff / μ)^3
        Δ = corrector_direction(w1, w2, w3, σ, μ, Δw2aff, Δw3aff; E=E, F=F, G=G, H=H, b=b, c=c)
        Δw1, Δw2, Δw3 = unpack(Δ, n=n, m=m)
        # αh = corrector_step_length(w2, w3, Δw2, Δw3; τ=1-0.3^k)
        # αh = corrector_step_length(w2, w3, Δw2, Δw3; τ=1-0.03^k)
        # αh = corrector_step_length(w2, w3, Δw2, Δw3; τ=τ)
        # w1 = w1 + αhaff * Δw1aff
        # w2 = w2 + αhaff * Δw2aff
        # w3 = w3 + αhaff * Δw3aff
        w1 = w1 + αh * Δw1
        w2 = w2 + αh * Δw2
        w3 = w3 + αh * Δw3
        res = residual(w1, w2, w3; E=E, F=F, G=G, H=H, b=b, c=c)
        println("res :     ", scn(res))
        if res < 1e-8
            break
        end
    end
    @show iter
    return w1, w2, w3
end

function corrector_step_length(w2, w3, Δw2, Δw3; τ=0.9995)
    m = length(w2)

    ατ_p = 1.0
    ατ_d = 1.0
    for i = 1:m
        if Δw2[i] < 0.0
            ατ_p = min(ατ_p, - τ * w2[i] / Δw2[i])
        end
        if Δw3[i] < 0.0
            ατ_d = min(ατ_d, - τ * w3[i] / Δw3[i])
        end
    end
    # @show ατ_p
    # @show ατ_d
    αh = min(ατ_p, ατ_d)
    return αh
end

function step_length(w2, w3, Δw2aff, Δw3aff)
    m = length(w2)

    αhaff = 1.0
    for i = 1:m
        if Δw2aff[i] < 0.0
            αhaff = min(αhaff, - w2[i] / Δw2aff[i])
        end
        if Δw3aff[i] < 0.0
            αhaff = min(αhaff, - w3[i] / Δw3aff[i])
        end
    end
    # @show αhaff
    μaff = (w2 + αhaff * Δw2aff)' * (w3 + αhaff * Δw3aff) / m
    return αhaff, μaff
end

function initial_state(; E=E, F=F, G=G, H=H, b=b, c=c)
    m, n = size(G)
    A = [E F zeros(n,m);
         G H Diagonal(ones(m))]
    wt = A' * ((A * A') \ [b; c])
    w1t, w2t, w3t = unpack(wt, n=n, m=m)

    δw2 = max(-1.5 * minimum(w2t), 0)
    δw3 = max(-1.5 * minimum(w3t), 0)

    w2h = w2t .+ δw2
    w3h = w3t .+ δw3

    δhw2 = 0.5 * w2h'*w3h / sum(w3h)
    δhw3 = 0.5 * w2h'*w3h / sum(w2h)

    w10 = w1t
    w20 = w2h .+ δhw2
    w30 = w3h .+ δhw3
    # w10 = ones(n)
    # w20 = ones(m)
    # w30 = ones(m)
    return w10, w20, w30
end

residual(w1, w2, w3)
initial_state(E=E, F=F, G=G, H=H, b=b, c=c)
pc_solve(newton_iter=30)
