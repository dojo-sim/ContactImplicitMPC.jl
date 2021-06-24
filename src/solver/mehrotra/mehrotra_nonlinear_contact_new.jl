using Random
using LinearAlgebra

mutable struct ProblemData13{T}
    n::Int
    m::Int
    E::AbstractMatrix{T}
    F::AbstractMatrix{T}
    G::AbstractMatrix{T}
    H::AbstractMatrix{T}
    J::AbstractMatrix{T}
    b::AbstractVector{T}
    c::AbstractVector{T}
    B::AbstractMatrix{T}
    idyn::AbstractVector{Int}
    irst::AbstractVector{Int}
    ibil::AbstractVector{Int}
    ix::AbstractVector{Int}
    iy1::AbstractVector{Int}
    iy2::AbstractVector{Int}
end

function get_data(s::Simulation, ref_traj::ContactTraj, t::Int; κ=1e-6)
    get_data(s, ref_traj.z[t], ref_traj.θ[t], κ=κ)
end

function get_data(s::Simulation, z, θ;
    κ=1e-6)
    model = s.model
    env = s.env
    # Dimensions and indices
    ix, iy1, iy2 = linearization_var_index(model, env)
    idyn, irst, ibil, ialt = linearization_term_index(model, env)
    nz = num_var(model, env)
    nθ = num_data(model)
    nx = length(ix)
    ny = length(iy1)

    r = zeros(nz)
    rz = zeros(nz, nz)
    rθ = zeros(nz, nθ)

    s.res.r!(r, z, θ, κ)
    s.res.rz!(rz, z, θ)
    s.res.rθ!(rθ, z, θ)

    # Problem data
    RZ = rz[[idyn; irst; ibil], [ix; iy1; iy2]]
    Rθ = rθ[[idyn; irst; ibil], :]

    E = rz[idyn, ix]
    F = rz[idyn, iy1]
    G = rz[irst, ix]
    H = rz[irst, iy1]
    J = rz[irst, iy2]
    b = r[idyn]
    c = r[irst]
    B = rθ[[idyn; irst], :]

    # Dimensions
    n = nx
    m = ny
    data = ProblemData13(n, m, E, F, G, H, J, b, c, B, idyn, irst, ibil, ix, iy1, iy2)
    return data
end

function residual(w1, w2, w3, z, δθ, data::ProblemData13; κ=0.0, verbose=false)
    δw1 = w1 - z[data.ix]
    δw2 = w2 - z[data.iy1]
    δw3 = w3 - z[data.iy2]
    # verbose && println("δw1: ", scn(norm(δw1)))
    # verbose && println("δw2: ", scn(norm(δw2)))
    # verbose && println("δw3: ", scn(norm(δw3)))
    # verbose && println("data.b: ", scn(norm(data.b, Inf)))
    # verbose && println("data.c: ", scn(norm(data.c, Inf)))
    # verbose && println("bil: ", scn(norm(z[data.iy1] .* z[data.iy2] .- κ, Inf)))
    n = data.n
    m = data.m

    A = [data.E data.F zeros(n,m);
         data.G data.H data.J]
    δwt = [δw1; δw2; δw3]
    x = [data.b; data.c] .- data.B*δθ
    r = [data.E*δw1 + data.F*δw2              - data.b;
         data.G*δw1 + data.H*δw2 + data.J*δw3 - data.c]
    r += data.B*δθ

    v1 = maximum(abs.(r))
    v2 = maximum(abs.(w2 .* w3 .- κ))
    v3 = maximum(abs.(min.(0, w2)))
    v4 = maximum(abs.(min.(0, w3)))
    # verbose && println("res: A", scn(norm(A), digits=4))
    # verbose && println("res: x", scn(norm(x), digits=4))
    # verbose && println("res: δwt", scn(norm(δwt), digits=4))
    # verbose && println("res: Aδwt", scn(norm(A*δwt), digits=4))
    # verbose && println("res: Aδwt-x", scn(norm(A*δwt .- x), digits=4))
    # verbose && println("res: r", scn(norm(r), digits=4))
    verbose && println("res: ", scn.([v1,v2,v3,v4]))
    return maximum([v1,v2,v3,v4])
end

function residual_light(w1, w2, w3, z, δθ, data::ProblemData13; κ=0.0, verbose=false)
    δw1 = w1 - z[data.ix]
    δw2 = w2 - z[data.iy1]
    δw3 = w3 - z[data.iy2]
    # verbose && println("δw1: ", scn(norm(δw1)))
    # verbose && println("δw2: ", scn(norm(δw2)))
    # verbose && println("δw3: ", scn(norm(δw3)))
    # verbose && println("data.b: ", scn(norm(data.b, Inf)))
    # verbose && println("data.c: ", scn(norm(data.c, Inf)))
    # verbose && println("bil: ", scn(norm(z[data.iy1] .* z[data.iy2] .- κ, Inf)))
    n = data.n
    m = data.m

    A = [data.E data.F zeros(n,m);
         data.G data.H data.J]
    δwt = [δw1; δw2; δw3]
    x = [data.b; data.c] .- data.B*δθ
    r = [data.E*δw1 + data.F*δw2              - data.b;
         data.G*δw1 + data.H*δw2 + data.J*δw3 - data.c]
    r += data.B*δθ


    v1 = maximum(abs.(r))
    v2 = maximum(abs.(w2 .* w3 .- κ))
    verbose && println("res: A", scn(norm(A), digits=4))
    verbose && println("res: x", scn(norm(x), digits=4))
    verbose && println("res: δwt", scn(norm(δwt), digits=4))
    verbose && println("res: Aδwt", scn(norm(A*δwt), digits=4))
    verbose && println("res: Aδwt-x", scn(norm(A*δwt .- x), digits=4))
    verbose && println("res: r", scn(norm(r), digits=4))
    verbose && println("res: ", scn.([v1,v2]))
    return maximum([v1,v2])
end

function affine_direction(w1, w2, w3, z, δθ, data::ProblemData13; κ=0.0, verbose = false)
    n = data.n
    m = data.m
    δw1 = w1 - z[data.ix]
    δw2 = w2 - z[data.iy1]
    δw3 = w3 - z[data.iy2]

    r_ = [- (data.E*δw1 + data.F*δw2 - data.b);
          - (data.G*δw1 + data.H*δw2 + data.J*δw3 - data.c);]
    r_ += - data.B*δθ
    r = [r_;
         - (w2 .* w3 .- κ)]

    A_ = [data.E data.F zeros(n,m);
          data.G data.H data.J;]

    M = [data.E data.F zeros(n,m);
         data.G data.H data.J;
         zeros(m,n) Diagonal(w3) Diagonal(w2)]

    verbose && println("**** A:", scn(norm(A_), digits=4))
    verbose && println("**** M:", scn(norm(M), digits=4))
    verbose && println("**** My1:", scn(norm(w2), digits=4))
    verbose && println("**** My2:", scn(norm(w3), digits=4))
    verbose && println("**** My1:", scn.(w2, digits=1))
    verbose && println("**** My2:", scn.(w3, digits=1))
    verbose && println("**** rdyn:", scn(norm(r[1:n]), digits=4))
    verbose && println("**** rrst:", scn(norm(r[n .+ (1:m)]), digits=4))
    verbose && println("**** rbil:", scn(norm(r[n+m .+ (1:m)]), digits=4))

    Δ = M \ r

    err = M * Δ .- r
    err2 = M * Δ .- [r[1:n]; zeros(2m)]
    verbose && println("**** err:", scn(norm(err), digits=4))
    verbose && println("**** err2:", scn(norm(err2), digits=4))
    verbose && println("**** errx:", scn(norm(err[data.ix]), digits=4))
    verbose && println("**** erry1:", scn(norm(err[data.iy1]), digits=4))
    verbose && println("**** erry2:", scn(norm(err[data.iy2]), digits=4))
    verbose && println("**** err:", scn.(err, digits=1))

    return Δ
end

function corrector_direction(w1, w2, w3, z, δθ, σ, μ, Δw2aff, Δw3aff, data::ProblemData13)
    n = data.n
    m = data.m
    δw1 = w1 - z[data.ix]
    δw2 = w2 - z[data.iy1]
    δw3 = w3 - z[data.iy2]

    r_ = [- (data.E*δw1 + data.F*δw2 - data.b);
          - (data.G*δw1 + data.H*δw2 + data.J*δw3 - data.c);]
    r_ += - data.B*δθ
    r = [r_;
         - w2 .* w3 - Δw2aff .* Δw3aff .+ σ * μ]

    M = [data.E data.F zeros(n,m);
         data.G data.H data.J;
         zeros(m,n) Diagonal(w3) Diagonal(w2)]
    Δ = M \ r
    return Δ
end

function unpack(Δ; n=n, m=m)
    off = 0
    w1 = Δ[off .+ (1:n)]; off += n
    w2 = Δ[off .+ (1:m)]; off += m
    w3 = Δ[off .+ (1:m)]; off += m
    return w1, w2, w3
end

function progress_(res; ϵ0 = 0.05)
    # @warn "modified progress"
    # ϵ = min(ϵ0, sqrt(res))
    ϵ = min(ϵ0, res^2)
    τ = 1 - ϵ
    return τ
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
    αh = min(ατ_p, ατ_d)
    return αh
end

function step_length_(w2, w3, Δw2aff, Δw3aff)
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

function initial_state(z, δθ, data::ProblemData13; verbose = false)
    n = data.n
    m = data.m

    A = [data.E data.F zeros(n,m);
         data.G data.H data.J]
    verbose && println("****  A:", scn(norm(A), digits=4))
    x = [data.b; data.c] - data.B*δθ
    verbose && println("****  x:", scn(norm(x), digits=4))
    wt = A' * ((A * A') \ x)
    verbose && println("**** wt:", scn(norm(wt), digits=4))
    verbose && println("**** Awt-x:", scn(norm(A*wt .- x), digits=4))

    w1t, w2t, w3t = unpack(wt, n=n, m=m)
    w1t += z[data.ix]
    w2t += z[data.iy1]
    w3t += z[data.iy2]
    # @warn "early return"
    # return w1t, w2t, w3t

    verbose && println("**** z+wt:", scn(norm([w1t; w2t; w3t]), digits=4))

    δw2 = max(-1.5 * minimum(w2t), 0)
    δw3 = max(-1.5 * minimum(w3t), 0)
    verbose && println("**** δw2:", scn(norm(δw2), digits=4))
    verbose && println("**** δw3:", scn(norm(δw3), digits=4))

    w2h = w2t .+ δw2
    w3h = w3t .+ δw3
    verbose && println("**** w2h:", scn.(w2h[1:3], digits=4))
    verbose && println("**** w3h:", scn.(w3h[1:3], digits=4))

    δhw2 = 0.5 * w2h'*w3h / sum(w3h)
    δhw3 = 0.5 * w2h'*w3h / sum(w2h)
    verbose && println("****  dot:", scn(norm(0.5 * w2h'*w3h), digits=4))
    verbose && println("**** δhw2:", scn(norm(δhw2), digits=4))
    verbose && println("**** δhw3:", scn(norm(δhw3), digits=4))

    w10 = w1t
    w20 = w2h .+ δhw2
    w30 = w3h .+ δhw3
    return w10, w20, w30
end

function initial_state2(z, Δaff, ix, iy1, iy2)

    w1t = z[ix]
    w2t = z[iy1]
    w3t = z[iy2]
    w2t = w2t + Δaff[iy1]
    w3t = w3t + Δaff[iy2]

    δw2 = max(-1.5 * minimum(w2t), 0)
    δw3 = max(-1.5 * minimum(w3t), 0)

    w2h = w2t .+ δw2
    w3h = w3t .+ δw3

    δhw2 = 0.5 * w2h'*w3h / sum(w3h)
    δhw3 = 0.5 * w2h'*w3h / sum(w2h)

    w10 = w1t
    w20 = w2h .+ δhw2
    w30 = w3h .+ δhw3
    return w10, w20, w30
end

function initial_reg(w1t, w2t, w3t)

    δw2 = max(-1.5 * minimum(w2t), 0)
    δw3 = max(-1.5 * minimum(w3t), 0)

    w2h = w2t .+ δw2
    w3h = w3t .+ δw3

    δhw2 = 0.5 * w2h'*w3h / sum(w3h)
    δhw3 = 0.5 * w2h'*w3h / sum(w2h)

    w10 = w1t
    w20 = w2h .+ δhw2
    w30 = w3h .+ δhw3
    return w10, w20, w30
end

function mehrotra_solve(s::Simulation, ref_traj::ContactTraj, t::Int;
        newton_iter::Int=100, res_tol=1e-8, verbose::Bool=false, verbose_res = false, θamp=1.0, ϵ0=0.05)
    κ = res_tol

    z = deepcopy(ref_traj.z[t])
    θ = deepcopy(ref_traj.θ[t])
    # δθ = θamp * θ
    @warn "changed δθ"
    # δθ = θamp * (0.5 .- rand(length(θ)-0))
    # δθ[end-1:end] .= max.(0.0, δθ[end-1:end])
    δθ = θamp * [ones(length(θ) - 2); zeros(2)]

    data = get_data(s, z, θ, κ=κ)
    ix, iy1, iy2 = linearization_var_index(s.model, s.env)

    n = data.n
    m = data.m

    res = residual(z[ix], z[iy1], z[iy2], z, δθ, data, verbose=verbose_res)
    verbose && println("res: ", scn(res))

    # w1, w2, w3 = ref_traj.z[t][ix], ref_traj.z[t][iy1], ref_traj.z[t][iy2]
    # w1, w2, w3 = initial_reg(w1, w2, w3)

    # Δaff = affine_direction(1e-0*ones(n), 1e-0*ones(m), 1e-0*ones(m), z, δθ, deepcopy(data), κ=0.0)
    # w1, w2, w3 = initial_state2(deepcopy(ref_traj.z[t]), Δaff, ix, iy1, iy2)


    verbose && println("****  δθ:", scn(norm(δθ), digits=4))
    verbose && println("****  δθ+θ:", scn(norm(δθ+θ), digits=4))
    verbose && println("****  θ:", scn(norm(θ), digits=4))
    verbose && println("****  z:", scn(norm(z), digits=4))

    resl = residual_light(z[ix], z[iy1], z[iy2], z, δθ, data, verbose=verbose_res)
    verbose && println("**** rl:", scn(resl, digits=4))

    w1, w2, w3 = initial_state(z, δθ, data, verbose = verbose)
    verbose && println("**** w1:", scn(norm(w1), digits=4))
    verbose && println("**** w2:", scn(norm(w2), digits=4))
    verbose && println("**** w3:", scn(norm(w3), digits=4))


    # w1, w2, w3 = unpack(1e-0ones(n + 2m), n=n, m=m)
    rinit = residual_light(w1, w2, w3, z, δθ, data, verbose=verbose_res)
    verbose && println("**** rinit:", scn(rinit, digits=4))

    iter = 0
    success = false
    for k = 1:newton_iter
        iter += 1

        Δaff = affine_direction(w1, w2, w3, z, δθ, data, κ=0.0, verbose = verbose)
        Δw1aff, Δw2aff, Δw3aff = unpack(Δaff, n=n, m=m)
        # verbose && println("**** κ:", scn(min(0.05, res^3), digits=4))
        verbose && println("**** Δaff1:", scn(norm(Δw1aff), digits=4))
        verbose && println("**** Δaff2:", scn(norm(Δw2aff), digits=4))
        verbose && println("**** Δaff3:", scn(norm(Δw3aff), digits=4))
        verbose && println("**** Δaff2:", scn.(Δw2aff, digits=1))
        verbose && println("**** Δaff3:", scn.(Δw3aff, digits=1))
        αhaff, μaff = step_length_(w2, w3, Δw2aff, Δw3aff)
        μ = w2'*w3 / m
        σ = (μaff / μ)^3
        Δ = corrector_direction(w1, w2, w3, z, δθ, σ, μ, Δw2aff, Δw3aff, data)
        Δw1, Δw2, Δw3 = unpack(Δ, n=n, m=m)
        res = residual(w1, w2, w3, z, δθ, data)
        τ = progress_(res, ϵ0=ϵ0)
        αh = corrector_step_length(w2, w3, Δw2, Δw3; τ=τ)
        w1 = w1 + αh * Δw1
        w2 = w2 + αh * Δw2
        w3 = w3 + αh * Δw3
        verbose && println("**** Δ1:", scn(norm(Δw1 * αh), digits=4))
        verbose && println("**** Δ2:", scn(norm(Δw2 * αh), digits=4))
        verbose && println("**** Δ3:", scn(norm(Δw3 * αh), digits=4))
        res = residual(w1, w2, w3, z, δθ, data, verbose=verbose_res)
        verbose && println("res: ", scn(res))
        # @show res
        if res < res_tol
            success = true
            break
        end
    end
    println("meh zdiff: ", scn(z_difference(z, w1, w2, w3, data)))
    res = residual(w1, w2, w3, z, δθ, data)
    verbose && println("iter: ", iter, " res: ", scn(res))
    return iter, success
end

function baseline_solve(s::Simulation, ref_traj::ContactTraj, t::Int; newton_iter::Int=100,
        res_tol=1e-8, verbose::Bool=false, verbose_res = false, θamp=1.0, ϵ0 = 0.05)
    κ = res_tol

    z = deepcopy(ref_traj.z[t])
    θ = deepcopy(ref_traj.θ[t])
    # δθ = θamp * θ
    @warn "changed δθ"
    δθ = θamp * (0.5 .- rand(length(θ)-0))
    δθ[end-1:end] .= max.(0.0, δθ[end-1:end])
    data = get_data(s, z, θ, κ=κ)
    ix, iy1, iy2 = linearization_var_index(s.model, s.env)
    n = data.n
    m = data.m

    res = residual(z[ix], z[iy1], z[iy2], z, δθ, data, verbose=verbose_res)
    verbose && println("res: ", scn(res))

    w1, w2, w3 = initial_state(z, δθ, data, verbose = verbose)
    # w1, w2, w3 = unpack(1e-0ones(n + 2m), n=n, m=m)

    rinit = residual_light(w1, w2, w3, z, δθ, data, verbose=verbose_res)
    verbose && println("**** rinit:", scn(rinit, digits=4))

    iter = 0
    success = false
    for k = 1:newton_iter
        iter += 1

        Δaff = affine_direction(w1, w2, w3, z, δθ, data, κ=κ, verbose = verbose)
        Δw1aff, Δw2aff, Δw3aff = unpack(Δaff, n=n, m=m)
        αhaff, μaff = step_length_(w2, w3, Δw2aff, Δw3aff)
        w1 = w1 + αhaff * Δw1aff
        w2 = w2 + αhaff * Δw2aff
        w3 = w3 + αhaff * Δw3aff
        verbose && println("**** Δ1:", scn(norm(Δw1aff * αhaff), digits=4))
        verbose && println("**** Δ2:", scn(norm(Δw2aff * αhaff), digits=4))
        verbose && println("**** Δ3:", scn(norm(Δw3aff * αhaff), digits=4))
        verbose && println("**** w1:", scn(norm(w1), digits=4))
        verbose && println("**** w2:", scn(norm(w2), digits=4))
        verbose && println("**** w3:", scn(norm(w3), digits=4))

        verbose && println("****    Δ2:", scn.(Δw2aff * αhaff, digits=1))
        verbose && println("****    Δ3:", scn.(Δw3aff * αhaff, digits=1))
        verbose && println("****    w2:", scn.(w2, digits=1))
        verbose && println("****    w3:", scn.(w3, digits=1))
        verbose && println("**** w2*w3:", scn.(w3.*w2, digits=1))

        res = residual(w1, w2, w3, z, δθ, data; κ=κ, verbose=verbose_res)
        verbose && println("res: ", scn(res))
        if res < res_tol
            success = true
            break
        end
    end
    println("bas zdiff: ", scn(z_difference(z, w1, w2, w3, data)))
    res = residual(w1, w2, w3, z, δθ, data; κ=κ, verbose=verbose_res)
    verbose && println("iter: ", iter, " res: ", scn(res))
    return iter, success
end

function z_difference(z, w1, w2, w3, data)
    return norm(z[[data.ix; data.iy1; data.iy2]] - [w1; w2; w3])
end

function evaluate(s::Simulation, ref_traj::ContactTraj;
        algorithm::Symbol=:mehrotra_solve, newton_iter=100, res_tol=1e-8, θamp=1.0, verbose::Bool = false, ϵ0=0.05)

    κ = res_tol
    H = ref_traj.H
    iter = []
    success = 0.0
    for t = 1:H
        it = newton_iter
        su = false
        it, su = eval(algorithm)(s, deepcopy(ref_traj), t,
            newton_iter = newton_iter,
            res_tol = res_tol,
            θamp = θamp,
            verbose = verbose, ϵ0 = ϵ0)
        try
            it, su = eval(algorithm)(s, deepcopy(ref_traj), t,
                newton_iter = newton_iter,
                res_tol = res_tol,
                θamp = θamp,
                verbose = verbose, ϵ0 = ϵ0)
        catch e
            @show e
        end
        push!(iter, it)
        success += su
    end

    μ_iter = mean(iter)
    σ_iter = (iter .- mean(iter)).^2 / H
    σ_iter = sqrt(sum(σ_iter))
    success_rate = success / H
    return μ_iter, σ_iter, success_rate
end




# # Test
# s = get_simulation("quadruped", "flat_2D_lc", "flat")
# ref_traj = deepcopy(get_trajectory(s.model, s.env,
#     joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
#     load_type = :split_traj_alt))
#
# data = get_data(s.model, s.env, ref_traj, t=2)
# @time mehrotra_solve(data, newton_iter=30)
# @time baseline_solve(data, newton_iter=30)

################################################################################
# Benchmark
################################################################################

# Test Quadruped
s1 = get_simulation("quadruped", "flat_2D_lc", "flat")
ref_traj1 = deepcopy(get_trajectory(s1.model, s1.env,
    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
    load_type = :split_traj_alt))

# Test Flamingo
s2 = get_simulation("flamingo", "flat_2D_lc", "flat")
ref_traj2 = deepcopy(get_trajectory(s2.model, s2.env,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))

# Test Hopper 2D
s3 = get_simulation("hopper_2D", "flat_2D_lc", "flat")
ref_traj3 = get_trajectory(s3.model, s3.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_in_place.jld2"),
    load_type=:joint_traj)

θamp = 1e-1
res_tol = 1e-8
μ_iter, σ_iter, success_rate = evaluate(s1, ref_traj1, algorithm = :mehrotra_solve,
                                        res_tol=res_tol, θamp=θamp)
μ_iter, σ_iter, success_rate = evaluate(s1, ref_traj1, algorithm = :baseline_solve,
                                        res_tol=res_tol, θamp=θamp)

μ_iter, σ_iter, success_rate = evaluate(s2, ref_traj2, algorithm = :mehrotra_solve,
                                        res_tol=res_tol, θamp=θamp)
μ_iter, σ_iter, success_rate = evaluate(s2, ref_traj2, algorithm = :baseline_solve,
                                        res_tol=res_tol, θamp=θamp)

μ_iter, σ_iter, success_rate = evaluate(s3, ref_traj3, algorithm = :mehrotra_solve,
                                        res_tol=res_tol, θamp=θamp, verbose=true)
μ_iter, σ_iter, success_rate = evaluate(s3, ref_traj3, algorithm = :baseline_solve,
                                        res_tol=res_tol, θamp=θamp, verbose=true)


function benchmark_mehrotra(ss, ref_trajs; plt = plot(), ϵ0 = 0.05, color=:red,
        ls=[:solid, :dash, :dot])
    n = length(ss)
    μ = [[] for i=1:n]
    σ = [[] for i=1:n]
    su = [[] for i=1:n]
    plot!(plt,
        xaxis=:log,
        legend=false,
        ylims=[0,20],
        xlabel="data pertubation magnitude",
        ylabel="iterations")
    θamps = 10 .^ (range(-4,stop=2,length=20))
    for θamp in θamps
        for i in eachindex(ss)
            s = ss[i]
            ref_traj = ref_trajs[i]
            # μ_iter, σ_iter, success_rate = evaluate(s, ref_traj, algorithm = :mehrotra_solve,
            μ_iter, σ_iter, success_rate = evaluate(s, ref_traj, algorithm = :baseline_solve,
                                                    res_tol=1e-5, θamp=θamp, ϵ0 = ϵ0)
            push!(μ[i], μ_iter)
            push!(σ[i], σ_iter)
            push!(su[i], success_rate)
        end
    end
    for i = 1:n
        plot!(plt, θamps, μ[i],        linewidth=5.0, color=color, linestyle=ls[i])
        # plot!(plt, θamps, μ[i] - σ[i], linewidth=1.0, color=color, linestyle=ls[i])
        # plot!(plt, θamps, μ[i] + σ[i], linewidth=1.0, color=color, linestyle=ls[i])
    end
    display(plt)
    return plt
end

ss = [s1, s2, s3]
ref_trajs = [ref_traj1, ref_traj2, ref_traj3]
plt = benchmark_mehrotra(ss, ref_trajs, ϵ0 = 0.005, color=:black)
plt = benchmark_mehrotra(ss, ref_trajs, ϵ0 = 0.01, color=:red, plt = plt)
plt = benchmark_mehrotra(ss, ref_trajs, ϵ0 = 0.02, color=:blue, plt = plt)
plt = benchmark_mehrotra(ss, ref_trajs, ϵ0 = 0.05, color=:yellow, plt = plt)
plot([1,2,3], linestyle=:dash)

data = get_data(s1, ref_traj1, 5, κ=1e-8)
# clearconsole()

# the residual should all be zeros
function residual_list(s::Simulation, ref_traj::ContactTraj; verbose::Bool = false)
    l = []
    for t = 1:ref_traj.H
        res = residual(s.model, s.env, ref_traj.z[t], ref_traj.θ[t], 0.0)
        push!(l, norm(res, Inf))
        verbose && println("res: ", scn(norm(res , Inf)))
    end
    return l
end
plot(log.(10, residual_list(s1, ref_traj1)), legend=false)
plot(log.(10, residual_list(s2, ref_traj2)), legend=false)
plot(log.(10, residual_list(s3, ref_traj3)), legend=false)






# Test
s = get_simulation("flamingo", "flat_2D_lc", "flat")
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))

Random.seed!(100)
@time mehrotra_solve(s, ref_traj, 1, newton_iter=30, res_tol=1e-8, verbose=true, θamp=1e-0, ϵ0=0.05)
@time baseline_solve(s, ref_traj, 1, newton_iter=30, res_tol=1e-8, verbose=true, θamp=0e-0, ϵ0=0.05)
mean([norm(θ, 1)/length(θ) for θ in ref_traj1.θ])


cnt = 0
itl = []
sul = 0
for t = 1:100
    Random.seed!(t)
    it, su = mehrotra_solve(s, ref_traj, 1, newton_iter=60, res_tol=1e-8, verbose=true, θamp=5e+1, ϵ0=0.05)
    # it, su = baseline_solve(s, ref_traj, 1, newton_iter=60, res_tol=1e-8, verbose=true, θamp=5e-0, ϵ0=0.05)
    sul += Int(su)
    push!(itl, it)
end
sul
mean(itl)

n = 10
m = 5
x = rand(n)
b = rand(m)
A = rand(m,n)
x0 = A \ b
x1 = A' * ((A * A') \ b )
x0 - x1
