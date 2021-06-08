using Random
using LinearAlgebra

mutable struct ProblemData{T}
    n::Int
    m::Int
    E::AbstractMatrix{T}
    F::AbstractMatrix{T}
    G::AbstractMatrix{T}
    H::AbstractMatrix{T}
    J::AbstractMatrix{T}
    b::AbstractVector{T}
    c::AbstractVector{T}
end

function get_data(s::Simulation, ref_traj::ContactTraj;
    κ=1e-6, t::Int=1, z::AbstractVector=ref_traj.z[t])
    model = s.model
    env = s.env
    # Dimensions and indices
    ix, iy1, iy2 = linearization_var_index(model, env)
    idyn, irst, ibil, ialt = linearization_term_index(model, env)
    nz = num_var(model, env)
    nθ = num_data(model)
    nx = length(ix)
    ny = length(iy1)

    # z = ref_traj.z[t]
    θ = ref_traj.θ[t]

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
    # b = (rz*z + r)[idyn]
    b = r[idyn]
    # c = (rz*z + r)[irst]
    c = r[irst]

    # Dimensions
    n = nx
    m = ny
    data = ProblemData(n, m, E, F, G, H, J, b, c)
    return data
end

function residual(w1, w2, w3, data::ProblemData; κ=0.0, verbose=false)
    v1 = maximum(abs.(data.E*w1 + data.F*w2 - data.b))
    v2 = maximum(abs.(data.G*w1 + data.H*w2 + data.J*w3 - data.c))
    v3 = maximum(abs.(w2 .* w3 .- κ))
    v4 = maximum(abs.(min.(0, w2)))
    v5 = maximum(abs.(min.(0, w3)))
    verbose && println("res: ", scn.([v1,v2,v3,v4,v5]))
    return maximum([v1,v2,v3,v4,v5])
end

function affine_direction(w1, w2, w3, data::ProblemData; κ=0.0)
    n = data.n
    m = data.m

    r = [- (data.E*w1 + data.F*w2 - data.b);
         - (data.G*w1 + data.H*w2 + data.J*w3 - data.c);
         - w2 .* w3 .+ κ]
    M = [data.E data.F zeros(n,m);
         data.G data.H data.J;
         zeros(m,n) Diagonal(w3) Diagonal(w2)]
    Δ = M \ r
    return Δ
end

function corrector_direction(w1, w2, w3, σ, μ, Δw2aff, Δw3aff, data::ProblemData)
    n = data.n
    m = data.m

    r = [- (data.E*w1 + data.F*w2 - data.b);
         - (data.G*w1 + data.H*w2 + data.J*w3 - data.c);
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

function progress(res)
    ϵ = min(0.1, res^2)
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

function initial_state(data::ProblemData)
    n = data.n
    m = data.m

    A = [data.E data.F zeros(n,m);
         data.G data.H data.J]
    wt = A' * ((A * A') \ [data.b; data.c])
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
    return w10, w20, w30
end

function mehrotra_solve(data::ProblemData; newton_iter::Int=100,
        τ=0.99995, res_tol=1e-8, verbose::Bool=false)
    n = data.n
    m = data.m

    w1, w2, w3 = initial_state(data)
    iter = 0
    success = false
    for k = 1:newton_iter
        iter += 1
        Δaff = affine_direction(w1, w2, w3, data, κ=0.0)
        Δw1aff, Δw2aff, Δw3aff = unpack(Δaff, n=n, m=m)
        αhaff, μaff = step_length(w2, w3, Δw2aff, Δw3aff)
        μ = w2'*w3 / m
        σ = (μaff / μ)^3
        Δ = corrector_direction(w1, w2, w3, σ, μ, Δw2aff, Δw3aff, data)
        Δw1, Δw2, Δw3 = unpack(Δ, n=n, m=m)
        res = residual(w1, w2, w3, data)
        τ = progress(res)
        αh = corrector_step_length(w2, w3, Δw2, Δw3; τ=τ)
        w1 = w1 + αh * Δw1
        w2 = w2 + αh * Δw2
        w3 = w3 + αh * Δw3
        res = residual(w1, w2, w3, data)
        verbose && println("res: ", scn(res))
        if res < res_tol
            success = true
            break
        end
    end
    verbose && println("iter: ", iter)
    return iter, success
end

function baseline_solve(data::ProblemData; newton_iter::Int=100,
        res_tol=1e-8, verbose::Bool=false)
    κ = res_tol

    n = data.n
    m = data.m

    w1, w2, w3 = initial_state(data)
    w1, w2, w3 = unpack(1e-3ones(n + 2m), n=n, m=m)
    iter = 0
    success = false
    for k = 1:newton_iter
        iter += 1
        Δaff = affine_direction(w1, w2, w3, data, κ=κ)
        Δw1aff, Δw2aff, Δw3aff = unpack(Δaff, n=n, m=m)
        αhaff, μaff = step_length(w2, w3, Δw2aff, Δw3aff)
        w1 = w1 + αhaff * Δw1aff
        w2 = w2 + αhaff * Δw2aff
        w3 = w3 + αhaff * Δw3aff
        res = residual(w1, w2, w3, data; κ=κ)
        verbose && println("res: ", scn(res))
        if res < res_tol
            success = true
            break
        end
    end
    verbose && println("iter: ", iter)
    return iter, success
end



function evaluate(s::Simulation, ref_traj::ContactTraj;
        algorithm::Symbol=:mehrotra_solve, newton_iter=100, res_tol=1e-8)

    κ = res_tol
    H = ref_traj.H
    iter = []
    success = 0.0
    for t = 1:H
        data = get_data(s, ref_traj, κ=κ, t=t)
        it = newton_iter
        su = false
        try
            it, su = eval(algorithm)(data, newton_iter = newton_iter, res_tol = res_tol)
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

μ_iter, σ_iter, success_rate = evaluate(s1, ref_traj1, algorithm = :mehrotra_solve, res_tol=1e-4)
μ_iter, σ_iter, success_rate = evaluate(s1, ref_traj1, algorithm = :baseline_solve, res_tol=1e-4)

μ_iter, σ_iter, success_rate = evaluate(s2, ref_traj2, algorithm = :mehrotra_solve, res_tol=1e-4)
μ_iter, σ_iter, success_rate = evaluate(s2, ref_traj2, algorithm = :baseline_solve, res_tol=1e-4)

μ_iter, σ_iter, success_rate = evaluate(s3, ref_traj3, algorithm = :mehrotra_solve, res_tol=1e-4)
μ_iter, σ_iter, success_rate = evaluate(s3, ref_traj3, algorithm = :baseline_solve, res_tol=1e-4)
