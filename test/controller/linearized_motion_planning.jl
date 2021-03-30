@testset "Controller: Linearized Motion Planning" begin
    model = get_model("quadruped")
    κ = 1.0e-4
    ref_traj = ContactControl.get_trajectory("quadruped", "gait1")
    ref_traj.κ .= κ
    H = ref_traj.H
    h = ref_traj.h
    nq = model.dim.q
    nu = model.dim.u
    nw = model.dim.w
    nc = model.dim.c
    nb = model.dim.b
    nd = nq + nc + nb
    nr = nq + nu + nc + nb + nd

    # Test Jacobian!
    cost = ContactControl.CostFunction(H, model.dim,
        q = [Diagonal(1.0e-2 *
            ([0.02, 0.02, 1.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15])) for t = 1:H],
        u = [Diagonal(3.0e-2 * ones(nu)) for t = 1:H],
        γ = [Diagonal(1.0e-6 * ones(nc)) for t = 1:H],
        b = [Diagonal(1.0e-6 * ones(nb)) for t = 1:H])

    im_traj = ContactControl.ImplicitTraj(ref_traj, model)

    num_pr_mp = H * (nq + nu + nc + nb)
    num_du_mp = H * nd
    num_var_mp = num_pr_mp + num_du_mp
    num_data_mp = 0

    τ0 = vcat([[ref_traj.q[t+2]; ref_traj.u[t]; ref_traj.γ[t]; ref_traj.b[t]] for t = 1:H]...)
    ν0 = zeros(H * nd)
    x0 = [τ0; ν0]

    opts = ContactControl.InteriorPointOptions(
            κ_init = 1.0,
            κ_tol = 1.0,
            r_tol = 1.0e-3,
            max_iter_outer = 1,
            res_norm = 2,
            reg = false,
            diff_sol = false)

    ip = ContactControl.interior_point(copy(x0), zeros(num_data_mp),
             idx_ineq = collect(1:0),
             idx_pr = collect(1:num_pr_mp),
             idx_du = collect(num_pr_mp .+ (1:num_du_mp)),
             opts = opts)

    mp_traj = ContactControl.trajectory_x(model, ip.z, ref_traj.q[1], ref_traj.q[2], H, h, ref_traj.κ)

    @test maximum([norm(mp_traj.q[t] - ref_traj.q[t]) for t = 1:length(ref_traj.q)]) < 1.0e-8
    @test maximum([norm(mp_traj.u[t] - ref_traj.u[t]) for t = 1:length(ref_traj.u)]) < 1.0e-8
    @test maximum([norm(mp_traj.γ[t] - ref_traj.γ[t]) for t = 1:length(ref_traj.γ)]) < 1.0e-8
    @test maximum([norm(mp_traj.b[t] - ref_traj.b[t]) for t = 1:length(ref_traj.b)]) < 1.0e-8

    res = ContactControl.LMPResidual(ip.r, H, model.dim)
    res_cand = ContactControl.LMPResidual(ip.r̄, H, model.dim)
    jac = ContactControl.LMPJacobian(ip.rz, H, model.dim)

    r_cache = ContactControl.LMPResidualCache(res, cost, mp_traj, ref_traj, im_traj, model)
    r̄_cache = ContactControl.LMPResidualCache(res_cand, cost, mp_traj, ref_traj, im_traj, model)
    rz_cache = ContactControl.LMPJacobianCache(jac, cost, im_traj)

    ip.methods.r! = ContactControl.r_lmp!
    ip.methods.rz! = ContactControl.rz_lmp!
    ip.r_cache = r_cache
    ip.r̄_cache = r̄_cache
    ip.rz_cache = rz_cache

    ip.z .= copy([τ0; ν0])
    ip.methods.r!(ip.r, ip.z, ip.θ, ip.κ[1], ip.r_cache)
    ip.methods.rz!(ip.rz, ip.z, ip.θ, ip.rz_cache)

    @test norm(ip.r_cache.traj.x - [τ0; ν0]) < 1.0e-8
    ContactControl.implicit_dynamics!(im_traj, model, mp_traj)
    @test norm(mp_traj.θ[1][1:2nq + nu + nw + 1] - ref_traj.θ[1][1:2nq + nu + nw + 1]) < 1.0e-8

    @test norm(r_cache.res.r - ip.r) < 1.0e-8
    ContactControl.r_lmp!(ip.r̄, ip.z, ip.θ, ip.κ[1], ip.r̄_cache)
    @test norm(r̄_cache.res.r - ip.r̄) < 1.0e-8

    ContactControl.rz_lmp!(ip.rz, ip.z, ip.θ, ip.rz_cache)
    @test norm(ip.rz_cache.jac.R - ip.rz) < 1.0e-8
    sz = size(ip.rz)
    @test rank(Array(ip.rz_cache.jac.R)) == sz[1]
    @test rank(ip.rz) == sz[1]

    ContactControl.r_lmp!(ip.r, ip.z, ip.θ, ip.κ[1], ip.r_cache)
    ContactControl.r_lmp!(ip.r̄, ip.z, ip.θ, ip.κ[1], ip.r̄_cache)
    ContactControl.rz_lmp!(ip.rz, ip.z, ip.θ, ip.rz_cache)
    norm(ip.r)
    norm(ip.r̄)
    ip.Δ .= ip.rz \ ip.r
    α = 1.0
    ip.z̄ .= ip.z - α * ip.Δ
    ContactControl.r_lmp!(ip.r̄, ip.z̄, ip.θ, ip.κ[1], ip.r̄_cache)
    norm(ip.r̄)
    ip.z .= ip.z̄

    ip.z .= copy([τ0; ν0])
    status = ContactControl.interior_point!(ip)
    @test status

    ip.z .= copy([τ0; ν0])
    ip.z[1:2] .+= 1.0e-2 * rand(2)
    status = ContactControl.interior_point!(ip)
    @test status

    solver = ContactControl.linear_motion_planning_solver(model, ref_traj, cost)
    status = ContactControl.lmp!(solver)
    @test status
end
