function indices_z(s::Simulation) 
    model = s.model
    nq = model.nq 
    nc = model.nc
    q = collect(1:nq) 
    γ = collect(index_γ1(s.model, s.env)) 
    sγ = collect(1:0)
    ψ = collect(1:0) 
    b = collect(index_b1(s.model, s.env)) 
    sψ = collect(1:0) 
    sb = collect(1:0) 
    IndicesZ(q, γ, sγ, ψ, b, sψ, sb)
end

function initialize_z!(z, model, idx::IndicesZ, q)
    z .= 1.0
    z[idx.q] .= q
end

function simulator(s, T; 
    h=0.01,
    policy=empty_policy(s.model), 
    dist=empty_disturbances(s.model), 
    f=friction_coefficients(s.model),
    residual=s.res.r!, 
    jacobian_z=s.res.rz!, 
    jacobian_θ=s.res.rθ!,
    diff_sol=false,
    solver_opts=InteriorPointOptions(
        undercut=Inf,
        γ_reg=0.1,
        r_tol=1e-8,
        κ_tol=1e-8,  
        max_ls=25,
        ϵ_min=0.25,
        diff_sol=diff_sol,
        verbose=false),
    stats=SimulatorStatistics(T),
    sim_opts=SimulatorOptions()
    )
     
    idx_z = indices_z(s)
    idx_θ = indices_θ(s.model, nf=length(f))
    idx_opt = IndicesOptimization(s.model, s.env)
    nz = num_var(s.model, s.env) 
    nθ = num_data(s.model)
    z0 = zeros(nz) 
    θ0 = zeros(nθ) 
    q0 = zeros(s.model.nq)
    initialize_z!(z0, s.model, idx_z, q0)
    initialize_θ!(θ0, s.model, idx_θ, q0, q0, zeros(s.model.nu), zeros(s.model.nw), f, h)
     
    ip = interior_point(
         z0,
         θ0,
         idx=idx_opt,
         r! = residual,
         rz! = jacobian_z,
         rθ! = jacobian_θ,
         rz=zeros(nz, nz),
         rθ=zeros(nz, nθ),
         opts=solver_opts)

    traj = Trajectory(s.model, T, nb=s.model.nc * friction_dim(s.env))
    grad = GradientTrajectory(s.model, T, nb=s.model.nc * friction_dim(s.env))
    
    return Simulator(s.model, policy, dist, traj, grad, ip, idx_z, idx_θ, f, h, stats, sim_opts)
end

"""
    control saturation
"""
control_saturation(u, uL, uU) = min.(max.(uL, u), uU)
