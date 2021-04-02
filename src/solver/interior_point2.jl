

function interior_point2(x, θ;
        num_var = length(x),
        num_data = length(θ),
        idx_ineq = collect(1:0),
        idx_pr = collect(1:num_var),
        idx_du = collect(1:0),
        r! = r!, rz! = rz!, rθ! = rθ!,
        # rz = spzeros(num_var, num_var),
        # rθ = spzeros(num_var, num_data),
        r  = zeros(num_var), #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        rz = zeros(num_var, num_var), #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        rθ = zeros(num_var, num_data), #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        reg_pr = 0.0, reg_du = 0.0,
        r_cache = NoCache(),
        r̄_cache = NoCache(),
        rz_cache = NoCache(),
        rθ_cache = NoCache(),
        opts = InteriorPointOptions()) where T

    # rz_!(rz, x, θ, rz_cache) # compute Jacobian for pre-factorization

    InteriorPoint(
        ResidualMethods(r!, rz!, rθ!),
        x,
        zeros(num_var),
        r,# zeros(num_var),#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        0.0,
        deepcopy(r),# zeros(num_var), #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        0.0,
        rz,
        rθ,
        zeros(num_var),
        idx_ineq,
        idx_pr,
        idx_du,
        zeros(num_var, num_data),
        θ,
        zeros(1),
        num_var,
        num_data,
        eval(opts.solver)(zeros(1,1)),
        view(zeros(num_var,num_var), CartesianIndex.(idx_pr, idx_pr)),#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        view(zeros(num_var,num_var), CartesianIndex.(idx_du, idx_du)),#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        reg_pr, reg_du,
        r_cache, r̄_cache, rz_cache, rθ_cache,
        opts)
end

# interior point solver
function interior_point2!(ip::InteriorPoint{T}) where T

    # # methods
    # r! = ip.methods.r!#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # rz! = ip.methods.rz!#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # rθ! = ip.methods.rθ!#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # options
    opts = ip.opts
    r_tol = opts.r_tol
    κ_tol = opts.κ_tol
    κ_init = opts.κ_init
    κ_scale = opts.κ_scale
    ls_scale = opts.ls_scale
    max_iter_inner = opts.max_iter_inner
    max_iter_outer = opts.max_iter_outer
    max_time = opts.max_time
    max_ls = opts.max_ls
    diff_sol = opts.diff_sol
    res_norm = opts.res_norm
    reg = opts.reg

    # unpack pre-allocated data
    z = ip.z
    z̄ = ip.z̄
    r = ip.r
    r_norm = ip.r_norm
    r̄ = ip.r̄
    r̄_norm = ip.r̄_norm
    rz = ip.rz
    Δ = ip.Δ
    idx_ineq = ip.idx_ineq
    θ = ip.θ
    κ = ip.κ
    v_pr = ip.v_pr
    v_du = ip.v_du
    reg_pr = ip.reg_pr
    reg_du = ip.reg_du
    r_cache = ip.r_cache
    r̄_cache = ip.r̄_cache
    rz_cache = ip.rz_cache
    solver = ip.solver

    # initialize barrier parameter
    κ[1] = κ_init

    # compute residual, residual Jacobian
    # r!(r, z, θ, κ[1], r_cache) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    r!(r, z, θ, κ[1])
    r_norm = norm(r, res_norm)

    elapsed_time = 0.0

    for i = 1:max_iter_outer
        elapsed_time >= max_time && break
        for j = 1:max_iter_inner
            elapsed_time >= max_time && break
            elapsed_time += @elapsed begin
                # check for converged residual
                if r_norm < r_tol
                    break
                end

                # compute residual Jacobian
                # rz!(rz, z, θ, rz_cache) ######$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                rz!(rz, z)

                # regularize (fixed, TODO: adaptive)
                # reg && regularize!(v_pr, v_du, reg_pr, reg_du) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

                # compute step
                # linear_solve!(solver, Δ, rz, r)#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                linear_solve!(Δ, rz, r)

                # initialize step length
                α = 1.0

                # candidate point
                z̄ .= z - α * Δ

                # check inequality constraints
                iter = 0
                while inequality_check(z̄, idx_ineq)
                    α *= ls_scale
                    z̄ .= z - α * Δ
                    iter += 1
                    if iter > max_ls
                        @error "backtracking line search fail"
                        return false
                    end
                end

                # reduce norm of residual
                # r!(r̄, z̄, θ, κ[1], r̄_cache) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                r!(r̄, z̄, θ, κ[1])
                r̄_norm = norm(r̄, res_norm)

                while r̄_norm >= (1.0 - 0.001 * α) * r_norm
                    α *= ls_scale
                    z̄ .= z - α * Δ

                    # r!(r̄, z̄, θ, κ[1], r̄_cache) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    r!(r̄, z̄, θ, κ[1])
                    r̄_norm = norm(r̄, Inf)

                    iter += 1
                    if iter > max_ls
                        @error "line search fail"
                        return false
                    end
                end

                # update
                z .= z̄
                # r .= r̄#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                r!(r, z, θ, κ[1])
                r_norm = r̄_norm
            end
        end

        if κ[1] <= κ_tol
            # differentiate solution
            diff_sol && differentiate_solution2!(ip)
            return true
        else
            # update barrier parameter
            κ[1] *= κ_scale

            # update residual
            # r!(r, z, θ, κ[1], r_cache) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            r!(r, z, θ, κ[1])
            r_norm = norm(r, res_norm)
        end
    end
end



function differentiate_solution2!(ip::InteriorPoint)
    z = ip.z
    θ = ip.θ
    rz = ip.rz
    rθ = ip.rθ
    δz = ip.δz
    rz_cache = ip.rz_cache
    rθ_cache = ip.rθ_cache
    κ = ip.κ

    # ip.methods.rz!(rz, z, θ, rz_cache) # maybe not needed #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # ip.methods.rθ!(rθ, z, θ, rθ_cache) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    rz!(rz, z)

    # linear_matrix_solve!(ip.solver, δz, rz, rθ) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    linear_solve!(δz, rz, rθ)
    @inbounds @views @. ip.δz .*= -1.0
    nothing
end
