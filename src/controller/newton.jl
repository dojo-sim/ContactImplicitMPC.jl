
# Trajectory Optimizer
function easy_lin_trajopt(probsize, model, q1_init, q2_init, q1_ref, q2_ref, lint0;
    ρ0=1e-3,
    β=1e1,
    outer_iter=1,
    res_tol=1e-8,
    solver_outer_iter::Int=3,
    solver_inner_iter::Int=10,
    utr_ref=[zeros(probsize.nu)  for k = 1:probsize.N-1],######################################################
    qtr_ref=[zeros(probsize.nq)  for k = 1:probsize.N-1],######################################################
    γtr_ref=[zeros(probsize.nγ)  for k = 1:probsize.N-1],######################################################
    btr_ref=[zeros(2probsize.nb) for k = 1:probsize.N-1],######################################################
    utr_wrm=utr_ref,######################################################
    qtr_wrm=qtr_ref,######################################################
    γtr_wrm=γtr_ref,######################################################
    btr_wrm=btr_ref,######################################################
    Qu=fill(Diagonal(1e-1*ones(probsize.nu)), probsize.N-1),#########################
    Qq=fill(Diagonal(1e+1*ones(probsize.nq)), probsize.N-1),#########################
    Qγ=fill(Diagonal(1e-7*ones(probsize.nγ)), probsize.N-1),#########################
    Qb=fill(Diagonal(1e-7*ones(2probsize.nb)), probsize.N-1),#########################
    u_amp=1e-1,
    live_plot::Bool=false,
    z_init=z_init,
    )
    N = probsize.N
    nq = probsize.nq
    nu = probsize.nu
    nγ = probsize.nγ
    nb = probsize.nb
    nz = nq+4nγ+4nb # residual of the step and variables of the step
    nθ = 2nq+nu # parameters of the step
    nλ = nq+nγ+2nb # dynamics lagrange multiplier
    m = nz+nu+nλ # primal dual vars in a step
    M = (N-1)*m # primal dual vars
    nr = 2nq+nu+2nγ+4nb # outer problem resiudal size

    traj = zeros(M)
    # Initilaization
    set_var_traj!(N,model, traj, utr_wrm, :u)
    set_var_traj!(N,model, traj, qtr_wrm, :q)
    set_var_traj!(N,model, traj, γtr_wrm, :γ)
    set_var_traj!(N,model, traj, btr_wrm, :b)

    function res(traj)
        r = zeros(eltype(traj), (N-1)*nr)
        rtr = [view(r, (k-1)*nr .+ (1:nr)) for k=1:N-1]
        # we organise the residual as follows [u,q,γ,b,constraint]
        indu = Vector(1:nu)
        indq = Vector(nu .+ (1:nq))
        indγ = Vector(nu+nq .+ (1:nγ))
        indb = Vector(nu+nq+nγ .+ (1:2nb))
        indc = Vector(nu+nq+nγ+2nb .+ (1:nλ))
        # indθ = vcat(indq-) # q1,q2,u2
        inds = vcat(indq, indγ, indb) # q3,γ2,b2

        utr = get_var_traj(N,model,traj,:u)
        qtr = get_var_traj(N,model,traj,:q)
        γtr = get_var_traj(N,model,traj,:γ)
        btr = get_var_traj(N,model,traj,:b)
        λdtr = get_var_traj(N,model,traj,:λd)
        for k = 1:N-1

            # optimality and constraints
            q1, q2 = get_init_config(N, [q1_init, q2_init], traj, k)
            zk, ∇k, flag = easy_lin_step(lint0[k],
                q1, q2, model, u2_ref=utr[k], ρ0=ρ0, outer_iter=outer_iter, z_init=z_init)
            ck = zk[1:nλ] - [qtr[k]; γtr[k]; btr[k]]

            rtr[k][indu] .+= Qu[k]*(utr[k] - utr_ref[k])
            rtr[k][indq] .+= Qq[k]*(qtr[k] - qtr_ref[k])
            rtr[k][indγ] .+= Qγ[k]*(γtr[k] - γtr_ref[k])
            rtr[k][indb] .+= Qb[k]*(btr[k] - btr_ref[k])
            rtr[k][indc] .+= ck

            # Constraint derivative
            # Minus Identity term #∇qk1, ∇γk, ∇bk
            rtr[k][nu .+ (1:nλ)] += -λdtr[k]
            # Implicit function theorem part #∇qk_1, ∇qk, ∇uk
            off = 0
            ∇q1 = ∇k[1:nλ, off .+ (1:nq)] # ∈ nλ×nθ
            off += nq
            ∇q2 = ∇k[1:nλ, off .+ (1:nq)] # ∈ nλ×nθ
            off += nq
            ∇u2 = ∇k[1:nλ, off .+ (1:nu)] # ∈ nλ×nθ
            off += nu

            k >= 3 ? rtr[k-2][indq] += ∇q1'*λdtr[k] : nothing
            k >= 2 ? rtr[k-1][indq] += ∇q2'*λdtr[k] : nothing
            rtr[k][indu] += ∇u2'*λdtr[k]

            # # AL terms
            # if k >= 3
            #     rtr[k][inds] += 0.5/ρ0 * -I * ck
            #     k >= 3 ? rtr[k-2][indq] += 0.5/ρ0 * ∇q1'*ck : nothing
            #     k >= 2 ? rtr[k-1][indq] += 0.5/ρ0 * ∇q2'*ck : nothing
            #     rtr[k][indu] += 0.5/ρ0 * ∇u2'*ck
            # end
        end
        return r
    end

    function jac(traj)
        jcb = spzeros((N-1)*nr,(N-1)*nr)
        indu = Vector(1:nu)
        indq = Vector(nu .+ (1:nq))
        indγ = Vector(nu+nq .+ (1:nγ))
        indb = Vector(nu+nq+nγ .+ (1:2nb))
        indc = Vector(nu+nq+nγ+2nb .+ (1:nλ))
        # indθ = vcat(indq-) # q1,q2,u2
        inds = vcat(indq, indγ, indb) # q3,γ2,b2
        indd = vcat(indq .- 2nr, indq .- nr, indu, indq, indγ, indb) # q1,q2,u2,q3,γ2,b2

        utr = get_var_traj(N,model,traj,:u)
        qtr = get_var_traj(N,model,traj,:q)
        γtr = get_var_traj(N,model,traj,:γ)
        btr = get_var_traj(N,model,traj,:b)
        λdtr = get_var_traj(N,model,traj,:λd)

        jcbtr = [view(jcb, (k-1)*nr .+ (1:nr), (k-1)*nr .+ (1:nr)) for k=1:N-1]
        for k = 1:N-1
            jcbtr[k][indu,indu] .+= Qu[k]
            jcbtr[k][indq,indq] .+= Qq[k]
            jcbtr[k][indγ,indγ] .+= Qγ[k]
            jcbtr[k][indb,indb] .+= Qb[k]
            jcbtr[k][indc,inds] += -I
            jcbtr[k][inds,indc] += -I

            # Linearization
            q1, q2 = get_init_config(N, [q1_init, q2_init], traj, k)
            zk, ∇k, flag = easy_lin_step(lint0[k],
                q1, q2, model, u2_ref=utr[k], ρ0=ρ0, outer_iter=outer_iter, z_init=z_init)

            # Implicit function theorem part #∇qk_1, ∇qk, ∇uk
            off = 0
            ∇q1 = ∇k[1:nλ, off .+ (1:nq)] # ∈ nλ×nθ
            off += nq
            ∇q2 = ∇k[1:nλ, off .+ (1:nq)] # ∈ nλ×nθ
            off += nq
            ∇u2 = ∇k[1:nλ, off .+ (1:nu)] # ∈ nλ×nθ
            off += nu

            ind_k_2 = (k-3)*nr .+ (1:nr)
            ind_k_1 = (k-2)*nr .+ (1:nr)
            ind_k   = (k-1)*nr .+ (1:nr)

            k >= 3 ? jcb_2  = view(jcb, ind_k, ind_k_2) : nothing
            k >= 3 ? jcb_2T = view(jcb, ind_k_2, ind_k) : nothing
            k >= 2 ? jcb_1  = view(jcb, ind_k, ind_k_1) : nothing
            k >= 2 ? jcb_1T = view(jcb, ind_k_1, ind_k) : nothing
            jcb_0 = view(jcb, ind_k, ind_k)

            k >= 3 ? jcb_2[indc, indq] .+= ∇q1 : nothing
            k >= 2 ? jcb_1[indc, indq] .+= ∇q2 : nothing
            jcb_0[indc, indu] .+= ∇u2

            k >= 3 ? jcb_2T[indq, indc] .+= ∇q1' : nothing
            k >= 2 ? jcb_1T[indq, indc] .+= ∇q2' : nothing
            jcb_0[indu, indc] .+= ∇u2'

            # Regularization
            jcbtr[k][indc, indc] += -β*ρ0*I

            # AL terms
            # ∇c = sparse([∇k[1:nλ, 1:nθ] -I(nλ)])
            # # @show size(∇c)
            # # @show size(∇c')
            # # @show size(∇c'*∇c)
            # k >= 3 ? (@show size(jcb[(k-1)*nr .+ indd, (k-1)*nr .+ indd] )) : nothing
            # k >= 3 ? jcb[(k-1)*nr .+ indd, (k-1)*nr .+ indd] += 1/ρ0 * ∇c'*∇c : nothing
        end
        return jcb
    end

    function easy_jac(traj; ϵ=1e-5)
        M = length(traj)
        indu = get_var_ind(N,model,traj,:u)
        indq = get_var_ind(N,model,traj,:q)
        indγ = get_var_ind(N,model,traj,:γ)
        indb = get_var_ind(N,model,traj,:b)
        indλd = get_var_ind(N,model,traj,:λd)
        masktr = [[indu[k]; indq[k]; indγ[k]; indb[k]; indλd[k]] for k=1:N-1]
        mask = vcat(masktr...)

        Nr = (N-1)*nr
        jcb = zeros(Nr, Nr)
        ep = zeros(M)
        em = zeros(M)
        for (i,j) in enumerate(mask)
            @show i/length(mask)
            ep[j] += ϵ
            em[j] -= ϵ
            jcb[:,i] = (res(traj+ep) - res(traj+em))./(2ϵ)
            ep *= 0.0
            em *= 0.0
        end
        return jcb
    end

    function solver(traj)
        indu = get_var_ind(N,model,traj,:u)
        indq = get_var_ind(N,model,traj,:q)
        indγ = get_var_ind(N,model,traj,:γ)
        indb = get_var_ind(N,model,traj,:b)
        indλd = get_var_ind(N,model,traj,:λd)
        masktr = [[indu[k]; indq[k]; indγ[k]; indb[k]; indλd[k]] for k=1:N-1]
        mask = vcat(masktr...)
        unmask = setdiff(1:length(traj), mask)

        traj_back = deepcopy(traj)
        traj_back[mask] .= 0.0

        for i = 1:solver_outer_iter
            live_plot ? visualize!(vis, model, get_var_traj(N,model,traj,:q), Δt=model.dt) : nothing
            # sleep(1.0)
            for l = 1:solver_inner_iter
                j = jac(traj)
                r = res(traj)
                Δ = j \ r
                if norm(r,1)/length(r) < res_tol
                    break
                end


                # line search the step direction
                α = 1.0
                iter = 0
                traj_trial = deepcopy(traj)
                traj_trial[mask] -= α * Δ

                # plt = plot()
                # plot!(reshape(vcat(get_var_traj(N,model,traj,:u)...), (nu,N-1))', label="u")
                # display(plt)

                while norm(res(traj_trial))^2.0 >= (1.0 - 0.001 * α) * norm(r)^2.0
                    α = 0.5 * α
                    # println("   α = $α")
                    iter += 1
                    if iter > 6
                        break
                    end
                    # if iter > 50
                    #     @error "line search fail"
                    #     flag = false
                    #     return traj, false
                    # end
                    traj_trial = deepcopy(traj)
                    traj_trial[mask] -= α * Δ
                end

                # update
                if iter > 6
                    β = min(β*1.3, 1e2)
                else
                    β = max(1e1, β/1.3)
                end
                println("ρ0: ", scn(ρ0, digits=0) ,
                    "     r: ", scn(norm(r,1)/length(r), digits=0),
                    "     Δ: ", scn(norm(Δ,1)/length(Δ), digits=0),
                    "     α: ", -Int(round(log(α))))
                traj[mask] -= α * Δ
            end
            ρ0 /= 10.0
        end
        return traj, flag
    end

    set_var_traj!(N,model, traj, [zeros(nλ) for k=1:N-1], :λd) ##?????????????????????? why
    traj, flag = solver(traj)
    return traj#, res(traj), jac(traj)#, easy_jac(traj)
end


function control!(ref_traj::ContactTraj, cost::CostFunction)

    return nothing
end


model = get_model("particle")
H = 10
nq = model.dim.q
nu = model.dim.u
cost0 = CostFunction(H, model.dim)
cost1 = CostFunction(H, model.dim,
    Qq=fill(Diagonal(1e-1*ones(SizedVector{nq})), H),
    Qu=fill(Diagonal(1e-1*ones(SizedVector{nu})), H),
    )

control!(cost1)
