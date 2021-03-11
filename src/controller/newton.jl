
# Trajectory Optimizer
function easy_lin_trajopt(probsize, model, q1_init, q2_init, q1_ref, q2_ref, lint0;
    ρ0=1e-3,#############################################
    β=1e1,#############################################
    outer_iter=1,#############################################
    res_tol=1e-8,#############################################
    solver_outer_iter::Int=3,#############################################
    solver_inner_iter::Int=10,#############################################
    utr_ref=[zeros(probsize.nu)  for k = 1:probsize.N-1],######################################################
    qtr_ref=[zeros(probsize.nq)  for k = 1:probsize.N-1],######################################################
    γtr_ref=[zeros(probsize.nγ)  for k = 1:probsize.N-1],######################################################
    btr_ref=[zeros(2probsize.nb) for k = 1:probsize.N-1],######################################################
    utr_wrm=utr_ref,############################################################################
    qtr_wrm=qtr_ref,############################################################################
    γtr_wrm=γtr_ref,############################################################################
    btr_wrm=btr_ref,############################################################################
    Qu=fill(Diagonal(1e-1*ones(probsize.nu)), probsize.N-1),##########################
    Qq=fill(Diagonal(1e+1*ones(probsize.nq)), probsize.N-1),##########################
    Qγ=fill(Diagonal(1e-7*ones(probsize.nγ)), probsize.N-1),##########################
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


function control!(model::ContactDynamicsModel, impl::ImplicitTraj11,
    ref_traj::ContactTraj, cost::CostFunction, s_opts::NewtonOptions{T}) where {T}

    # we have an initial state q0 q1


    return nothing
end


function residual!(model::ContactDynamicsModel, impl::ImplicitTraj11,
    ref_traj::ContactTraj, cost::CostFunction, s_opts::NewtonOptions{T}) where {T}
    # we have an initial state q0 q1


    return nothing
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
    end
    return jcb
end




mutable struct Newton25{T,nq,nu,nc,nb,n1,n2,n3,Vq,Vu,Vγ,Vb,VI,VIT,Vq0,Vq0T,Vq1,Vq1T,Vu1,Vu1T}
    jcb::Any                     # Jacobian
    r::Vector{T}                 # residual
    r̄::Vector{T}                 # candidate residual
    ν::Vector{SizedVector{n1,T}} # implicit dynamics lagrange multiplier
    H::Int                       # horizon
    nd::Int                      # implicit dynamics constraint size
    nr::Int                      # size of a one-time-step block
    iq::SizedVector{nq,Int}      # configuration indices
    iu::SizedVector{nu,Int}      # control indices
    iγ::SizedVector{nc,Int}      # impact indices
    ib::SizedVector{nb,Int}      # linear friction indices
    iν::SizedVector{n1,Int}      # implicit dynamics lagrange multiplier
    iz::SizedVector{n2,Int}      # IP solver solution [q2, γ1, b1]
    iθ::SizedVector{n3,Int}      # IP solver data [q0, q1, u1]
    Qq::Vector{Vq}              # cost function views
    Qu::Vector{Vu}              # cost function views
    Qγ::Vector{Vγ}              # cost function views
    Qb::Vector{Vb}              # cost function views
    IV::Vector{VI}              # dynamics -I views
    ITV::Vector{VIT}            # dynamics -I views transposed
    q0::Vector{Vq0}             # dynamics q0 views
    q0T::Vector{Vq0T}           # dynamics q0 views transposed
    q1::Vector{Vq1}             # dynamics q1 views
    q1T::Vector{Vq1T}           # dynamics q1 views transposed
    u1::Vector{Vu1}             # dynamics u1 views
    u1T::Vector{Vu1T}           # dynamics u1 views transposed
end


function Newton25(H::Int, dim::Dimensions)
    nq = dim.q # configuration
    nu = dim.u # control
    nc = dim.c # contact
    nb = dim.b # linear friction
    nd = nq + nc + nb # implicit dynamics constraint
    nr = nq+nu+nc+nb+nd # size of a one-time-step block

    off = 0
    iq = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q2
    iu = SizedVector{nu}(off .+ (1:nu)); off += nu # index of the control u1
    iγ = SizedVector{nc}(off .+ (1:nc)); off += nc # index of the impact γ1
    ib = SizedVector{nb}(off .+ (1:nb)); off += nb # index of the linear friction b1
    iν = SizedVector{nd}(off .+ (1:nd)); off += nd # index of the dynamics lagrange multiplier ν1
    iz = vcat(iq, iγ, ib) # index of the IP solver solution [q2, γ1, b1]
    iθ = vcat(iq .- 2nr, iq .- nr, iu) # index of the IP solver data [q0, q1, u1]

    jcb = spzeros(H*nr,H*nr)
    jcbV = [view(jcb, (t-1)*nr .+ (1:nr), (t-1)*nr .+ (1:nr)) for t=1:H]
    Qq  = [view(jcb, (t-1)*nr .+ iq, (t-1)*nr .+ iq) for t=1:H]
    Qu  = [view(jcb, (t-1)*nr .+ iu, (t-1)*nr .+ iu) for t=1:H]
    Qγ  = [view(jcb, (t-1)*nr .+ iγ, (t-1)*nr .+ iγ) for t=1:H]
    Qb  = [view(jcb, (t-1)*nr .+ ib, (t-1)*nr .+ ib) for t=1:H]
    IV  = [view(jcb, (t-1)*nr .+ iz, (t-1)*nr .+ iν) for t=1:H]
    ITV = [view(jcb, (t-1)*nr .+ iν, (t-1)*nr .+ iz) for t=1:H]
    q0  = [view(jcb, (t-1)*nr .+ iν, (t-3)*nr .+ iq) for t=3:H]
    q0T = [view(jcb, (t-3)*nr .+ iq, (t-1)*nr .+ iν) for t=3:H]
    q1  = [view(jcb, (t-1)*nr .+ iν, (t-2)*nr .+ iq) for t=2:H]
    q1T = [view(jcb, (t-2)*nr .+ iq, (t-1)*nr .+ iν) for t=2:H]
    u1  = [view(jcb, (t-1)*nr .+ iν, (t-1)*nr .+ iu) for t=1:H]
    u1T = [view(jcb, (t-1)*nr .+ iu, (t-1)*nr .+ iν) for t=1:H]
    ν = [zeros(SizedVector{nd}) for t=1:H]
    r = zeros(H*nr)
    r̄ = zeros(H*nr)
    return Newton25{T,nq,nu,nc,nb,nd,nd,2nq+nu,
        eltype.((Qq, Qu, Qγ, Qb))...,
        eltype.((IV, ITV, q0, q0T, q1, q1T, u1, u1T))...,
        }(
        jcb, r, r̄, ν, H, nd, nr,
        iq, iu, iγ, ib, iν, iz, iθ,
        Qq, Qu, Qγ, Qb, IV, ITV, q0, q0T, q1, q1T, u1, u1T)
end


function sparse_zero!(spm::SparseMatrixCSC)
	n = length(spm.nzval)
	for i = 1:n
		spm.nzval[i] = 0.0
	end
	return nothing
end

function jacobian!(core::Newton25, impl::ImplicitTraj11,
    cost::CostFunction, s_opts::NewtonOptions{T}) where {T}
    # unpack
    H = core.H
    jcb = core.jcb
    jcb = sparse_zero!(jcb)

    for t=1:H
        # Cost function
        core.Qq[t] .= cost.Qq[t]
        core.Qu[t] .= cost.Qu[t]
        core.Qγ[t] .= cost.Qγ[t]
        core.Qb[t] .= cost.Qb[t]
        # Implicit dynamics
        core.IV[t][diagind(core.IV[t])]   .= - 1.0
        core.ITV[t][diagind(core.ITV[t])] .= - 1.0
        if t >=3
            core.q0[t-2]  .= impl.δq0[t]
            core.q0T[t-2] .= impl.δq0[t]'
        end
        if t >= 2
            core.q1[t-1]  .= impl.δq1[t]
            core.q1T[t-1] .= impl.δq1[t]'
        end
        core.u1[t]  .= impl.δu1[t]
        core.u1T[t] .= impl.δu1[t]'
    end
    return nothing
end

T = Float64
H = 50
h = 0.03
κ = 1e-3
model = get_model("particle")
nq = model.dim.q
nu = model.dim.u
nγ = model.dim.c
nb = model.dim.b

core0 = Newton25(H, model.dim)
impl0 = ImplicitTraj11(H, model)
cost0 = CostFunction(H, model.dim,
    Qq=fill(Diagonal(1e-1*ones(SizedVector{nq})), H),
    Qu=fill(Diagonal(1e-2*ones(SizedVector{nu})), H),
    Qγ=fill(Diagonal(1e-3*ones(SizedVector{nγ})), H),
    Qb=fill(Diagonal(1e-4*ones(SizedVector{nb})), H),
    )
n_opts = NewtonOptions()
@allocated jacobian!(core0, impl0, cost0, n_opts)
@allocated jacobian!(core0, impl0, cost0, n_opts)
@allocated jacobian!(core0, impl0, cost0, n_opts)
@allocated jacobian!(core0, impl0, cost0, n_opts)
@allocated jacobian!(core0, impl0, cost0, n_opts)
@benchmark jacobian!(core0, impl0, cost0, n_opts)
core0.jcb





plot(Gray.(Matrix((1e3.*core0.jcb))))
core0.IV[1][diagind(core0.IV[1])] .= 1.0


a = spzeros(3,3)
a[1,1] = 10
a = similar(a, T)

A = [10 20 30 ; 40 50 60 ; 70 80 90]
LinearIndices(A)[2:3, 2]
A[2:3, 2]

view(A, LinearIndices(A)[diagind(A)])
LinearIndices(A)
A[diagind(A)]

LinearIndices(A)[diagind(A)]

diagind(A)

A

# interior-point solver options
@with_kw mutable struct NewtonOptions{T}
    r_tol::T = 1.0e-5            # primal dual residual tolerance
    solver_outer_iter::Int = 3   # outer iter on the κ parameter
    solver_inner_iter::Int = 10  # primal dual iter
    κ_init::T = 1e-3             # inner solver intialization
    κ_scale::T = 0.1             # inner solver scaling
    κ_tol::T = 2.0e-3            # inner solver tolerance
    β::T = 1e1                   # dual regularization
end




T = Float64
H = 10
h = 0.03
κ = 1e-3
model = get_model("particle")
nq = model.dim.q
nu = model.dim.u

q0 = SVector{nq,T}([0.0, 0.0, 0.2])
q1 = SVector{nq,T}([0.1, 0.1, 0.2])


ip_opts = InteriorPointOptions(κ_init=κ, κ_tol=κ*2, r_tol=1e-5)
sim0 = simulator2_base(model, q0, q1, h, H;
    u = [@SVector zeros(model.dim.u) for t = 1:H],
    w = [@SVector zeros(model.dim.w) for t = 1:H],
    ip_opts = ip_opts,
    sim_opts = SimulatorOptions{T}())

simulate!(sim0; verbose = false)
ref_traj0 = deepcopy(sim0.traj)
cost0 = CostFunction(H, model.dim,
    Qq=fill(Diagonal(1e-1*ones(SizedVector{nq})), H),
    Qu=fill(Diagonal(1e-1*ones(SizedVector{nu})), H),
    )
n_opts = NewtonOptions()
impl = ImplicitTraj11(H, model)

control!(model, impl, ref_traj0, cost1, n_opts)
