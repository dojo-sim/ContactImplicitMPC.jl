# Newton solver options
@with_kw mutable struct NewtonOptions{T}
    r_tol::T = 1.0e-5            # primal dual residual tolerance
    solver_outer_iter::Int = 3   # outer iter on the κ parameter
    solver_inner_iter::Int = 10  # primal dual iter
    κ_init::T = 1e-3             # inner solver intialization
    κ_scale::T = 0.1             # inner solver scaling
    κ_tol::T = 2.0e-3            # inner solver tolerance
    β_init::T = 1e1              # initial dual regularization
    β::T = 1e1                   # dual regularization
    live_plots::Bool = false     # visualize the trajectory during the solve
end

mutable struct Residual{T,
    vq2,vu1,vγ1,vb1,
    vd,vI,vq0,vq1,
    }
    r::Vector{T}                           # residual

    q2::Vector{vq2}                        # rsd cost function views
    u1::Vector{vu1}                        # rsd cost function views
    γ1::Vector{vγ1}                        # rsd cost function views
    b1::Vector{vb1}                        # rsd cost function views

    rd::Vector{vd}                         # rsd dynamics lagrange multiplier views
    rI::Vector{vI}                         # rsd dynamics -I views [q2, γ1, b1]
    q0::Vector{vq0}                        # rsd dynamics q0 views
    q1::Vector{vq1}                        # rsd dynamics q1 views
 end

function Residual(H::Int, dim::Dimensions)
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

    r = zeros(H*nr)

    q2  = [view(r, (t-1)*nr .+ iq) for t=1:H]
    u1  = [view(r, (t-1)*nr .+ iu) for t=1:H]
    γ1  = [view(r, (t-1)*nr .+ iγ) for t=1:H]
    b1  = [view(r, (t-1)*nr .+ ib) for t=1:H]

    rd  = [view(r, (t-1)*nr .+ iν) for t=1:H]
    rI  = [view(r, (t-1)*nr .+ iz) for t=1:H]
    q0  = [view(r, (t-3)*nr .+ iq) for t=3:H]
    q1  = [view(r, (t-2)*nr .+ iq) for t=2:H]
    T = eltype(r)
    return Residual{T,
        eltype.((q2, u1, γ1, b1))...,
        eltype.((rd, rI, q0, q1))...,
        }(
        r,
        q2, u1, γ1, b1,
        rd, rI, q0, q1,
        )
end

mutable struct Jacobian{T,
    Vq,Vu,Vγ,Vb,
    VI,VIT,Vq0,Vq0T,Vq1,Vq1T,Vu1,Vu1T,Vreg,
    }
    j::SparseMatrixCSC{T,Int}                 # jacobian

    Qq2::Vector{Vq}                          # jcb cost function views
    Qu1::Vector{Vu}                          # jcb cost function views
    Qγ1::Vector{Vγ}                          # jcb cost function views
    Qb1::Vector{Vb}                          # jcb cost function views

    IV::Vector{VI}                          # jcb dynamics -I views [q2, γ1, b1]
    ITV::Vector{VIT}                        # jcb dynamics -I views [q2, γ1, b1] transposed
    q0::Vector{Vq0}                         # jcb dynamics q0 views
    q0T::Vector{Vq0T}                       # jcb dynamics q0 views transposed
    q1::Vector{Vq1}                         # jcb dynamics q1 views
    q1T::Vector{Vq1T}                       # jcb dynamics q1 views transposed
    u1::Vector{Vu1}                         # jcb dynamics u1 views
    u1T::Vector{Vu1T}                       # jcb dynamics u1 views transposed
    reg::Vector{Vreg}                       # jcb dual regularization views
end

function Jacobian(H::Int, dim::Dimensions)
    nq = dim.q # configuration
    nu = dim.u # control
    nw = dim.w # disturbance
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

    j = spzeros(H*nr,H*nr)

    Qq2  = [view(j, (t-1)*nr .+ iq, (t-1)*nr .+ iq) for t=1:H]
    Qu1  = [view(j, (t-1)*nr .+ iu, (t-1)*nr .+ iu) for t=1:H]
    Qγ1  = [view(j, (t-1)*nr .+ iγ, (t-1)*nr .+ iγ) for t=1:H]
    Qb1  = [view(j, (t-1)*nr .+ ib, (t-1)*nr .+ ib) for t=1:H]

    IV  = [view(j, (t-1)*nr .+ iz, (t-1)*nr .+ iν) for t=1:H]
    ITV = [view(j, (t-1)*nr .+ iν, (t-1)*nr .+ iz) for t=1:H]
    q0  = [view(j, (t-1)*nr .+ iν, (t-3)*nr .+ iq) for t=3:H]
    q0T = [view(j, (t-3)*nr .+ iq, (t-1)*nr .+ iν) for t=3:H]
    q1  = [view(j, (t-1)*nr .+ iν, (t-2)*nr .+ iq) for t=2:H]
    q1T = [view(j, (t-2)*nr .+ iq, (t-1)*nr .+ iν) for t=2:H]
    u1  = [view(j, (t-1)*nr .+ iν, (t-1)*nr .+ iu) for t=1:H]
    u1T = [view(j, (t-1)*nr .+ iu, (t-1)*nr .+ iν) for t=1:H]
    reg = [view(j, (t-1)*nr .+ iν, (t-1)*nr .+ iν) for t=1:H]

    return Jacobian{eltype(j),
        eltype.((Qq2, Qu1, Qγ1, Qb1))...,
        eltype.((IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg))...,
        }(
        j,
        Qq2, Qu1, Qγ1, Qb1,
        IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg,
        )
end

mutable struct CoreIndex{nq,nu,nc,nb,n1,n2,n3,I<:Int}
    H::I                                    # horizon
    nd::I                                   # implicit dynamics constraint
    nr::I                                   # size of a one-time-step block
    iq::SizedArray{Tuple{nq},I,1,1}         # configuration indices
    iu::SizedArray{Tuple{nu},I,1,1}         # control indices
    iγ::SizedArray{Tuple{nc},I,1,1}         # impact indices
    ib::SizedArray{Tuple{nb},I,1,1}         # linear friction indices
    iν::SizedArray{Tuple{n1},I,1,1}         # implicit dynamics lagrange multiplier
    iz::SizedArray{Tuple{n2},I,1,1}         # IP solver solution [q2, γ1, b1]
    iθ::SizedArray{Tuple{n3},I,1,1}         # IP solver data [q0, q1, u1]
    Iq::Vector{SizedArray{Tuple{nq},I,1,1}} # configuration indices
    Iu::Vector{SizedArray{Tuple{nu},I,1,1}} # control indices
    Iγ::Vector{SizedArray{Tuple{nc},I,1,1}} # impact indices
    Ib::Vector{SizedArray{Tuple{nb},I,1,1}} # linear friction indices
    Iν::Vector{SizedArray{Tuple{n1},I,1,1}} # implicit dynamics lagrange multiplier
    Iz::Vector{SizedArray{Tuple{n2},I,1,1}} # IP solver solution [q2, γ1, b1]
    Iθ::Vector{SizedArray{Tuple{n3},I,1,1}} # IP solver data [q0, q1, u1]
end

function CoreIndex(H::Int, dim::Dimensions)
    nq = dim.q # configuration
    nu = dim.u # control
    nw = dim.w # disturbance
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

    Iq = [(t-1)*nr .+ iq for t=1:H]
    Iu = [(t-1)*nr .+ iu for t=1:H]
    Iγ = [(t-1)*nr .+ iγ for t=1:H]
    Ib = [(t-1)*nr .+ ib for t=1:H]
    Iν = [(t-1)*nr .+ iν for t=1:H]
    Iz = [(t-1)*nr .+ iz for t=1:H]
    Iθ = [(t-1)*nr .+ iθ for t=1:H]

    return CoreIndex{nq,nu,nc,nb,nd,nd,2nq+nu,Int}(
        H, nd, nr,
        iq, iu, iγ, ib, iν, iz, iθ,
        Iq, Iu, Iγ, Ib, Iν, Iz, Iθ,
        )
end

mutable struct Newton{T,nq,nu,nw,nc,nb,nz,nθ,n1,n2,n3,
    }
    j::Jacobian{T}                                  # Jacobian
    r::Residual{T}                                  # residual
    r̄::Residual{T}                                  # candidate residual
    Δ::Residual{T}                                  # step direction in the Newton solve, it contains: q2-qH+1, u1-uH, γ1-γH, b1-bH, λd1-λdH
    ν::Vector{SizedArray{Tuple{n1},T,1,1}}          # implicit dynamics lagrange multiplier
    ν_::Vector{SizedArray{Tuple{n1},T,1,1}}         # candidate implicit dynamics lagrange multiplier
    traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}       # optimized trajectory
    trial_traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ} # trial trajectory used in line search
    Δq::Vector{SizedArray{Tuple{nq},T,1,1}}         # difference between the traj and ref_traj
    Δu::Vector{SizedArray{Tuple{nu},T,1,1}}         # difference between the traj and ref_traj
    Δγ::Vector{SizedArray{Tuple{nc},T,1,1}}         # difference between the traj and ref_traj
    Δb::Vector{SizedArray{Tuple{nb},T,1,1}}         # difference between the traj and ref_traj
    H::Int                                          # horizon
    h::T                                            # time step
    nd::Int                                         # implicit dynamics constraint size
    nr::Int                                         # size of a one-time-step block
    ind::CoreIndex{nq,nu,nc,nb,n1,n2,n3}            # indices of a one-time-step block
    cost::CostFunction{T,nq,nu,nc,nb}               # cost function
    n_opts::NewtonOptions{T}                        # Newton solver options
end

function Newton(H::Int, h::T, model::ContactDynamicsModel;
    cost::CostFunction=CostFunction(H,model.dim),
    n_opts::NewtonOptions=NewtonOptions()) where {T}
    dim = model.dim
    nq = dim.q # configuration
    nu = dim.u # control
    nw = dim.w # disturbance
    nc = dim.c # contact
    nb = dim.b # linear friction
    nz = num_var(dim)
    nθ = num_data(dim)
    nd = nq + nc + nb # implicit dynamics constraint
    nr = nq+nu+nc+nb+nd # size of a one-time-step block

    j = Jacobian(H,dim)
    r = Residual(H,dim)
    r̄ = Residual(H,dim)
    Δ = Residual(H,dim)
    ν = [zeros(SizedVector{nd,T}) for t=1:H]
    ν_ = deepcopy(ν)
    trial_traj = contact_trajectory(H, h, model)
    traj = contact_trajectory(H, h, model)
    Δq  = [zeros(SizedVector{nq,T}) for t=1:H]
    Δu  = [zeros(SizedVector{nu,T}) for t=1:H]
    Δγ  = [zeros(SizedVector{nc,T}) for t=1:H]
    Δb  = [zeros(SizedVector{nb,T}) for t=1:H]
    ind = CoreIndex(H, dim)
    return Newton{T,nq,nu,nw,nc,nb,nz,nθ,nd,nd,2nq+nu,
        }(
        j, r, r̄, Δ, ν, ν_, traj, trial_traj, Δq, Δu, Δγ, Δb, H, h, nd, nr, ind, cost, n_opts,
        )
end

function jacobian!(model::ContactDynamicsModel, core::Newton, jcb::Jacobian,
    impl::ImplicitTraj{T}) where {T}

    # @show "Recompute implicit dynamics"
    # implicit_dynamics!(model, core.traj, impl; κ=core.traj.κ)

    # unpack
    H = core.H
    cost = core.cost
    n_opts = core.n_opts
    fill!(jcb.j, 0.0)

    for t=1:H
        # Cost function
        jcb.Qq2[t] .+= cost.Qq[t]
        jcb.Qu1[t] .+= cost.Qu[t]
        jcb.Qγ1[t] .+= cost.Qγ[t]
        jcb.Qb1[t] .+= cost.Qb[t]
        # Implicit dynamics
        jcb.IV[t][diagind(jcb.IV[t])]   .+= - 1.0
        jcb.ITV[t][diagind(jcb.ITV[t])] .+= - 1.0
        if t >=3
            jcb.q0[t-2]  .+= impl.δq0[t]
            jcb.q0T[t-2] .+= impl.δq0[t]'
        end
        if t >= 2
            jcb.q1[t-1]  .+= impl.δq1[t]
            jcb.q1T[t-1] .+= impl.δq1[t]'
        end
        jcb.u1[t]  .+= impl.δu1[t]
        jcb.u1T[t] .+= impl.δu1[t]'
        # Dual regularization
        jcb.reg[t][diagind(jcb.reg[t])] .+= -n_opts.β * impl.lin[t].κ0 # TODO sort the κ stuff, maybe make it a prameter of this function
    end
    return nothing
end

function residual!(model::ContactDynamicsModel, core::Newton, res::Residual,
    ν::Vector{D},#Vector{SizedArray{Tuple{nd},T,1,1,Array{T,1}}},
    impl::ImplicitTraj{T}, traj::ContactTraj{T,nq,nu,nc,nb}, ref_traj::ContactTraj{T,nq,nu,nc,nb},
    ) where {T,nq,nu,nc,nb,D}

    # unpack
    n_opts = core.n_opts
    cost = core.cost
    res.r .= 0.0

    # @show "Recompute implicit dynamics"
    # implicit_dynamics!(model, traj, impl; κ=traj.κ)
    # plt = plot()
    # plot!(hcat(Vector.(traj.q)...)', label="q_residual!")
    # display(plt)
    for t in eachindex(ν)
        # Cost function
        delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2])
        delta!(core.Δu[t], traj.u[t], ref_traj.u[t])
        delta!(core.Δγ[t], traj.γ[t], ref_traj.γ[t])
        delta!(core.Δb[t], traj.b[t], ref_traj.b[t])
        # mul!(res.q2[t], cost.Qq[t], core.Δq[t])
        # mul!(res.u1[t], cost.Qu[t], core.Δu[t])
        # mul!(res.γ1[t], cost.Qγ[t], core.Δγ[t])
        # mul!(res.b1[t], cost.Qb[t], core.Δb[t])

        res.q2[t] .+= cost.Qq[t]*core.Δq[t]
        res.u1[t] .+= cost.Qu[t]*core.Δu[t]
        res.γ1[t] .+= cost.Qγ[t]*core.Δγ[t]
        res.b1[t] .+= cost.Qb[t]*core.Δb[t]

        # Implicit dynamics
        # set!(res.rd[t], impl.d[t])
        res.rd[t] .+= impl.d[t]

        # Minus Identity term #∇qk1, ∇γk, ∇bk
        # setminus!(res.rI[t], ν[t])
        res.rI[t] .+= - ν[t]
        # Implicit function theorem part #∇qk_1, ∇qk, ∇uk
        # t >= 3 ? mul!(res.q0[t-2], impl.δq0[t]', ν[t]) : nothing
        # t >= 2 ? mul!(res.q1[t-1], impl.δq1[t]', ν[t]) : nothing
        # mul!(res.u1[t], impl.δu1[t]', ν[t])
        t >= 3 ? res.q0[t-2] .+= impl.δq0[t]'*ν[t] : nothing
        t >= 2 ? res.q1[t-1] .+= impl.δq1[t]'*ν[t] : nothing
        res.u1[t] .+= impl.δu1[t]'*ν[t]
    end
    return nothing
end

function delta!(Δx::SizedArray{Tuple{nx},T,1,1}, x::SizedArray{Tuple{nx},T,1,1},
    x_ref::SizedArray{Tuple{nx},T,1,1}) where {nx,T}
    Δx .= x
    Δx .-= x_ref
    return nothing
end

function delta!(Δx::SizedArray{Tuple{nx},T,1,1}, x::AbstractArray{T},
    x_ref::AbstractArray{T}) where {nx,T}
    Δx .= x
    Δx .-= x_ref
    return nothing
end

function set!(s::SubArray, x::SizedArray{Tuple{nx},T,1,1}) where {nx,T}
    s .= x
    return nothing
end

function setminus!(s::SubArray, x::SizedArray{Tuple{nx},T,1,1}) where {nx,T}
    s .= -1.0.*x
    return nothing
end

function set_traj!(target::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
        source::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
        νtarget::Vector{SizedArray{Tuple{n1},T,1,1}},
        νsource::Vector{SizedArray{Tuple{n1},T,1,1}},
        Δ::Residual{T},
        α::T,
        ) where {T,nq,nu,nw,nc,nb,nz,nθ,n1}
    # Check that trajectory propoerties match
    H = target.H
    @assert H == source.H
    @assert (target.h - source.h)/target.h < 1e-4
    @assert (target.κ[1] - source.κ[1])/(target.κ[1]+1e-10) < 1e-4

    for t = 1:H
        target.q[t+2] .= source.q[t+2] .+ α.*Δ.q2[t]
        target.u[t] .= source.u[t] .+ α.*Δ.u1[t]
        # target.w[t] .= source.w[t] + α.*Δ.w1[t]
        target.γ[t] .= source.γ[t] .+ α.*Δ.γ1[t]
        target.b[t] .= source.b[t] .+ α.*Δ.b1[t]
        # target.z[t] .= source.z[t] + α.*Δ.z[t]
        # target.θ[t] .= source.θ[t] + α.*Δ.θ[t]
        νtarget[t]  .= νsource[t] .+ α.*Δ.rd[t]
    end
    update_z!(target)
    update_θ!(target)
    return nothing
end

function copy_traj!(target::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
        source::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}, H::Int) where {T,nq,nu,nw,nc,nb,nz,nθ,n1}
    Ht = target.H
    Hs = target.H
    @assert Hs >= H
    @assert Ht >= H
    # target.h = source.h
    target.κ .= source.κ
    for t in eachindex(1:H+2)
        target.q[t] .= source.q[t]
    end
    for t in eachindex(1:H)
        target.u[t] .= source.u[t]
        target.w[t] .= source.w[t]
        target.γ[t] .= source.γ[t]
        target.b[t] .= source.b[t]
        target.z[t] .= source.z[t]
        target.θ[t] .= source.θ[t]
    end
    return nothing
end

function reset!(core::Newton, ref_traj::ContactTraj; warm_start::Bool=false,
	initial_offset::Bool=false, q0=ref_traj.q[1], q1=ref_traj.q[2]) where {T}
    n_opts = core.n_opts
	# Reset β value
    n_opts.β = n_opts.β_init
	if !warm_start
		# Reset duals
		for t = 1:core.H
			core.ν[t] .= 0.0
			core.ν_[t] .= 0.0
		end
		# Set up trajectory
        copy_traj!(core.traj, ref_traj, core.traj.H)
		if initial_offset
		    rd = -0.03*[[1,1]; ones(model.dim.q-2)]
		    core.traj.q[1] .+= rd
		    core.traj.q[2] .+= rd
			update_θ!(core.traj, 1)
			update_θ!(core.traj, 2)
		end
	end
    for t = 1:core.H
        core.ν[t] .= 0.0
        core.ν_[t] .= 0.0
    end
    core.traj.q[1] .= deepcopy(q0)
    core.traj.q[2] .= deepcopy(q1)
    update_θ!(core.traj, 1)
    update_θ!(core.traj, 2)
	# Set up traj trial
	core.trial_traj = deepcopy(core.traj)
	return nothing
end

function newton_solve!(model::ContactDynamicsModel, core::Newton, impl::ImplicitTraj{T},
    ref_traj::ContactTraj{T,nq,nu,nc,nb} ; warm_start::Bool=false, initial_offset::Bool=false,
    q0=ref_traj.q[1], q1=ref_traj.q[2]) where {T,nq,nu,nc,nb}

	reset!(core, ref_traj; warm_start=warm_start, initial_offset=initial_offset, q0=q0, q1=q1)
    # for i = 1:core.n_opts.solver_outer_iter
        # (core.n_opts.live_plot) && (visualize!(vis, model, traj.q, Δt=h))
    for l = 1:core.n_opts.solver_inner_iter
		# Compute implicit dynamics about traj
		implicit_dynamics!(model, core.traj, impl; κ=core.traj.κ)
		# plt = plot(legend=false)
		# 	plot!([log(10, norm(d)) for d in impl.d], ylims=log.(10,(1e-11,1e1)), label="ν")
		# display(plt)

        # Compute residual
        residual!(model, core, core.r, core.ν, impl, core.traj, ref_traj)
        # Compute Jacobian
        jacobian!(model, core, core.j, impl)
        core.Δ.r .= - core.j.j \ core.r.r

		println("res:", scn(norm(core.r.r,1)/length(core.r.r), digits=3))
        if norm(core.r.r,1)/length(core.r.r) < core.n_opts.r_tol
            break
        end

        # line search the step direction
        α = 1.0
        iter = 0
        # plt = plot(legend=false)
        # # plot!([norm(d) for d in impl.d], label="ν")
        # # plot!(hcat(Vector.(core.ν)...)', label="ν")
        # # plot!(hcat(Vector.(core.ν_)...)', linewidth=3.0, linestyle=:dot, label="ν_")
        # plot!(hcat(Vector.(core.traj.q)...)', label="q")
        # plot!(hcat(Vector.(ref_traj.q)...)', linewidth=3.0, linestyle=:dot, label="q")
        # # plot!(hcat(Vector.(core.traj.u)...)', label="u")
        # # plot!(hcat(Vector.(ref_traj.u)...)', linewidth=3.0, linestyle=:dot, label="u")
        # display(plt)

        set_traj!(core.trial_traj, core.traj, core.ν_, core.ν, core.Δ, α)
		# Compute implicit dynamics about trial_traj
		implicit_dynamics!(model, core.trial_traj, impl; κ=core.trial_traj.κ)
        residual!(model, core, core.r̄, core.ν_, impl, core.trial_traj, ref_traj)
        while norm(core.r̄.r)^2.0 >= (1.0 - 0.001 * α) * norm(core.r.r)^2.0
            α = 0.5 * α
            iter += 1
            if iter > 6
                break
            end
            set_traj!(core.trial_traj, core.traj, core.ν_, core.ν, core.Δ, α)
			# Compute implicit dynamics about trial_traj
			implicit_dynamics!(model, core.trial_traj, impl; κ=core.trial_traj.κ)
            residual!(model, core, core.r̄, core.ν_, impl, core.trial_traj, ref_traj)
        end

        # update # maybe not useful never activated
        if iter > 6
            core.n_opts.β = min(core.n_opts.β*1.3, 1e2)
        else
            core.n_opts.β = max(1e1, core.n_opts.β/1.3)
        end
        println(" l: ", l ,
            "     r̄: ", scn(norm(core.r̄.r,1)/length(core.r̄.r), digits=0),
            "     r: ", scn(norm(core.r.r,1)/length(core.r.r), digits=0),
            "     Δ: ", scn(norm(core.Δ.r,1)/length(core.Δ.r), digits=0),
            "     α: ", -Int(round(log(α))),
            "     κ: ", scn(core.traj.κ[1], digits=0) ,
            )
        set_traj!(core.traj, core.traj, core.ν, core.ν, core.Δ, α)
    end
    #     κ /= 10.0
    # end
    return nothing
end




model = get_model("quadruped")
κ = 1e-4
ref_traj0 = get_trajectory("quadruped", "gait1")
ref_traj0.κ .= κ
H = ref_traj0.H
h = 0.1
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq+nc+nb
nr = nq+nu+nc+nb+nd

# Test Jacobian!
cost0 = CostFunction(H, model.dim,
    Qq=fill(Diagonal(1e-2*SizedVector{nq}([0.02,0.02,1,.15,.15,.15,.15,.15,.15,.15,.15,])), H),
    Qu=fill(Diagonal(3e-2*ones(SizedVector{nu})), H),
    Qγ=fill(Diagonal(1e-6*ones(SizedVector{nc})), H),
    Qb=fill(Diagonal(1e-6*ones(SizedVector{nb})), H),
    )
n_opts0 = NewtonOptions(r_tol=3e-4)
core0 = Newton(H, h, model, cost=cost0, n_opts=n_opts0)
impl0 = ImplicitTraj(H, model)
linearization!(model, ref_traj0, impl0)
implicit_dynamics!(model, ref_traj0, impl0)

newton_solve!(model, core0, impl0, ref_traj0)
