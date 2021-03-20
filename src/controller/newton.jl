# Newton23 solver options
@with_kw mutable struct Newton23Options{T}
    r_tol::T = 1.0e-5            # primal dual residual tolerance
    solver_outer_iter::Int = 3   # outer iter on the κ parameter
    solver_inner_iter::Int = 10  # primal dual iter
    κ_init::T = 1e-3             # inner solver intialization
    κ_scale::T = 0.1             # inner solver scaling
    κ_tol::T = 2.0e-3            # inner solver tolerance
    β_init::T = 1e1              # dual regularization
    β::T = 1e1                   # dual regularization
    live_plots::Bool = false    # Visulaize the trajectory during the solve
end

mutable struct Residual11{T,
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

function Residual11(H::Int, dim::Dimensions)
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

    return Residual11{T,
        eltype.((q2, u1, γ1, b1))...,
        eltype.((rd, rI, q0, q1))...,
        }(
        r,
        q2, u1, γ1, b1,
        rd, rI, q0, q1,
        )
end

mutable struct Jacobian12{T,
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

function Jacobian12(H::Int, dim::Dimensions)
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

    return Jacobian12{eltype(j),
        eltype.((Qq2, Qu1, Qγ1, Qb1))...,
        eltype.((IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg))...,
        }(
        j,
        Qq2, Qu1, Qγ1, Qb1,
        IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg,
        )
end

mutable struct CoreIndex14{nq,nu,nc,nb,n1,n2,n3,I<:Int}
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

function CoreIndex14(H::Int, dim::Dimensions)
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

    return CoreIndex14{nq,nu,nc,nb,nd,nd,2nq+nu,Int}(
        H, nd, nr,
        iq, iu, iγ, ib, iν, iz, iθ,
        Iq, Iu, Iγ, Ib, Iν, Iz, Iθ,
        )
end

mutable struct Newton23{T,nq,nu,nw,nc,nb,nz,nθ,n1,n2,n3,
    }
    j::Jacobian12{T}                                  # Jacobian12
    r::Residual11{T}                                  # residual
    r̄::Residual11{T}                                  # candidate residual
    Δ::Residual11{T}                                  # step direction in the Newton solve, it contains: q2-qH+1, u1-uH, γ1-γH, b1-bH, λd1-λdH
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
    ind::CoreIndex14{nq,nu,nc,nb,n1,n2,n3}            # indices of a one-time-step block
end

function Newton23(H::Int, h::T, model::ContactDynamicsModel) where {T}
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

    j = Jacobian12(H,dim)
    r = Residual11(H,dim)
    r̄ = Residual11(H,dim)
    Δ = Residual11(H,dim)
    ν = [zeros(SizedVector{nd,T}) for t=1:H]
    ν_ = deepcopy(ν)
    trial_traj = contact_trajectory(H, h, model)
    traj = contact_trajectory(H, h, model)
    Δq  = [zeros(SizedVector{nq,T}) for t=1:H]
    Δu  = [zeros(SizedVector{nu,T}) for t=1:H]
    Δγ  = [zeros(SizedVector{nc,T}) for t=1:H]
    Δb  = [zeros(SizedVector{nb,T}) for t=1:H]
    ind = CoreIndex14(H, dim)
    return Newton23{T,nq,nu,nw,nc,nb,nz,nθ,nd,nd,2nq+nu,
        }(
        j, r, r̄, Δ, ν, ν_, traj, trial_traj, Δq, Δu, Δγ, Δb, H, h, nd, nr, ind,
        )
end

function jacobian!(model::ContactDynamicsModel, core::Newton23, jcb::Jacobian12,
    impl::ImplicitTraj{T}, cost::CostFunction, n_opts::Newton23Options{T}) where {T}

    # @show "Recompute implicit dynamics"
    # implicit_dynamics!(model, core.traj, impl; κ=core.traj.κ)

    # unpack
    H = core.H
    # jcb = core.j
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

function residual!(model::ContactDynamicsModel, core::Newton23, res::Residual11,
    ν::Vector{D},#Vector{SizedArray{Tuple{nd},T,1,1,Array{T,1}}},
    impl::ImplicitTraj{T}, cost::CostFunction,
    traj::ContactTraj{T,nq,nu,nc,nb}, ref_traj::ContactTraj{T,nq,nu,nc,nb},
    n_opts::Newton23Options{T}) where {T,nq,nu,nc,nb,D}#nd}
    # unpack
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
        mul!(res.q2[t], cost.Qq[t], core.Δq[t])
        mul!(res.u1[t], cost.Qu[t], core.Δu[t])
        mul!(res.γ1[t], cost.Qγ[t], core.Δγ[t])
        mul!(res.b1[t], cost.Qb[t], core.Δb[t])

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
        Δ::Residual11{T},
        α::T,
        ) where {T,nq,nu,nw,nc,nb,nz,nθ,n1,I}
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
    return nothing
end
