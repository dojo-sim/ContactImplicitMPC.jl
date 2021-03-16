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

    q2::Vector{vq2}                         # rsd cost function views
    u1::Vector{vu1}                         # rsd cost function views
    γ1::Vector{vγ1}                         # rsd cost function views
    b1::Vector{vb1}                         # rsd cost function views

    rd::Vector{vd}                         # rsd dynamics lagrange multiplier views
    rI::Vector{vI}                         # rsd dynamics -I views [q2, γ1, b1]
    q0::Vector{vq0}                       # rsd dynamics q0 views
    q1::Vector{vq1}                       # rsd dynamics q1 views
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

function Newton23(H::Int, h::T, dim::Dimensions) where {T}
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
    trial_traj = ContactTraj(H, h, dim)
    traj = ContactTraj(H, h, dim)
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

function jacobian!(core::Newton23, jcb::Jacobian12, impl::ImplicitTraj{T},
    cost::CostFunction, n_opts::Newton23Options{T}) where {T}
    # unpack
    H = core.H
    # jcb = core.j
    fill!(jcb.j, 0.0)

    for t=1:H
        # Cost function
        jcb.Qq2[t] .= cost.Qq[t]
        jcb.Qu1[t] .= cost.Qu[t]
        jcb.Qγ1[t] .= cost.Qγ[t]
        jcb.Qb1[t] .= cost.Qb[t]
        # Implicit dynamics
        jcb.IV[t][diagind(jcb.IV[t])]   .= - 1.0
        jcb.ITV[t][diagind(jcb.ITV[t])] .= - 1.0
        if t >=3
            jcb.q0[t-2]  .= impl.δq0[t]
            jcb.q0T[t-2] .= impl.δq0[t]'
        end
        if t >= 2
            jcb.q1[t-1]  .= impl.δq1[t]
            jcb.q1T[t-1] .= impl.δq1[t]'
        end
        jcb.u1[t]  .= impl.δu1[t]
        jcb.u1T[t] .= impl.δu1[t]'
        # Dual regularization
        jcb.reg[t][diagind(jcb.reg[t])] .= -n_opts.β * impl.lin[t].κ0 # TODO sort the κ stuff, maybe make it a prameter of this function
    end
    return nothing
end

function residual!(model::ContactDynamicsModel, core::Newton23, res::Residual11,
    impl::ImplicitTraj{T}, cost::CostFunction,
    traj::ContactTraj{T,nq,nu,nc,nb}, ref_traj::ContactTraj{T,nq,nu,nc,nb},
    n_opts::Newton23Options{T}) where {T,nq,nu,nc,nb}
    # unpack
    res.r .= 0.0

    # @show "Recompute implicit dynamics"
    implicit_dynamics!(model, traj, impl; κ=traj.κ)

    for t in eachindex(core.ν)
        # Cost function
        delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2])
        delta!(core.Δu[t], traj.u[t], ref_traj.u[t])
        delta!(core.Δγ[t], traj.γ[t], ref_traj.γ[t])
        delta!(core.Δb[t], traj.b[t], ref_traj.b[t])
        mul!(res.q2[t], cost.Qq[t], core.Δq[t])
        mul!(res.u1[t], cost.Qu[t], core.Δu[t])
        mul!(res.γ1[t], cost.Qγ[t], core.Δγ[t])
        mul!(res.b1[t], cost.Qb[t], core.Δb[t])
        # Implicit dynamics
        set!(res.rd[t], impl.d[t])
        # Minus Identity term #∇qk1, ∇γk, ∇bk
        setminus!(res.rI[t], core.ν[t])
        # Implicit function theorem part #∇qk_1, ∇qk, ∇uk
        t >= 3 ? mul!(res.q0[t-2], impl.δq0[t]', core.ν[t]) : nothing
        t >= 2 ? mul!(res.q1[t-1], impl.δq1[t]', core.ν[t]) : nothing
        mul!(res.u1[t], impl.δu1[t]', core.ν[t])
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

function sparse_zero!(spm::SparseMatrixCSC)
	n = length(spm.nzval)
	for i = 1:n
		spm.nzval[i] = 0.0
	end
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
    @assert (target.κ - source.κ)/(target.κ+1e-10) < 1e-4

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

function newton_solve!(model::ContactDynamicsModel, core::Newton23, impl::ImplicitTraj{T},
    cost::CostFunction, ref_traj::ContactTraj{T,nq,nu,nc,nb},
    n_opts::Newton23Options{T}) where {T,nq,nu,nc,nb}

    H = core.H
    n_opts.β = n_opts.β_init
    # Reset
    # set Δ, ν to 0
    core.Δ.r .= 0.0
    for t = 1:H
        core.ν[t] .= 0.0
        core.ν_[t] .= 0.0
    end
    # Initialization
    core.traj = deepcopy(ref_traj)
    rd = 0.1*rand(model.dim.q)
    core.traj.q[1] .+= rd
    core.traj.q[2] .+= rd
    core.trial_traj = deepcopy(core.traj)

    # for i = 1:n_opts.solver_outer_iter
        # (n_opts.live_plot) && (visualize!(vis, model, traj.q, Δt=h))
    for l = 1:n_opts.solver_inner_iter
        @show l
        # Compute Jacobian12
        jacobian!(core, core.j, impl, cost, n_opts)
        # Compute residual
        residual!(model, core, core.r, impl, cost, core.traj, ref_traj, n_opts)

        core.Δ.r .= - core.j.j \ core.r.r

        if norm(core.r.r,1)/length(core.r.r) < n_opts.r_tol
            break
        end

        # line search the step direction
        α = 1.0
        iter = 0
        plt = plot()
        plot!(hcat(Vector.(core.ν)...)', label="ν")
        plot!(hcat(Vector.(core.ν_)...)', linewidth=3.0, linestyle=:dot, label="ν_")
        # plot!(hcat(Vector.(core.traj.q)...)', label="q")
        # plot!(hcat(Vector.(ref_traj.q)...)', linewidth=3.0, linestyle=:dot, label="q")
        # plot!(hcat(Vector.(core.traj.u)...)', label="u")
        # plot!(hcat(Vector.(ref_traj.u)...)', linewidth=3.0, linestyle=:dot, label="u")
        display(plt)

        set_traj!(core.trial_traj, core.traj, core.ν_, core.ν, core.Δ, α)

        residual!(model, core, core.r̄, impl, cost, core.trial_traj, ref_traj, n_opts)
        while norm(core.r̄.r)^2.0 >= (1.0 - 0.001 * α) * norm(core.r.r)^2.0
            println("     r̄: ", scn(norm(core.r̄.r,1)/length(core.r̄.r), digits=4))
            α = 0.5 * α
            # println("   α = $α")
            iter += 1
            if iter > 6
                break
            end
            set_traj!(core.trial_traj, core.traj, core.ν_, core.ν, core.Δ, α)
            residual!(model, core, core.r̄, impl, cost, core.trial_traj, ref_traj, n_opts)
        end

        # update # maybe not useful never activated
        if iter > 6
            n_opts.β = min(n_opts.β*1.3, 1e2)
        else
            n_opts.β = max(1e1, n_opts.β/1.3)
        end
        println(" κ: ", scn(core.traj.κ, digits=0) ,
            "     r̄: ", scn(norm(core.r̄.r,1)/length(core.r̄.r), digits=0),
            "     r: ", scn(norm(core.r.r,1)/length(core.r.r), digits=0),
            "     Δ: ", scn(norm(core.Δ.r,1)/length(core.Δ.r), digits=0),
            "     α: ", -Int(round(log(α))))
        set_traj!(core.traj, core.traj, core.ν, core.ν, core.Δ, α)
    end
    #     κ /= 10.0
    # end
    return nothing
end

# vis = Visualizer()
# open(vis)


n_opts.r_tol = 1e-6
core1 = Newton23(H, h, model.dim)
# @code_warntype newton_solve!(model, core1, impl0, cost0, ref_traj0, n_opts)
@time newton_solve!(model, core1, impl0, cost0, ref_traj0, n_opts)


visualize!(vis, model, ref_traj0.q, Δt=10*h)
visualize!(vis, model, core0.traj.q, Δt=10*h)
visualize!(vis, model, core0.trial_traj.q, Δt=10*h)

plot(hcat(Vector.(ref_traj0.q)...)')
plot(hcat(Vector.(core0.traj.q)...)')
plot(hcat(Vector.(core0.trial_traj.q)...)')


residual!(model, core0, core0.r, impl0, cost0, core0.traj, ref_traj0, n_opts)
norm(impl0.d[1])
norm(core0.r.r)
residual!(model, core0, core0.r, impl0, cost0, core0.trial_traj, ref_traj0, n_opts)
norm(impl0.d[1])
norm(core0.r.r)
residual!(model, core0, core0.r, impl0, cost0, ref_traj0, ref_traj0, n_opts)
norm(impl0.d[1])
norm(core0.r.r)


# @profiler fff(20)
function fff(M)
    for k = 1:M
        newton_solve!(model, core0, impl0, cost0, ref_traj0, n_opts)
    end
    return nothing
end

core0.Δ
core0.Δ.r .+= 100
core0.Δ.r
core0.Δ.q2

res = Residual11(H, model.dim)
res.r .= 100
res.q2

core1 = Newton23(H,h,model.dim)
core1.Δ.r .+= 233
core1.Δ.q2

a = 10
a = 10
a = 10
a = 10

# # Test set_traj!
# T = Float64
# H = 10
# h = 0.1
# # model = get_model("quadruped")
# target = ContactTraj(H, h, model.dim)
# source = ContactTraj(H, h, model.dim)
# nd = model.dim.q + model.dim.c + model.dim.b
# νtarget = [-30*ones(SizedVector{nd,T}) for t=1:H]
# νsource = [+30*ones(SizedVector{nd,T}) for t=1:H]
# for t = 1:H
#     source.q[t] .= +1.0
#     source.u[t] .= +2.0
#     source.w[t] .= +3.0
#     source.γ[t] .= +4.0
#     source.b[t] .= +5.0
#     source.z[t] .= +6.0
#     source.θ[t] .= +7.0
#
#     target.q[t] .= -1.0
#     target.u[t] .= -2.0
#     target.w[t] .= -3.0
#     target.γ[t] .= -4.0
#     target.b[t] .= -5.0
#     target.z[t] .= -6.0
#     target.θ[t] .= -7.0
# end
# Δ0 = Residual11(H, model.dim)
# Δ0.r .+= 100*ones(ind0.nr*H)
# set_traj!(target, source, νtarget, νsource, Δ0, 2.0)
# @code_warntype set_traj!(target, source, νtarget, νsource, Δ0, 2.0)
# @benchmark set_traj!(target, source, νtarget, νsource, Δ0, 2.0)
# source.q[1] .+= 1000.0
# source.u[1] .+= 1000.0
# @test target.q[1][1] == 201.0
# @test target.u[1][1] == 202.0
# @test target.γ[1][1] == 204.0
# @test target.b[1][1] == 205.0
# @test νtarget[1][1]  == 230.0

# ρ0=1e-3,
# β=1e1,
# outer_iter=1,
# res_tol=1e-8,
# solver_outer_iter::Int=3,
# solver_inner_iter::Int=10,
# utr_ref=[zeros(probsize.nu)  for k = 1:probsize.N-1],
# qtr_ref=[zeros(probsize.nq)  for k = 1:probsize.N-1],
# γtr_ref=[zeros(probsize.nγ)  for k = 1:probsize.N-1],
# btr_ref=[zeros(2probsize.nb) for k = 1:probsize.N-1],
# utr_wrm=utr_ref,
# qtr_wrm=qtr_ref,
# γtr_wrm=γtr_ref,
# btr_wrm=btr_ref,
# Qu=fill(Diagonal(1e-1*ones(probsize.nu)), probsize.N-1),
# Qq=fill(Diagonal(1e+1*ones(probsize.nq)), probsize.N-1),
# Qγ=fill(Diagonal(1e-7*ones(probsize.nγ)), probsize.N-1),
# Qb=fill(Diagonal(1e-7*ones(2probsize.nb)), probsize.N-1),
# u_amp=1e-1,
# live_plot::Bool=false,
# z_init=z_init,
# )
# N = probsize.N
# nq = probsize.nq
# nu = probsize.nu
# nγ = probsize.nγ
# nb = probsize.nb
# nz = nq+4nγ+4nb # residual of the step and variables of the step
# nθ = 2nq+nu # parameters of the step
# nλ = nq+nγ+2nb # dynamics lagrange multiplier
# m = nz+nu+nλ # primal dual vars in a step
# M = (N-1)*m # primal dual vars
# nr = 2nq+nu+2nγ+4nb # outer problem resiudal size
#
# traj = zeros(M)
# # Initilaization
# set_var_traj!(N,model, traj, utr_wrm, :u)
# set_var_traj!(N,model, traj, qtr_wrm, :q)
# set_var_traj!(N,model, traj, γtr_wrm, :γ)
# set_var_traj!(N,model, traj, btr_wrm, :b)
#
#
# function solver(traj)
#     # indu = get_var_ind(N,model,traj,:u)
#     # indq = get_var_ind(N,model,traj,:q)
#     # indγ = get_var_ind(N,model,traj,:γ)
#     # indb = get_var_ind(N,model,traj,:b)
#     # indλd = get_var_ind(N,model,traj,:λd)
#     # masktr = [[indu[k]; indq[k]; indγ[k]; indb[k]; indλd[k]] for k=1:N-1]
#     # mask = vcat(masktr...)
#     # unmask = setdiff(1:length(traj), mask)
#
#     for i = 1:n_opts.solver_outer_iter
#         live_plot ? visualize!(vis, model, traj.q, Δt=model.dt) : nothing
#         # sleep(1.0)
#         for l = 1:solver_inner_iter
#             j = jac(traj)
#             r = res(traj)
#             Δ = j \ r
#             if norm(r,1)/length(r) < res_tol
#                 break
#             end
#
#
#             # line search the step direction
#             α = 1.0
#             iter = 0
#             traj_trial = deepcopy(traj)
#             traj_trial[mask] -= α * Δ
#
#             # plt = plot()
#             # plot!(reshape(vcat(get_var_traj(N,model,traj,:u)...), (nu,N-1))', label="u")
#             # display(plt)
#
#             while norm(res(traj_trial))^2.0 >= (1.0 - 0.001 * α) * norm(r)^2.0
#                 α = 0.5 * α
#                 # println("   α = $α")
#                 iter += 1
#                 if iter > 6
#                     break
#                 end
#
#                 traj_trial = deepcopy(traj)
#                 traj_trial[mask] -= α * Δ
#             end
#
#             # update
#             if iter > 6
#                 β = min(β*1.3, 1e2)
#             else
#                 β = max(1e1, β/1.3)
#             end
#             println("ρ0: ", scn(ρ0, digits=0) ,
#                 "     r: ", scn(norm(r,1)/length(r), digits=0),
#                 "     Δ: ", scn(norm(Δ,1)/length(Δ), digits=0),
#                 "     α: ", -Int(round(log(α))))
#             traj[mask] -= α * Δ
#         end
#         ρ0 /= 10.0
#     end
#     return traj, flag
# end



T = Float64
κ = 1e-4
# model = get_model("quadruped")
@load joinpath(pwd(), "src/dynamics/quadruped/gaits/gait1.jld2") z̄ x̄ ū h̄ q u γ b
# time
h = h̄
H = length(u)

nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
core0 = Newton23(H, h, model.dim)
impl0 = ImplicitTraj(H, model)

# initial conditions
q0 = SVector{model.dim.q}(q[1])
q1 = SVector{model.dim.q}(q[2])

function z_initialize!(z, model::Quadruped, q1)
	nq = model.dim.q
    z .= 1.0e-1
    z[1:nq] = q1
end

sim0 = simulator2(model, q0, q1, h, H,
    u = [SVector{model.dim.u}(h * ut) for ut in u],
    r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
    rz = model.spa.rz_sp,
    rθ = model.spa.rθ_sp,
    ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 2κ, κ_init = κ),
    sim_opts = SimulatorOptions(warmstart = true))

simulate!(sim0; verbose = false)
ref_traj0 = deepcopy(sim0.traj)
traj0 = deepcopy(sim0.traj)
linearization!(model, ref_traj0, impl0)
implicit_dynamics!(model, ref_traj0, impl0, κ=κ)

cost0 = CostFunction(H, model.dim,
    Qq=fill(Diagonal(1e-2*SizedVector{nq}([0.02,0.02,1,.15,.15,.15,.15,.15,.15,.15,.15,])), H),
    Qu=fill(Diagonal(3e-2*ones(SizedVector{nu})), H),
    Qγ=fill(Diagonal(1e-100*ones(SizedVector{nc})), H),
    Qb=fill(Diagonal(1e-100*ones(SizedVector{nb})), H),
    )
n_opts = Newton23Options()
@allocated residual!(core0, core0.r, impl0, cost0, traj0, ref_traj0, n_opts)
@allocated residual!(core0, core0.r, impl0, cost0, traj0, ref_traj0, n_opts)
# @allocated residual!(core0, core0.r, impl0, cost0, traj0, ref_traj0, n_opts)
# @code_warntype residual!(core0, core0.r, impl0, cost0, traj0, ref_traj0, n_opts)
# @benchmark residual!(core0, core0.r, impl0, cost0, traj0, ref_traj0, n_opts)

@allocated jacobian!(core0, core0.j, impl0, cost0, n_opts)
@allocated jacobian!(core0, core0.j, impl0, cost0, n_opts)
# @allocated jacobian!(core0, core0.j, impl0, cost0, n_opts)
# @code_warntype jacobian!(core0, core0.j, impl0, cost0, n_opts)
# @benchmark jacobian!(core0, core0.j, impl0, cost0, n_opts)

newton_solve!(model, core0, impl0, cost0, ref_traj0, n_opts)
visualize!(vis, model, ref_traj0.q, Δt=5*h)









T = Float64
κ = 1e-4
# model = get_model("quadruped")
@load joinpath(pwd(), "src/dynamics/quadruped/gaits/gait1.jld2") z̄ x̄ ū h̄ q u γ b
# time
h = h̄
H = length(u)

# initial conditions
q0 = SVector{model.dim.q}(q[1])
q1 = SVector{model.dim.q}(q[2])

function z_initialize!(z, model::Quadruped, q1)
	nq = model.dim.q
    z .= 1.0e-1
    z[1:nq] = q1
end

sim0 = simulator2(model, q0, q1, h, H,
    u = [SVector{model.dim.u}(h * ut) for ut in u],
    r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
    rz = model.spa.rz_sp,
    rθ = model.spa.rθ_sp,
    ip_opts = InteriorPointOptions(r_tol = 1.0e-8, κ_tol = 2κ, κ_init = κ),
    sim_opts = SimulatorOptions(warmstart = true))
simulate!(sim0; verbose = false)

ref_traj0 = deepcopy(sim0.traj)
ref_traj1 = deepcopy(ref_traj0)
for t = 1:H+2
    ref_traj1.q[t] .+= 1e-0*rand(model.dim.q)
end

impl0 = ImplicitTraj(H, model)
impl1 = ImplicitTraj(H, model)
linearization!(model, ref_traj0, impl0, ref_traj0.κ)
linearization!(model, ref_traj0, impl1, ref_traj0.κ)

implicit_dynamics!(model, ref_traj0, impl0)
mean(norm.(impl0.d, 2)) < 5e-3
plot(norm.(impl0.d, 2))
implicit_dynamics!(model, ref_traj1, impl1)
mean(norm.(impl1.d, 2)) > 5e-1
plot(norm.(impl1.d, 2))



plot(Gray.(Matrix((1e3.*core0.j.j[1:300, 1:300]))))


# Solves Ax = b using LDL factors for A.
# Solves in place (x replaces b)
function solve!(F::QDLDLFactorisation,b)

    #bomb if logical factorisation only
    if F.logical
        error("Can't solve with logical factorisation only")
    end

    #permute b
    tmp = F.perm == nothing ? b : permute!(F.workspace.fwork,b,F.perm)

    QDLDL_solve!(F.workspace.Ln,
                 F.workspace.Lp,
                 F.workspace.Li,
                 F.workspace.Lx,
                 F.workspace.Dinv,
                 tmp)

    #inverse permutation
    b = F.perm == nothing ? tmp : ipermute!(b,F.workspace.fwork,F.perm)

    return nothing
end




A = deepcopy(core0.jcb)
b = deepcopy(core0.r.r) + 0.1*rand(length(core0.r.r))
A1 = deepcopy(A)
A1[1:10, 1:10] += Diagonal(rand(10))
A2 = deepcopy(A1)
A2[1:1000, 1:1000] += Diagonal(rand(1000))

rank(A)
rank(A1)
rank(A2)


F = qdldl(A)
x1 = deepcopy(b)
QDLDL.solve!(F,x1)
norm(A*x1 - b, Inf)

G = QDLDLFactorisationAF(A,F)
qdldl!(A,G)
x2 = deepcopy(b)
QDLDL.solve!(G.F,x2)
norm(b, Inf)
norm(A*x2 - b, Inf)


qdldl!(A1,G)
x4 = deepcopy(b)
QDLDL.solve!(G.F,x4)
norm(A1*x4 - b, Inf)

@benchmark qdldl!(A2,G)
x5 = deepcopy(b)
@benchmark QDLDL.solve!(G.F,x5)
norm(A*x5 - b, Inf)
norm(A1*x5 - b, Inf)
norm(A2*x5 - b, Inf)
