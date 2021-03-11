# Newton solver options
@with_kw mutable struct NewtonOptions{T}
    r_tol::T = 1.0e-5            # primal dual residual tolerance
    solver_outer_iter::Int = 3   # outer iter on the κ parameter
    solver_inner_iter::Int = 10  # primal dual iter
    κ_init::T = 1e-3             # inner solver intialization
    κ_scale::T = 0.1             # inner solver scaling
    κ_tol::T = 2.0e-3            # inner solver tolerance
    β::T = 1e1                   # dual regularization
end

mutable struct Residual11{T,
    vq,vu,vγ,vb,
    vd,vI,vq0,vq1,vu1,
    }
    r::Vector{T}                           # residual

    qq::Vector{vq}                         # rsd cost function views
    qu::Vector{vu}                         # rsd cost function views
    qγ::Vector{vγ}                         # rsd cost function views
    qb::Vector{vb}                         # rsd cost function views

    rd::Vector{vd}                         # rsd dynamics lagrange multiplier views
    rI::Vector{vI}                         # rsd dynamics -I views
    rq0::Vector{vq0}                       # rsd dynamics q0 views
    rq1::Vector{vq1}                       # rsd dynamics q1 views
    ru1::Vector{vu1}                       # rsd dynamics u1 views
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

    qq  = [view(r, (t-1)*nr .+ iq) for t=1:H]
    qu  = [view(r, (t-1)*nr .+ iu) for t=1:H]
    qγ  = [view(r, (t-1)*nr .+ iγ) for t=1:H]
    qb  = [view(r, (t-1)*nr .+ ib) for t=1:H]

    rd  = [view(r, (t-1)*nr .+ iν) for t=1:H]
    rI  = [view(r, (t-1)*nr .+ iz) for t=1:H]
    rq0 = [view(r, (t-3)*nr .+ iq) for t=3:H]
    rq1 = [view(r, (t-2)*nr .+ iq) for t=2:H]
    ru1 = [view(r, (t-1)*nr .+ iu) for t=1:H]

    return Residual11{T,
        eltype.((qq, qu, qγ, qb))...,
        eltype.((rd, rI, rq0, rq1, ru1))...,
        }(
        r,
        qq, qu, qγ, qb,
        rd, rI, rq0, rq1, ru1,
        )
end


mutable struct Newton31{T,nq,nu,nc,nb,n1,n2,n3,
    Vq,Vu,Vγ,Vb,
    # vq,vu,vγ,vb,
    VI,VIT,Vq0,Vq0T,Vq1,Vq1T,Vu1,Vu1T,Vreg,
    # vd,vI,vq0,vq1,vu1,
    }
    jcb::Any                               # Jacobian
    # r::Vector{T}                           # residual
    r::Residual11{T}                           # residual
    # r̄::Vector{T}                           # candidate residual
    r̄::Residual11{T}                       # candidate residual
    ν::Vector{SizedArray{Tuple{n1},T,1,1}} # implicit dynamics lagrange multiplier
    Δq::Vector{SizedArray{Tuple{nq},T,1,1}} # difference between the traj and ref_traj
    Δu::Vector{SizedArray{Tuple{nu},T,1,1}} # difference between the traj and ref_traj
    Δγ::Vector{SizedArray{Tuple{nγ},T,1,1}} # difference between the traj and ref_traj
    Δb::Vector{SizedArray{Tuple{nb},T,1,1}} # difference between the traj and ref_traj
    H::Int                                 # horizon
    nd::Int                                # implicit dynamics constraint size
    nr::Int                                # size of a one-time-step block
    iq::SizedVector{nq,Int}                # configuration indices
    iu::SizedVector{nu,Int}                # control indices
    iγ::SizedVector{nc,Int}                # impact indices
    ib::SizedVector{nb,Int}                # linear friction indices
    iν::SizedVector{n1,Int}                # implicit dynamics lagrange multiplier
    iz::SizedVector{n2,Int}                # IP solver solution [q2, γ1, b1]
    iθ::SizedVector{n3,Int}                # IP solver data [q0, q1, u1]

    Qq::Vector{Vq}                         # jcb cost function views
    Qu::Vector{Vu}                         # jcb cost function views
    Qγ::Vector{Vγ}                         # jcb cost function views
    Qb::Vector{Vb}                         # jcb cost function views

    # qq::Vector{vq}                         # rsd cost function views
    # qu::Vector{vu}                         # rsd cost function views
    # qγ::Vector{vγ}                         # rsd cost function views
    # qb::Vector{vb}                         # rsd cost function views

    IV::Vector{VI}                         # jcb dynamics -I views
    ITV::Vector{VIT}                       # jcb dynamics -I views transposed
    q0::Vector{Vq0}                        # jcb dynamics q0 views
    q0T::Vector{Vq0T}                      # jcb dynamics q0 views transposed
    q1::Vector{Vq1}                        # jcb dynamics q1 views
    q1T::Vector{Vq1T}                      # jcb dynamics q1 views transposed
    u1::Vector{Vu1}                        # jcb dynamics u1 views
    u1T::Vector{Vu1T}                      # jcb dynamics u1 views transposed
    reg::Vector{Vreg}                      # jcb dual regularization views

    # rd::Vector{vd}                         # rsd dynamics lagrange multiplier views
    # rI::Vector{vI}                         # rsd dynamics -I views
    # rq0::Vector{vq0}                       # rsd dynamics q0 views
    # rq1::Vector{vq1}                       # rsd dynamics q1 views
    # ru1::Vector{vu1}                       # rsd dynamics u1 views
end

function Newton31(H::Int, dim::Dimensions)
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
    r = zeros(H*nr)
    r̄ = zeros(H*nr)

    ν = [zeros(SizedVector{nd}) for t=1:H]
    Δq  = [zeros(SizedVector{nq}) for t=1:H]
    Δu  = [zeros(SizedVector{nu}) for t=1:H]
    Δγ  = [zeros(SizedVector{nγ}) for t=1:H]
    Δb  = [zeros(SizedVector{nb}) for t=1:H]

    Qq  = [view(jcb, (t-1)*nr .+ iq, (t-1)*nr .+ iq) for t=1:H]
    Qu  = [view(jcb, (t-1)*nr .+ iu, (t-1)*nr .+ iu) for t=1:H]
    Qγ  = [view(jcb, (t-1)*nr .+ iγ, (t-1)*nr .+ iγ) for t=1:H]
    Qb  = [view(jcb, (t-1)*nr .+ ib, (t-1)*nr .+ ib) for t=1:H]

    qq  = [view(r, (t-1)*nr .+ iq) for t=1:H]
    qu  = [view(r, (t-1)*nr .+ iu) for t=1:H]
    qγ  = [view(r, (t-1)*nr .+ iγ) for t=1:H]
    qb  = [view(r, (t-1)*nr .+ ib) for t=1:H]

    IV  = [view(jcb, (t-1)*nr .+ iz, (t-1)*nr .+ iν) for t=1:H]
    ITV = [view(jcb, (t-1)*nr .+ iν, (t-1)*nr .+ iz) for t=1:H]
    q0  = [view(jcb, (t-1)*nr .+ iν, (t-3)*nr .+ iq) for t=3:H]
    q0T = [view(jcb, (t-3)*nr .+ iq, (t-1)*nr .+ iν) for t=3:H]
    q1  = [view(jcb, (t-1)*nr .+ iν, (t-2)*nr .+ iq) for t=2:H]
    q1T = [view(jcb, (t-2)*nr .+ iq, (t-1)*nr .+ iν) for t=2:H]
    u1  = [view(jcb, (t-1)*nr .+ iν, (t-1)*nr .+ iu) for t=1:H]
    u1T = [view(jcb, (t-1)*nr .+ iu, (t-1)*nr .+ iν) for t=1:H]
    reg = [view(jcb, (t-1)*nr .+ iν, (t-1)*nr .+ iν) for t=1:H]

    rd  = [view(r, (t-1)*nr .+ iν) for t=1:H]
    rI  = [view(r, (t-1)*nr .+ iz) for t=1:H]
    rq0 = [view(r, (t-3)*nr .+ iq) for t=3:H]
    rq1 = [view(r, (t-2)*nr .+ iq) for t=2:H]
    ru1 = [view(r, (t-1)*nr .+ iu) for t=1:H]

    return Newton31{T,nq,nu,nc,nb,nd,nd,2nq+nu,
        eltype.((Qq, Qu, Qγ, Qb))...,
        # eltype.((qq, qu, qγ, qb))...,
        eltype.((IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg))...,
        # eltype.((rd, rI, rq0, rq1, ru1))...,
        }(
        jcb, r, r̄, ν, Δq, Δu, Δγ, Δb, H, nd, nr,
        iq, iu, iγ, ib, iν, iz, iθ,
        Qq, Qu, Qγ, Qb,
        # qq, qu, qγ, qb,
        IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg,
        # rd, rI, rq0, rq1, ru1,
        )
end

function jacobian!(core::Newton31, impl::ImplicitTraj11{T},
    cost::CostFunction, n_opts::NewtonOptions{T}) where {T}
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
        # Dual regularization
        core.reg[t][diagind(core.reg[t])] .= -n_opts.β * impl.lin[t].κ0 # TODO sort the κ stuff, maybe make it a prameter of this function
    end
    return nothing
end

function residual!(core::Newton31, impl::ImplicitTraj11{T},
    cost::CostFunction, traj::ContactTraj{T,nq,nu,nc,nb}, ref_traj::ContactTraj{T,nq,nu,nc,nb},
    n_opts::NewtonOptions{T}) where {T,nq,nu,nc,nb}
    # unpack
    r = core.r
    r .= 0.0

    for t in eachindex(core.ν)
        # Cost function
        delta!(core.Δq[t], traj.q[t], ref_traj.q[t])
        delta!(core.Δu[t], traj.u[t], ref_traj.u[t])
        delta!(core.Δγ[t], traj.γ[t], ref_traj.γ[t])
        delta!(core.Δb[t], traj.b[t], ref_traj.b[t])
        mul!(core.qq[t], cost.Qq[t], core.Δq[t])
        mul!(core.qu[t], cost.Qu[t], core.Δu[t])
        mul!(core.qγ[t], cost.Qγ[t], core.Δγ[t])
        mul!(core.qb[t], cost.Qb[t], core.Δb[t])
        # Implicit dynamics
        set!(core.rd[t], impl.d[t])
        # Minus Identity term #∇qk1, ∇γk, ∇bk
        setminus!(core.rI[t], core.ν[t])
        # Implicit function theorem part #∇qk_1, ∇qk, ∇uk
        t >= 3 ? mul!(core.rq0[t-2], impl.δq0[t]', core.ν[t]) : nothing
        t >= 2 ? mul!(core.rq1[t-1], impl.δq1[t]', core.ν[t]) : nothing
        mul!(core.ru1[t], impl.δu1[t]', core.ν[t])
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







T = Float64
H = 50
h = 0.03
κ = 1e-3
model = get_model("quadruped")
nq = model.dim.q
nu = model.dim.u
nγ = model.dim.c
nb = model.dim.b
core0 = Newton31(H, model.dim)
impl0 = ImplicitTraj11(H, model)

q0 = SVector{nq,T}([2.0, 2.0, zeros(nq-2)...])
q1 = SVector{nq,T}([2.0, 2.0, zeros(nq-2)...])
ip_opts = InteriorPointOptions(κ_init=κ, κ_tol=κ*2, r_tol=1e-5)
sim0 = simulator2_base(model, q0, q1, h, H;
    u = [@SVector zeros(model.dim.u) for t = 1:H],
    w = [@SVector zeros(model.dim.w) for t = 1:H],
    ip_opts = ip_opts,
    sim_opts = SimulatorOptions{T}())
simulate!(sim0; verbose = false)
ref_traj0 = deepcopy(sim0.traj)
traj0 = deepcopy(sim0.traj)
linearization!(model, ref_traj0, impl0)
implicit_dynamics!(model, ref_traj0, impl0, κ=κ)

cost0 = CostFunction(H, model.dim,
    Qq=fill(Diagonal(1e-1*ones(SizedVector{nq})), H),
    Qu=fill(Diagonal(1e-2*ones(SizedVector{nu})), H),
    Qγ=fill(Diagonal(1e-3*ones(SizedVector{nγ})), H),
    Qb=fill(Diagonal(1e-4*ones(SizedVector{nb})), H),
    )
n_opts = NewtonOptions()
@allocated residual!(core0, impl0, cost0, traj0, ref_traj0, n_opts)
@allocated residual!(core0, impl0, cost0, traj0, ref_traj0, n_opts)
@allocated residual!(core0, impl0, cost0, traj0, ref_traj0, n_opts)
@code_warntype residual!(core0, impl0, cost0, traj0, ref_traj0, n_opts)
@benchmark residual!(core0, impl0, cost0, traj0, ref_traj0, n_opts)

@allocated jacobian!(core0, impl0, cost0, n_opts)
@allocated jacobian!(core0, impl0, cost0, n_opts)
@allocated jacobian!(core0, impl0, cost0, n_opts)
@allocated jacobian!(core0, impl0, cost0, n_opts)
@benchmark jacobian!(core0, impl0, cost0, n_opts)





plot(Gray.(Matrix((1e3.*core0.jcb))))
