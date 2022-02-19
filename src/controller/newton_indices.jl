# Configurations and forces
mutable struct NewtonIndicesConfigurationForce{nq,nu,nc,nb,n1,n2,n3} <: NewtonIndices
    nd::Int                                   # implicit dynamics constraint
    nr::Int                                   # size of a one-time-step block
    iq::SizedArray{Tuple{nq},Int,1,1}         # configuration indices
    iu::SizedArray{Tuple{nu},Int,1,1}         # control indices
    iγ::SizedArray{Tuple{nc},Int,1,1}         # impact indices
    ib::SizedArray{Tuple{nb},Int,1,1}         # linear friction indices
    iν::SizedArray{Tuple{n1},Int,1,1}         # implicit dynamics lagrange multiplier
    iz::SizedArray{Tuple{n2},Int,1,1}         # IP solver solution [q2, γ1, b1]
    iθ::SizedArray{Tuple{n3},Int,1,1}         # IP solver data [q0, q1, u1]
    Iq::Vector{SizedArray{Tuple{nq},Int,1,1}} # configuration indices
    Iu::Vector{SizedArray{Tuple{nu},Int,1,1}} # control indices
    Iγ::Vector{SizedArray{Tuple{nc},Int,1,1}} # impact indices
    Ib::Vector{SizedArray{Tuple{nb},Int,1,1}} # linear friction indices
    Iν::Vector{SizedArray{Tuple{n1},Int,1,1}} # implicit dynamics lagrange multiplier
    Iz::Vector{SizedArray{Tuple{n2},Int,1,1}} # IP solver solution [q2, γ1, b1]
    Iθ::Vector{SizedArray{Tuple{n3},Int,1,1}} # IP solver data [q0, q1, u1]
end

function NewtonIndicesConfigurationForce(model::ContactModel, env::Environment, H::Int)
    dim = model.dim
    nq = dim.q # configuration
    nu = dim.u # control
    nw = dim.w # disturbance
    nc = dim.c # contact
    nb = nc * friction_dim(env) # linear friction
    nd = nq + nc + nb # implicit dynamics constraint
    nr = nq + nu + nc + nb + nd # size of a one-time-step block

    off = 0
    iq = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q2
    iu = SizedVector{nu}(off .+ (1:nu)); off += nu # index of the control u1
    iγ = SizedVector{nc}(off .+ (1:nc)); off += nc # index of the impact γ1
    ib = SizedVector{nb}(off .+ (1:nb)); off += nb # index of the linear friction b1
    iν = SizedVector{nd}(off .+ (1:nd)); off += nd # index of the dynamics lagrange multiplier ν1
    iz = vcat(iq, iγ, ib) # index of the IP solver solution [q2, γ1, b1]
    iθ = vcat(iq .- 2nr, iq .- nr, iu) # index of the IP solver data [q0, q1, u1]

    Iq = [(t - 1) * nr .+ iq for t = 1:H]
    Iu = [(t - 1) * nr .+ iu for t = 1:H]
    Iγ = [(t - 1) * nr .+ iγ for t = 1:H]
    Ib = [(t - 1) * nr .+ ib for t = 1:H]
    Iν = [(t - 1) * nr .+ iν for t = 1:H]
    Iz = [(t - 1) * nr .+ iz for t = 1:H]
    Iθ = [(t - 1) * nr .+ iθ for t = 1:H]

    return NewtonIndicesConfigurationForce{nq,nu,nc,nb,nd,nd,2nq+nu}(
        nd, nr,
        iq, iu, iγ, ib, iν, iz, iθ,
        Iq, Iu, Iγ, Ib, Iν, Iz, Iθ)
end

mutable struct NewtonIndicesConfiguration{nq,nu,n1,n2,n3} <: NewtonIndices
    nd::Int                                   # implicit dynamics constraint
    nr::Int                                   # size of a one-time-step block
    iq::SizedArray{Tuple{nq},Int,1,1}         # configuration indices
    iu::SizedArray{Tuple{nu},Int,1,1}         # control indices
    iν::SizedArray{Tuple{n1},Int,1,1}         # implicit dynamics lagrange multiplier
    iz::SizedArray{Tuple{n2},Int,1,1}         # IP solver solution [q2]
    iθ::SizedArray{Tuple{n3},Int,1,1}         # IP solver data [q0, q1, u1]
    Iq::Vector{SizedArray{Tuple{nq},Int,1,1}} # configuration indices
    Iu::Vector{SizedArray{Tuple{nu},Int,1,1}} # control indices
    Iν::Vector{SizedArray{Tuple{n1},Int,1,1}} # implicit dynamics lagrange multiplier
    Iz::Vector{SizedArray{Tuple{n2},Int,1,1}} # IP solver solution [q2]
    Iθ::Vector{SizedArray{Tuple{n3},Int,1,1}} # IP solver data [q0, q1, u1]
end

function NewtonIndicesConfiguration(model::ContactModel, env::Environment, H::Int)
    dim = model.dim
    nq = dim.q # configuration
    nu = dim.u # control
    nw = dim.w # disturbance
    nc = dim.c # contact
    nd = nq # implicit dynamics constraint
    nr = nq + nu + nd # size of a one-time-step block

    off = 0
    iq = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q2
    iu = SizedVector{nu}(off .+ (1:nu)); off += nu # index of the control u1
    iν = SizedVector{nd}(off .+ (1:nd)); off += nd # index of the dynamics lagrange multiplier ν1
    iz = copy(iq)                                  # index of the IP solver solution [q2]
    iθ = vcat(iq .- 2nr, iq .- nr, iu) # index of the IP solver data [q0, q1, u1]

    Iq = [(t - 1) * nr .+ iq for t = 1:H]
    Iu = [(t - 1) * nr .+ iu for t = 1:H]
    Iν = [(t - 1) * nr .+ iν for t = 1:H]
    Iz = [(t - 1) * nr .+ iz for t = 1:H]
    Iθ = [(t - 1) * nr .+ iθ for t = 1:H]

    return NewtonIndicesConfiguration{nq,nu,nd,nd,2nq+nu}(
        nd, nr,
        iq, iu, iν, iz, iθ,
        Iq, Iu, Iν, Iz, Iθ)
end

function NewtonIndices(model::ContactModel, env::Environment, H::Int;
    mode = :configurationforce)

    if mode == :configurationforce
        return NewtonIndicesConfigurationForce(model, env, H)
    elseif mode == :configuration
        NewtonIndicesConfiguration(model, env, H)
    else
        @error "mode not implemented"
    end
end
