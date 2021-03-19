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
