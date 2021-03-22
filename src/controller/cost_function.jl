mutable struct CostFunction{T,nq,nu,nc,nb}
    H::Int
    Qq::Vector{Diagonal{T,SizedArray{Tuple{nq},T,1,1,Array{T,1}}}}
    Qu::Vector{Diagonal{T,SizedArray{Tuple{nu},T,1,1,Array{T,1}}}}
    Qγ::Vector{Diagonal{T,SizedArray{Tuple{nc},T,1,1,Array{T,1}}}}
    Qb::Vector{Diagonal{T,SizedArray{Tuple{nb},T,1,1,Array{T,1}}}}
end

function CostFunction(H::Int, dim::Dimensions;
    Qq::Vector{Diagonal{T,SizedArray{Tuple{nq},T,1,1,Array{T,1}}}}=fill(Diagonal(zeros(SizedVector{dim.q})), H),
    Qu::Vector{Diagonal{T,SizedArray{Tuple{nu},T,1,1,Array{T,1}}}}=fill(Diagonal(zeros(SizedVector{dim.u})), H),
    Qγ::Vector{Diagonal{T,SizedArray{Tuple{nc},T,1,1,Array{T,1}}}}=fill(Diagonal(zeros(SizedVector{dim.c})), H),
    Qb::Vector{Diagonal{T,SizedArray{Tuple{nb},T,1,1,Array{T,1}}}}=fill(Diagonal(zeros(SizedVector{dim.b})), H),
    ) where {T,nq,nu,nc,nb}
    return CostFunction{T,nq,nu,nc,nb}(H, Qq, Qu, Qγ, Qb)
end
