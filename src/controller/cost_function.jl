
mutable struct CostFunction{T,nq,nu,nc,nb}
    H::Int
    Qq::Vector{Diagonal{T,SizedArray{Tuple{nq},T,1,1}}}
    Qu::Vector{Diagonal{T,SizedArray{Tuple{nu},T,1,1}}}
    Qγ::Vector{Diagonal{T,SizedArray{Tuple{nc},T,1,1}}}
    Qb::Vector{Diagonal{T,SizedArray{Tuple{nb},T,1,1}}}
end

function CostFunction(H::Int, dim::Dimensions;
    Qq::Vector{Diagonal{T,SizedArray{Tuple{nq},T,1,1}}}=fill(Diagonal(zeros(SizedVector{dim.q})), H),
    Qu::Vector{Diagonal{T,SizedArray{Tuple{nu},T,1,1}}}=fill(Diagonal(zeros(SizedVector{dim.u})), H),
    Qγ::Vector{Diagonal{T,SizedArray{Tuple{nc},T,1,1}}}=fill(Diagonal(zeros(SizedVector{dim.c})), H),
    Qb::Vector{Diagonal{T,SizedArray{Tuple{nb},T,1,1}}}=fill(Diagonal(zeros(SizedVector{dim.b})), H),
    ) where {T,nq,nu,nc,nb}
    return CostFunction{T,nq,nu,nc,nb}(H, Qq, Qu, Qγ, Qb)
end

# function CostFunction(H::Int, dim::Dimensions;
#     Qq::Diagonal{T,SizedArray{Tuple{nq},T,1,1}}=Diagonal(zeros(SizedVector{dim.q})),
#     Qu::Diagonal{T,SizedArray{Tuple{nu},T,1,1}}=Diagonal(zeros(SizedVector{dim.u})),
#     Qγ::Diagonal{T,SizedArray{Tuple{nc},T,1,1}}=Diagonal(zeros(SizedVector{dim.c})),
#     Qb::Diagonal{T,SizedArray{Tuple{nb},T,1,1}}=Diagonal(zeros(SizedVector{dim.b})),
#     ) where {T,nq,nu,nc,nb}
#     return CostFunction{T,nq,nu,nc,nb}(H, fill(Qq,H), fill(Qu,H), fill(Qγ,H), fill(Qb,H))
# end
