mutable struct DMGSData{n,T} <: GSData{n,T}
    n::Int
    A::Matrix{T}
    qs::SVector{n,Vector{T}}
    rs::Vector{Array{T,1}}
end

function DMGSData!(A::AbstractMatrix{T}, n::Int) where {T}
    qs = SVector{n,Vector{T}}([zeros(n) for i=1:n])
    rs = [zeros(T,1) for i=1:Int((n+1)*n/2)]
    return DMGSData{n,T}(n, A, qs, rs)
end

function dmgs!(mgs_data::DMGSData{n,T}, A::Matrix{T}) where {n,T}
    # Unpack
    qs = mgs_data.qs
    rs = mgs_data.rs

    off = 1
    for j in eachindex((1:n))
        # qi
        @inbounds @. @views qs[j] .= A[:,j]
        for k in eachindex((1:j-1))
            # rk
            # @inbounds @views rs[off] .= transpose(qs[j])*qs[k]
            @inbounds @views rs[off] .= LinearAlgebra.BLAS.dot(qs[j], qs[k])
            # qu
            @inbounds @. @views qs[j] .-= qs[k] .* rs[off]
            off += 1
        end
        # re
        # @inbounds @views rs[off] .= norm(qs[j],2)
        @inbounds @views rs[off] .= LinearAlgebra.BLAS.nrm2(qs[j])
        # q
        @inbounds @. @views qs[j] ./= rs[off]
        off += 1
    end
    return nothing
end


const ContactControl = Main



#
# T = Float64
# n = 80
# Random.seed!(100)
# A = sprand(n, n, 0.16)
# while rank(A) < n
#     A = sprand(n, n, 0.16)
# end
# Ad = Matrix(A)
# dgs_data = DMGSData!(Ad, n)
#
# @benchmark dmgs!(dgs_data, Ad)
# @time dmgs!(dgs_data, Ad)
# @test norm(A - hcat(dgs_data.qs...) * triangularize(dgs_data.rs, n), Inf) < 1e-10
# @benchmark qr_solve!(dgs_data, c, b)
#
# gs_data = MGSData!(A, n)
# @benchmark mgs!(gs_data, A)
# @test norm(A - hcat(gs_data.qs...) * triangularize(gs_data.rs, n), Inf) < 1e-10
#
# b = rand(SizedVector{n,T})
# c = zeros(SizedVector{n,T})
# qr_solve!(gs_data, c, b)
# @benchmark x = A1 \ b
# x = A1 \ b
# @test norm(c - x, Inf) < 1e-10
# @benchmark lu!($Ad, check = false)
# b = rand(T,n)
# c = zeros(T,n)
# x = zeros(T,n)
# @benchmark ldiv!(x, lu!(Ad, check = false), b)
#
#
#
# function ffff(M::Int)
#     for k = 1:M
#         dmgs!(dgs_data, Ad)
#     end
#     return nothing
# end
#
# @profiler ffff(10000)
