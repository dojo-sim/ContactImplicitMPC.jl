abstract type GSData{n,T}
end

function triangularize(rs::AbstractVector,n::Int)
    R = zeros(n,n)
    off = 1
    for j = 1:n
        for k = 1:j
            R[k,j] = rs[off][1]
            off += 1
        end
    end
    return UpperTriangular(R)
end

function triu_perm(k::Int, j::Int)
    ind = Int((j-1)*j/2)+k
    return ind::Int
end

function qr_solve!(gs_data::GSData{n,T}, x::SizedVector{n,T}, b::SizedVector{n,T}) where {n,T}
    qr_solve!(gs_data.qs, gs_data.rs, x, b)
end

function qr_solve!(qs::SVector{n,V}, rs::Vector{V}, x::SizedVector{n,T}, b::SizedVector{n,T}) where {n,T,V<:Vector{T}}
    for j = 1:n
        x[j] = transpose(qs[j])*b
    end
    for j = n:-1:1
        for k = j+1:n
            x[j] -= rs[triu_perm(j,k)][1]*x[k]
        end
        x[j] /= rs[triu_perm(j,j)][1]
    end
    return nothing
end
