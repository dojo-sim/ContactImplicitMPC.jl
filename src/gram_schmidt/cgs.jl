mutable struct CGSData{n,T} <: GSData{n,T}
    n::Int
    A::SparseMatrixCSC{T,Int}
    qs::SVector{n,Vector{T}}
    rs::Vector{Array{T,1}}
    qi_expr::Vector{Any}
    r_expr::Vector{Any}
    qu_expr::Any
    re_expr::Any
    q_expr::Any
end

function CGSData!(A::AbstractMatrix{T}, n::Int) where {T}
    qi_expr = []
    r_expr  = []
    for j = 1:n
        # Initialize qu
        @variables Af[1:nnz(A)]
        Ae = similar(A, eltype(Af))
        Ae.nzval .= Af
        a = Ae[:,j]
        qu = zeros(n) + Vector(a)
        push!(qi_expr, eval(build_function(simplify.(qu), Af)[2]))

        # Compute rk,j
        @variables Af[1:nnz(A)]
        Ae = similar(A, eltype(Af))
        Ae.nzval .= Af
        a = Ae[:,j]
        @variables qk[1:n]
        r = [transpose(qk)*a]
        push!(r_expr, eval(build_function(simplify.(r), qk, Af)[2]))
    end

    # Compute q_unnormalized
    @variables qu[1:n]
    @variables qk[1:n]
    @variables r[1:1]
    qu_ = qu - r[1]*qk
    qu_expr = eval(build_function(simplify.(qu_), qu, qk, r)[2])

    # Compute rj,j
    @variables qu[1:n]
    re = [sqrt(transpose(qu)*qu)]
    re_expr = eval(build_function(simplify.(re), qu)[2])

    # Compute q_normalized
    @variables qu[1:n]
    @variables re[1:1]
    q = qu ./ re[1]
    q_expr = eval(build_function(simplify.(q), qu, re)[2])

    # return qi_expr, r_expr, qu_expr, re_expr, q_expr
    qs = SVector{n,Vector{T}}([zeros(n) for i=1:n])
    rs = [zeros(T,1) for i=1:Int((n+1)*n/2)]
    return CGSData{n,T}(n, A, qs, rs, qi_expr, r_expr, qu_expr, re_expr, q_expr)
end

function cgs!(cgs_data::CGSData{n,T}, A::SparseMatrixCSC{T,Int}) where {n,T}
    # Unpack
    qs = cgs_data.qs
    rs = cgs_data.rs
    qi_expr = cgs_data.qi_expr
    r_expr = cgs_data.r_expr
    qu_expr = cgs_data.qu_expr
    re_expr = cgs_data.re_expr
    q_expr = cgs_data.q_expr

    off = 1
    for j = 1:n
        # qi
        qi_expr[j](qs[j], A.nzval)
        for k in eachindex((1:j-1))
            # rk
            r_expr[j](rs[off], qs[k], A.nzval)
            # qu
            qu_expr(qs[j], qs[j], qs[k], rs[off])
            off += 1
        end
        # re
        re_expr(rs[off], qs[j])
        # q
        q_expr(qs[j], qs[j], rs[off])
        off += 1
    end
    return nothing
end
