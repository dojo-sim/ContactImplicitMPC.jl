using Pkg
Pkg.activate(joinpath(@__DIR__, "../src/mgs/"))

using Symbolics
using LinearAlgebra
using Test
using SparseArrays
using Random
using BenchmarkTools
using Plots

function mgs!(A::AbstractMatrix{T}, R::UpperTriangular{T}) where {T}
    n = size(A)[1]
    @assert n == size(A)[2] == size(R)[1] == size(R)[2]
    for j = 1:n
        for k = 1:j-1
            R[k,j] = A[:,k]'*A[:,j]
            A[:,j] = A[:,j] - R[k,j]*A[:,k]
        end
        R[j,j] = norm(A[:,j], 2)
        A[:,j] ./= R[j,j]
    end
    return A,R
end

function mgs_flat!(A::AbstractMatrix{T}, R::AbstractVector{T}) where {T}
    for j = 1:n
        for k = 1:j-1
            R[triu_ind(k,j,n)] = A[:,k]'*A[:,j]
            A[:,j] = A[:,j] - R[triu_ind(k,j,n)]*A[:,k]
        end
        R[triu_ind(j,j,n)] = sqrt(transpose(A[:,j]) * A[:,j])
        A[:,j] ./= R[triu_ind(j,j,n)]
    end
    return A,R
end

function mgs_sparse!(A::AbstractMatrix{T}, R::AbstractVector{T}) where {T}
    Q = copy(A)
    for j = 1:n
        Q[:,j] = A[:,j]
        for k = 1:j-1
            R[triu_ind(k,j,n)] = A[:,k]'*q
            A[:,j] .= A[:,j] - R[triu_ind(k,j,n)]*A[:,k]
            A[:,j] = simplify(A[:,j])
        end
        R[triu_ind(j,j,n)] = sqrt(transpose(A[:,j]) * A[:,j])
        A[:,j] ./= R[triu_ind(j,j,n)]
        A[:,j] = simplify(A[:,j])
    end
    return A,R
end

function cgs_sparse!(A::AbstractMatrix{T}, Q::AbstractMatrix{T},
    R::AbstractVector{T}) where {T}
    for j = 1:n
        Q[:,j] = A[:,j]
        for k = 1:j-1
            R[triu_ind(k,j,n)] = Q[:,k]'*A[:,j]
            Q[:,j] .= Q[:,j] - R[triu_ind(k,j,n)]*Q[:,k]
            Q[:,j] = simplify(Q[:,j])
        end
        R[triu_ind(j,j,n)] = sqrt(transpose(Q[:,j]) * Q[:,j])
        Q[:,j] ./= R[triu_ind(j,j,n)]
        Q[:,j] = simplify(Q[:,j])
    end
    return A,Q,R
end

function triu_ind(i::Int, j::Int, n::Int)
    @assert i <= j
    return Int((n+n-i+1)*i/2) -n + j
end

function triangularize(Rf::AbstractVector{T}, n::Int) where {T}
    R = zeros(n,n)
    for i = 1:n
        for j = i:n
            R[i,j] = Rf[triu_ind(i,j,n)]
        end
    end
    return R
end

function triangularize_vect(rs)
    n = length(rs)
    R = zeros(n,n)
    for j = 1:n
        R[1:j,j] .= rs[j]
    end
    return UpperTriangular(R)
end

function test_triu(n::Int)
    R = UpperTriangular(zeros(n,n))
    for i = 1:n
        for j = i:n
            R[i,j] = triu_ind(i,j,n)
        end
    end
    return R
end

function mgs_code_gen(n::Int, m::Int, a::SparseMatrixCSC)
    @variables Ae[1:m]
    Q = zeros(eltype(Ae), n, n)
    A = similar(a, eltype(Ae))
    A.nzval .= Ae
    A += zeros(n,n)
    @variables R[1:Int((n+1)*n/2)]
    A,Q,R = cgs_sparse!(A,Q,R)

    Q_expr = build_function(Q, Ae)[2]
    R_expr = build_function(R, Ae)[2]
    return Q_expr, R_expr
end



n = 2
Random.seed!(1)
a = sprand(n,n,0.16)
while rank(a) < n
    a = sprand(n,n,0.16)
end
rank(a) == n
Q_expr, R_expr = mgs_code_gen(n,nnz(a),a)

Q_fast! = eval(Q_expr)
R_fast! = eval(R_expr)

q = zeros(n,n)
@time Q_fast!(q, a.nzval)
qs = deepcopy(sparse(q))
fill!(qs, 0.0)
@time Q_fast!(qs, a.nzval)
r = zeros(Int((n+1)*n/2))
@time R_fast!(r, a.nzval)
rs = triangularize(r, n)


b = rand(n)
x0 = a\b

x1 = rs\(qs'*b)

norm(a*x0 - b, Inf)
norm(a*x1 - b, Inf)
norm(qs*rs*x1 - b, Inf)

qs*rs - Matrix(a)
qs*rs



at, qt, rt = cgs_sparse!(a, zeros(n,n), zeros(Int(n*(n+1)/2)))
at - qt*triangularize(rt,n)






function mgs_expr!(A::AbstractMatrix{T}, R::AbstractVector{T}) where {T}
    for j = 1:n
        for k = 1:j-1
            R[triu_ind(k,j,n)] = A[:,k]'*A[:,j]
            A[:,j] .= A[:,j] - R[triu_ind(k,j,n)]*A[:,k]
            A[:,j] = simplify(A[:,j])
        end
        R[triu_ind(j,j,n)] = sqrt(transpose(A[:,j]) * A[:,j])
        A[:,j] ./= R[triu_ind(j,j,n)]
        A[:,j] = simplify(A[:,j])
    end
    return A,R
end



function mgs_best!(A::AbstractMatrix{T}, R::AbstractVector{T}, n::Int) where {T}
    q_expr = []
    r_expr = []
    for j = 1:n
        # # j = 1
        # @variables a[1:n]
        # q = a
        # r = zeros(eltype(q), 1)
        # r[end] = sqrt(transpose(q)*q)
        # q = q ./ r[end]
        #
        # q_expr = build_function(q, a)[2]
        # r_expr = build_function(r, a)[2]
        #
        # # j = 2
        # @variables a[1:n]
        # @variables q1[1:n]
        # q = a
        # r = zeros(eltype(a),2)
        # r[1] = transpose(q1)*q
        # q = q - r[1]*q1
        # r[end] = sqrt(transpose(q)*q)
        # q = q ./ r[end]

        # j = j
        @variables a[1:n]
        @variables qp[1:n*(j-1)]
        qs = [qp[(i-1)*n .+ (1:n)] for i=1:j-1]
        q = a
        r = zeros(eltype(a), j)
        for i = 1:j-1
            r[i] = transpose(qs[i])*q
            q = q - r[i]*qs[i]
        end
        r[end] = sqrt(transpose(q)*q)
        q = q ./ r[end]

        push!(q_expr, eval(build_function(q, a, qp)[2]))
        push!(r_expr, eval(build_function(r, a, qp)[2]))
    end
    return q_expr, r_expr
end



function cgs_best!(A::AbstractMatrix{T}, R::AbstractVector{T}, n::Int) where {T}
    q_expr = []
    qu_expr = []
    r_expr = []
    re_expr = []
    for j = 1:n
        # # j = 1
        # @variables a[1:n]
        # q = a
        # r = zeros(eltype(q), 1)
        # r[end] = sqrt(transpose(q)*q)
        # q = q ./ r[end]
        #
        # q_expr = build_function(q, a)[2]
        # r_expr = build_function(r, a)[2]
        #
        # # j = 2
        # @variables a[1:n]
        # @variables q1[1:n]
        # q = a
        # r = zeros(eltype(a),2)
        # r[1] = transpose(q1)*q
        # q = q - r[1]*q1
        # r[end] = sqrt(transpose(q)*q)
        # q = q ./ r[end]

        # j = j
        add_expr!(q_expr, qu_expr, r_expr, re_expr, j, n)
    end
    return q_expr, qu_expr, r_expr, re_expr
end

function add_expr!(q_expr, qu_expr, r_expr, re_expr, j::Int, n::Int)
    @variables a[1:n]
    @variables qp[1:n*(j-1)]
    qs = [qp[(i-1)*n .+ (1:n)] for i=1:j-1]
    q = a
    r = zeros(eltype(a), j-1)
    for i = 1:j-1
        r[i] = transpose(qs[i])*q
    end
    push!(r_expr, eval(build_function(r, a, qp)[2]))


    @variables a[1:n]
    @variables qp[1:n*(j-1)]
    @variables r[1:j-1]
    qs = [qp[(i-1)*n .+ (1:n)] for i=1:j-1]
    qu = a
    for i = 1:j-1
        qu = qu - r[i]*qs[i]
    end
    re = [sqrt(transpose(qu)*qu)]

    push!(qu_expr, eval(build_function(qu, a, qp, r)[2]))
    push!(re_expr, eval(build_function(re, a, qp, r)[2]))

    # Final
    @variables qu[1:n]
    @variables re[1:1]
    q = qu ./ re[1]
    push!(q_expr, eval(build_function(q, qu, re)[2]))

    return nothing
end


n = 43
A = sprand(n,n,0.16)
while rank(A) < n
    A = sprand(n,n,0.16)
end
Q = deepcopy(A)
b = rand(n)
A\b
R = zeros(Int((n+1)*n/2))
q_expr, qu_expr, r_expr, re_expr = cgs_best!(A,R,n)
@show q_expr[1]
@show q_expr[2]
@show q_expr[3]
@show q_expr[n]

function qr_fast(q_expr, qu_expr, r_expr, re_expr, A, qp, rs, n::Int)
    for j = 1:n
        # r
        qs = [qp[(i-1)*n .+ (1:n)] for i=1:j-1]
        qpj = qp[1:n*(j-1)]
        fr! = r_expr[j]
        rsj = rs[j]
        fr!(rsj, A[:,j], qpj)
        rs[j] = rsj

        # qu, re
        qj = qp[(j-1)*n .+ (1:n)]
        fqu! = qu_expr[j]
        fqu!(qj, A[:,j], qpj, rs[j][1:end-1])
        qp[(j-1)*n .+ (1:n)] = qj

        re = rs[j][end:end]
        fre! = re_expr[j]
        fre!(re, A[:,j], qpj, rs[j][1:end-1])
        rs[j][end:end] .= re

        # q
        qj = qp[(j-1)*n .+ (1:n)]
        fq! = q_expr[j]
        fq!(qj, qj, rs[j][end:end])
        qp[(j-1)*n .+ (1:n)] = qj

    end
    return qp, rs
end

qp = zeros(n*n)
rs = [zeros(i) for i=1:n]
@btime qr_fast(q_expr, qu_expr, r_expr, re_expr, A, qp, rs, n)
qp
rs
triangularize_vect(rs)


Matrix(A) - reshape(qp, (n,n))*triangularize_vect(rs)






plot([nnz(Q[:,i]) for i=1:n])
plot!([nnz(A[:,i]) for i=1:n])
sparse(R)





a = 11111
a = 11111
a = 11111




# Basic test
n = 43
R = UpperTriangular(Float64.(rand(-n:n,n,n)))
A_ = rand(-n:n, n,n)
A = Float64.(A_'*A_)
b = Float64.(rand(-n:n,n))

Q = deepcopy(A)
mgs!(Q,R)
@test norm(A - Q*R, Inf) < 1e-8

Q = deepcopy(A)
Rf = zeros(Int((n+1)*n/2))
mgs_flat!(Q,Rf)
R = triangularize(Rf,n)
@test norm(A - Q*R, Inf) < 1e-8


@benchmark A\b # 160 μs
Q = deepcopy(A)
R = UpperTriangular(zeros(n,n))
@benchmark mgs!(Q,R) # 600μs
@benchmark R\(Q'b) # 1.5μs
