"""
    LU solver
"""
function lu_solver(A::SparseMatrixCSC{T,Int}) where T
    A = Array(A) 
    lu_solver(A)
end

function linear_solve!(solver::LUSolver{T}, x::Vector{T}, A::SparseMatrixCSC{T,Int},
        b::Vector{T}; reg::T = 0.0, fact::Bool = true) where T
    linear_solve!(solver, x, Array(A), b, reg=reg, fact=fact)
end

## Sparse LU solver 






using SparseArrays, SuiteSparse, LinearAlgebra 
n = 20
S = sprand(n, n, 0.4)
b = rand(n) 

xsol1 = S \ b

F = lu(S)


lu!(F, S)
@benchmark lu!($F, $S)

@benchmark SparseArrays.getcolptr($S)
@benchmark SparseArrays.rowvals($S)

UMFPACK_INFO = 90
UMFPACK_CONTROL = 20
umf_ctrl = Vector{Float64}(undef, UMFPACK_CONTROL)
umf_info = Vector{Float64}(undef, UMFPACK_INFO)
tmp = Vector{Ptr{Cvoid}}(undef, 1)

function umfpack_numeric!(U, tmp, umf_ctrl, umf_info)
    status = ccall((:umfpack_di_numeric, :libumfpack), Int64,
                   (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Cvoid}, Ptr{Cvoid},
                    Ptr{Float64}, Ptr{Float64}),
                   U.colptr, U.rowval, U.nzval, U.symbolic, tmp,
                   umf_ctrl, umf_info)
    return nothing
end

@benchmark umfpack_numeric!($F, $tmp, $umf_ctrl, $umf_info)
F.colptr
typeof(F)
S
function my_lu!(F::SuiteSparse.UMFPACK.UmfpackLU{T,I}, S::SparseMatrixCSC{T,I}, tmp, umf_ctrl, umf_info) where {T,I} 
    F.m = size(S, 1)
    F.n = size(S, 2)
    F.colptr = SparseArrays.getcolptr(S)
    F.rowval = SparseArrays.rowvals(S)
    F.nzval .= S.nzval  
    umfpack_numeric!(F, tmp, umf_ctrl, umf_info)
end

my_lu!(F, S, tmp, umf_ctrl, umf_info)
xsol2 = F \ b
@benchmark my_lu!($F, $S, $tmp, $umf_ctrl, $umf_info)

# F.m = size(S, 1)
# F.n = size(S, 2)
# F.colptr .= SparseArrays.getcolptr(S)# zerobased ? copy(getcolptr(S)) : decrement(getcolptr(S))
# F.rowval .= SparseArrays.rowvals(S)#zerobased ? copy(rowvals(S)) : decrement(rowvals(S))
# F.nzval .= S.nzval#copy(nonzeros(S))
# umfpack_numeric!(F, reuse_numeric = false)


# umfpack_numeric!(F, reuse_numeric = false)
# check && (issuccess(F) || throw(LinearAlgebra.SingularException(0)))

function lu!(F::UmfpackLU, S::SparseMatrixCSC{<:UMFVTypes,<:UMFITypes}; check::Bool=true)
    zerobased = getcolptr(S)[1] == 0
    F.m = size(S, 1)
    F.n = size(S, 2)
    F.colptr = zerobased ? copy(getcolptr(S)) : decrement(getcolptr(S))
    F.rowval = zerobased ? copy(rowvals(S)) : decrement(rowvals(S))
    F.nzval = copy(nonzeros(S))

    umfpack_numeric!(F, reuse_numeric = false)
    check && (issuccess(F) || throw(LinearAlgebra.SingularException(0)))
    return F
end


using SparseArrays, SuiteSparse, LinearAlgebra, BenchmarkTools


const UMFPACK_INFO = 90
const UMFPACK_CONTROL = 20
const umf_ctrl = Vector{Float64}(undef, UMFPACK_CONTROL)
const umf_info = Vector{Float64}(undef, UMFPACK_INFO)
tmp = Vector{Ptr{Cvoid}}(undef, 1)

n = 350
S = sprand(n, n, 0.4)
b = rand(n) 

xsol1 = S \ b

F = lu(copy(S))

# @benchmark SuiteSparse.UMFPACK.lu!($F, $S)
function _umfpack_numeric!(U::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64}, tmp; reuse_numeric = true)
    # if (reuse_numeric && U.numeric != C_NULL) return U end
    # if U.symbolic == C_NULL umfpack_symbolic!(U) end
    
    status = ccall((:umfpack_dl_numeric, :libumfpack), Int64,
                   (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Cvoid}, Ptr{Cvoid},
                    Ptr{Float64}, Ptr{Float64}),
                   U.colptr, U.rowval, U.nzval, U.symbolic, tmp,
                   umf_ctrl, umf_info)
    return nothing
    # U.status = status
    # if status != SuiteSparse.UMFPACK.UMFPACK_WARNING_singular_matrix
    #     umferror(status)
    # end
    # U.numeric != C_NULL && SuiteSparse.UMFPACK.umfpack_free_numeric(U)
    # U.numeric = tmp[1]
    # return U
end

_umfpack_numeric!(F, tmp, reuse_numeric=false)

function _lu!(F::SuiteSparse.UMFPACK.UmfpackLU, S::SparseMatrixCSC{<:SuiteSparse.UMFPACK.UMFVTypes,<:SuiteSparse.UMFPACK.UMFITypes}, tmp; check::Bool=true)
    zerobased = SparseArrays.getcolptr(S)[1] == 0
    F.m = size(S, 1)
    F.n = size(S, 2)
    F.colptr .= SparseArrays.getcolptr(S)
    F.rowval .= SparseArrays.rowvals(S) 
    if !zerobased 
        F.colptr .-= 1 
        F.rowval .-= 1 
    end
    F.nzval .= nonzeros(S)

    _umfpack_numeric!(F, tmp, reuse_numeric = false)

    # SuiteSparse.UMFPACK.umfpack_numeric!(F, reuse_numeric = false)
    # check && (issuccess(F) || throw(LinearAlgebra.SingularException(0)))
    # return F
end

_lu!(F, S, tmp)
@benchmark _lu!($F, $S, $tmp)

xsol2 = F \ b
norm(xsol1 - xsol2)

xsol3 = zeros(n)
# ldiv!(xsol3, F, b)
# @benchmark ldiv!($xsol3, $F, $b)
# norm(xsol3 - xsol1)

function _solve!(x, lu, b)
    ccall((:umfpack_dl_solve, :libumfpack), Int64,
        (Int64, Ptr{Int64}, Ptr{Int64}, Ptr{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}, Ptr{Float64},
         Ptr{Float64}),
        0, lu.colptr, lu.rowval, lu.nzval,
        x, b, lu.numeric, umf_ctrl,
        umf_info)
    return x
end

_solve!(xsol3, F, b)

norm(xsol3 - xsol1, Inf)