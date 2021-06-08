"""
    QDLDL inplace functionality
"""
mutable struct LDLSolver{Tf<:AbstractFloat,Ti<:Integer} <: LinearSolver
    # QDLDL Factorization
    F::QDLDL.QDLDLFactorisation{Tf,Ti}
    # Allocate memory
    Pr::Vector{Ti}
    Pc::Vector{Ti}
    Pv::Vector{Tf}
    num_entries::Vector{Ti}
end

function LDLSolver(A::SparseMatrixCSC{Tv,Ti}, F::QDLDL.QDLDLFactorisation{Tv,Ti}) where {Tv<:AbstractFloat,Ti<:Integer}
    Pr = zeros(Ti, nnz(A))
    Pc = zeros(Ti, size(A, 1) + 1)
    Pv = zeros(Tv, nnz(A))
    num_entries = zeros(Ti, size(A, 2))
    return LDLSolver{Tv,Ti}(F, Pr, Pc, Pv, num_entries)
end

function qdldl!(A::SparseMatrixCSC{Tv,Ti},
                s::LDLSolver{Tv,Ti};
              ) where {Tv<:AbstractFloat, Ti<:Integer}
    # Reset the pre-allocated fields
    s.Pr .= 0
    s.Pc .= 0
    s.Pv .= 0.0
    s.num_entries .= 0

    # Triangularize the matrix with the allocation-free method.
    A = permute_symmetricAF(A, s.F.iperm, s.Pr, s.Pc, s.Pv, s.num_entries)  #returns an upper triangular matrix

    # Update the workspace, triuA is the only field we need to update
    s.F.workspace.triuA.nzval .= A.nzval

    #factor the matrix
    QDLDL.factor!(s.F.workspace, s.F.logical)

    return nothing
end

function permute_symmetricAF(A::SparseMatrixCSC{Tv, Ti}, iperm::AbstractVector{Ti},
    Pr::AbstractVector{Ti}, Pc::AbstractVector{Ti}, Pv::AbstractVector{Tv},
    num_entries::AbstractVector{Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    # 1. count number of entries that each column of P will have
    n = size(A, 2)
    # num_entries = zeros(Ti, n)
    Ar = A.rowval
    Ac = A.colptr
    Av = A.nzval
    # count the number of upper-triangle entries in columns of P, keeping in mind the row permutation
    for colA = 1:n
        colP = iperm[colA]
        # loop over entries of A in column A...
        for row_idx = Ac[colA]:Ac[colA+1]-1
            rowA = Ar[row_idx]
            rowP = iperm[rowA]
            # ...and check if entry is upper triangular
            if rowA <= colA
                # determine to which column the entry belongs after permutation
                col_idx = max(rowP, colP)
                num_entries[col_idx] += one(Ti)
            end
        end
    end
    # 2. calculate permuted Pc = P.colptr from number of entries
    Pc[1] = one(Ti)
    @inbounds for k = 1:n
        Pc[k + 1] = Pc[k] + num_entries[k]

        # reuse this vector memory to keep track of free entries in rowval
        num_entries[k] = Pc[k]
    end
    # use alias
    row_starts = num_entries

    # 3. permute the row entries and position of corresponding nzval
    for colA = 1:n
        colP = iperm[colA]
        # loop over rows of A and determine where each row entry of A should be stored
        for rowA_idx = Ac[colA]:Ac[colA+1]-1
            rowA = Ar[rowA_idx]
            # check if upper triangular
            if rowA <= colA
                rowP = iperm[rowA]
                # determine column to store the entry
                col_idx = max(colP, rowP)

                # find next free location in rowval (this results in unordered columns in the rowval)
                rowP_idx = row_starts[col_idx]

                # store rowval and nzval
                Pr[rowP_idx] = min(colP, rowP)
                Pv[rowP_idx] = Av[rowA_idx]

                # increment next free location
                row_starts[col_idx] += 1
            end
        end
    end
    P = SparseMatrixCSC{Tv, Ti}(n, n, Pc, Pr, Pv)
    # order row indices within P.rowcal[P.colptr[k]:P.colptr[k+1]-1]
    return (P')'
end


"""
    LDL solver
"""
function ldl_solver(A::SparseMatrixCSC{T,Int}) where T
    LDLSolver(A, qdldl(A))
end

ldl_solver(A::Array{T, 2}) where T = ldl_solver(sparse(A))

function linear_solve!(solver::LDLSolver{Tv,Ti}, x::Vector{Tv}, A::SparseMatrixCSC{Tv,Ti}, b::Vector{Tv}) where {Tv<:AbstractFloat,Ti<:Integer}
    qdldl!(A, solver) # factorize
    x .= b
    QDLDL.solve!(solver.F, x) # solve
end

function linear_solve!(solver::LDLSolver{Tv,Ti}, X::Matrix{Tv}, A::AbstractMatrix{Tv}, B::AbstractMatrix{Tv}) where {Tv<:AbstractFloat,Ti<:Integer}
    X .= A \ B # TODO: fix
end
