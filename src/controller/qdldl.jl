module QDLDL

export qdldl, \, solve, solve!, update_diagonal!, positive_inertia

using AMD, SparseArrays
using LinearAlgebra: istriu, triu, Diagonal

const QDLDL_UNKNOWN = -1;
const QDLDL_USED   = true;
const QDLDL_UNUSED = false;


struct QDLDLWorkspace{Tf<:AbstractFloat,Ti<:Integer}

    #internal workspace data
    etree::Vector{Ti}
      Lnz::Vector{Ti}
    iwork::Vector{Ti}
    bwork::Vector{Bool}
    fwork::Vector{Tf}

    #L matrix row indices and data
    Ln::Int         #always Int since SparseMatrixCSC does it this way
    Lp::Vector{Ti}
    Li::Vector{Ti}
    Lx::Vector{Tf}

    #D and its inverse
    D::Vector{Tf}
    Dinv::Vector{Tf}

    #number of positive values in D
    positive_inertia::Base.RefValue{Ti}

    #The upper triangular matrix factorisation target
    triuA::SparseMatrixCSC{Tf,Ti}
end

function QDLDLWorkspace(triuA::SparseMatrixCSC{Tf,Ti}) where {Tf<:AbstractFloat,Ti<:Integer}

    #A should be an upper triangular matrix input

    etree  = Vector{Ti}(undef,triuA.n)
    Lnz    = Vector{Ti}(undef,triuA.n)
    iwork  = Vector{Ti}(undef,triuA.n*3)
    bwork  = Vector{Bool}(undef,triuA.n)
    fwork  = Vector{Tf}(undef,triuA.n)

    #compute elimination gree using QDLDL converted code
    sumLnz = QDLDL_etree!(triuA.n,triuA.colptr,triuA.rowval,iwork,Lnz,etree)

    if(sumLnz < 0)
        error("Input matrix is not upper triangular or has an empty column")
    end

    #allocate space for the L matrix row indices and data
    Ln = triuA.n
    Lp = Vector{Ti}(undef,triuA.n + 1)
    Li = Vector{Ti}(undef,sumLnz)
    Lx = Vector{Tf}(undef,sumLnz)

    #allocate for D and D inverse
    D  = Vector{Tf}(undef,triuA.n)
    Dinv = Vector{Tf}(undef,triuA.n)

    #allocate for positive inertia count.  -1 to
    #start since we haven't counted anything yet
    positive_inertia = Base.RefValue{Ti}(-1)

    QDLDLWorkspace(etree,Lnz,iwork,bwork,fwork,Ln,Lp,Li,Lx,D,Dinv,positive_inertia,triuA)

end

struct QDLDLFactorisation{Tf<:AbstractFloat,Ti<:Integer}

    #permutation vector (nothing if no permutation)
    perm::Union{Nothing,Vector{Ti}}
    #inverse permutation (nothing if no permutation)
    iperm::Union{Nothing,Vector{Ti}}
    #lower triangular factor
    L::SparseMatrixCSC{Tf,Ti}
    #Inverse of D matrix in ldl
    Dinv::Diagonal{Tf,Vector{Tf}}
    #workspace data
    workspace::QDLDLWorkspace{Tf,Ti}
    #is it logical factorisation only?
    logical::Bool
end




# Usage :
# qdldl(A) uses the default AMD ordering
# qdldl(A,perm=p) uses a caller specified ordering
# qdldl(A,perm = nothing) factors without reordering
#
# qdldl(A,logical=true) produces a logical factorisation only

function qdldl(A::SparseMatrixCSC{Tv,Ti};
               perm::Union{Array{Ti},Nothing}=amd(A),
               logical::Bool=false
              ) where {Tv<:AbstractFloat, Ti<:Integer}

    #store the inverse permutation to enable matrix updates
    iperm = perm == nothing ? nothing : invperm(perm)

    #permute using symperm, producing a triu matrix to factor
    if perm != nothing
        A = permute_symmetric(A, iperm)  #returns an upper triangular matrix
    else
        if(!istriu(A))
            A = triu(A);
        end
    end

    #allocate workspace
    workspace = QDLDLWorkspace(A)

    #factor the matrix
    factor!(workspace,logical)

    #make user-friendly factors
    L = SparseMatrixCSC(workspace.Ln,
                        workspace.Ln,
                        workspace.Lp,
                        workspace.Li,
                        workspace.Lx)
    Dinv = Diagonal(workspace.Dinv)

    return QDLDLFactorisation(perm, iperm, L, Dinv, workspace, logical)

end

function positive_inertia(F::QDLDLFactorisation)
    F.workspace.positive_inertia[]
end


function update_diagonal!(F::QDLDLFactorisation,indices,scalarValue::Real)
    update_diagonal!(F,indices,[scalarValue])
end


function update_diagonal!(F::QDLDLFactorisation,indices,values)

    (length(values) != length(indices) && length(values) != 1 ) &&
        throw(DimensionMismatch("Index and value arrays must be the same size, or values must be a scalar."))

    triuA = F.workspace.triuA
    invp  = F.iperm
    nvals = length(values)

    #triuA is full rank and  upper triangular, so the diagonal element
    #in each column will always be the last nonzero
    for i in 1:length(indices)
         idx = invp[indices[i]]
         val = nvals == 1 ? values[1] : values[i]
         triuA.nzval[triuA.colptr[idx+1]-1] = val
    end

    #force a refactorisation
    refactor!(F)

end


function Base.:\(F::QDLDLFactorisation,b)
    return solve(F,b)
end


function refactor!(F::QDLDLFactorisation)
    factor!(F.workspace,F.logical)
end


function factor!(workspace::QDLDLWorkspace{Tf,Ti},logical) where {Tf<:AbstractFloat,Ti<:Integer}

    if(logical)
        workspace.Lx   .= 1
        workspace.D    .= 1
        workspace.Dinv .= 1
    end

    #factor using QDLDL converted code
    A = workspace.triuA
    posDCount = QDLDL_factor!(A.n,A.colptr,A.rowval,A.nzval,
                              workspace.Lp,
                              workspace.Li,
                              workspace.Lx,
                              workspace.D,
                              workspace.Dinv,
                              workspace.Lnz,
                              workspace.etree,
                              workspace.bwork,
                              workspace.iwork,
                              workspace.fwork,
                              logical)

    if(posDCount < 0)
        error("Zero entry in D (matrix is not quasidefinite)")
    end

    workspace.positive_inertia[] = posDCount

    return nothing

end


# Solves Ax = b using LDL factors for A.
# Returns x, preserving b
function solve(F::QDLDLFactorisation,b)
    x = copy(b)
    solve!(F,x)
    return x
end

# Solves Ax = b using LDL factors for A.
# Solves in place (x replaces b)
function solve!(F::QDLDLFactorisation,b)

    #bomb if logical factorisation only
    if F.logical
        error("Can't solve with logical factorisation only")
    end

    #permute b
    tmp = F.perm == nothing ? b : permute!(F.workspace.fwork,b,F.perm)

    QDLDL_solve!(F.workspace.Ln,
                 F.workspace.Lp,
                 F.workspace.Li,
                 F.workspace.Lx,
                 F.workspace.Dinv,
                 tmp)

    #inverse permutation
    b = F.perm == nothing ? tmp : ipermute!(b,F.workspace.fwork,F.perm)

    return nothing
end



# Compute the elimination tree for a quasidefinite matrix
# in compressed sparse column form.

function QDLDL_etree!(n,Ap,Ai,work,Lnz,etree)

    @inbounds for i = 1:n
        # zero out Lnz and work.  Set all etree values to unknown
        work[i]  = 0
        Lnz[i]   = 0
        etree[i] = QDLDL_UNKNOWN

        #Abort if A doesn't have at least one entry
        #one entry in every column
        if(Ap[i] == Ap[i+1])
            return -1
        end
    end

    @inbounds for j = 1:n
        work[j] = j
        @inbounds for p = Ap[j]:(Ap[j+1]-1)
            i = Ai[p]
            if(i > j)
                return -1
            end
            @inbounds while(work[i] != j)
                if(etree[i] == QDLDL_UNKNOWN)
                    etree[i] = j
                end
                Lnz[i] += 1        #nonzeros in this column
                work[i] = j
                i = etree[i]
            end
        end #end for p
    end

    #tally the total nonzeros
    sumLnz = sum(Lnz)

    return sumLnz
end




function QDLDL_factor!(n,Ap,Ai,Ax,Lp,Li,Lx,D,Dinv,Lnz,etree,bwork,iwork,fwork,logicalFactor)


    positiveValuesInD = 0

    #partition working memory into pieces
    yMarkers        = bwork
    yIdx            = view(iwork,      1:n)
    elimBuffer      = view(iwork,  (n+1):2*n)
    LNextSpaceInCol = view(iwork,(2*n+1):3*n)
    yVals           = fwork;


    Lp[1] = 1 #first column starts at index one / Julia is 1 indexed

    @inbounds for i = 1:n

        #compute L column indices
        Lp[i+1] = Lp[i] + Lnz[i]   #cumsum, total at the end

        # set all Yidx to be 'unused' initially
        #in each column of L, the next available space
        #to start is just the first space in the column
        yMarkers[i]  = QDLDL_UNUSED
        yVals[i]     = 0.0
        D[i]         = 0.0
        LNextSpaceInCol[i] = Lp[i]
    end

    if(!logicalFactor)
        # First element of the diagonal D.
        D[1]     = Ax[1]
        if(D[1] == 0.0) return -1 end
        if(D[1]  > 0.0) positiveValuesInD += 1 end
        Dinv[1] = 1/D[1];
    end

    #Start from 1 here. The upper LH corner is trivially 0
    #in L b/c we are only computing the subdiagonal elements
    @inbounds for k = 2:n

        #NB : For each k, we compute a solution to
        #y = L(0:(k-1),0:k-1))\b, where b is the kth
        #column of A that sits above the diagonal.
        #The solution y is then the kth row of L,
        #with an implied '1' at the diagonal entry.

        #number of nonzeros in this row of L
        nnzY = 0  #number of elements in this row

        #This loop determines where nonzeros
        #will go in the kth row of L, but doesn't
        #compute the actual values
        @inbounds for i = Ap[k]:(Ap[k+1]-1)

            bidx = Ai[i]   # we are working on this element of b

            #Initialize D[k] as the element of this column
            #corresponding to the diagonal place.  Don't use
            #this element as part of the elimination step
            #that computes the k^th row of L
            if(bidx == k)
                D[k] = Ax[i];
                continue
            end

            yVals[bidx] = Ax[i]   # initialise y(bidx) = b(bidx)

            # use the forward elimination tree to figure
            # out which elements must be eliminated after
            # this element of b
            nextIdx = bidx

            if(yMarkers[nextIdx] == QDLDL_UNUSED)  #this y term not already visited

                yMarkers[nextIdx] = QDLDL_USED     #I touched this one
                elimBuffer[1]     = nextIdx  # It goes at the start of the current list
                nnzE              = 1         #length of unvisited elimination path from here

                nextIdx = etree[bidx];

                @inbounds while(nextIdx != QDLDL_UNKNOWN && nextIdx < k)
                    if(yMarkers[nextIdx] == QDLDL_USED) break; end

                    yMarkers[nextIdx] = QDLDL_USED;   #I touched this one
                    #NB: Julia is 1-indexed, so I increment nnzE first here,
                    #no after writing into elimBuffer as in the C version
                    nnzE += 1                   #the list is one longer than before
                    elimBuffer[nnzE] = nextIdx; #It goes in the current list
                    nextIdx = etree[nextIdx];   #one step further along tree

                end #end while

                # now I put the buffered elimination list into
                # my current ordering in reverse order
                @inbounds while(nnzE != 0)
                    #NB: inc/dec reordered relative to C because
                    #the arrays are 1 indexed
                    nnzY += 1;
                    yIdx[nnzY] = elimBuffer[nnzE];
                    nnzE -= 1;
                end #end while
            end #end if

        end #end for i

        #This for loop places nonzeros values in the k^th row
        @inbounds for i = nnzY:-1:1

            #which column are we working on?
            cidx = yIdx[i]

            # loop along the elements in this
            # column of L and subtract to solve to y
            tmpIdx = LNextSpaceInCol[cidx];

            #don't compute Lx for logical factorisation
            #this is not implemented in the C version
            if(!logicalFactor)
                yVals_cidx = yVals[cidx]
                @inbounds for j = Lp[cidx]:(tmpIdx-1)
                    yVals[Li[j]] -= Lx[j]*yVals_cidx
                end

                #Now I have the cidx^th element of y = L\b.
                #so compute the corresponding element of
                #this row of L and put it into the right place
                Lx[tmpIdx] = yVals_cidx *Dinv[cidx]

                #D[k] -= yVals[cidx]*yVals[cidx]*Dinv[cidx];
                D[k] -= yVals_cidx*Lx[tmpIdx]
            end

            #also record which row it went into
            Li[tmpIdx] = k

            LNextSpaceInCol[cidx] += 1

            #reset the yvalues and indices back to zero and QDLDL_UNUSED
            #once I'm done with them
            yVals[cidx]     = 0.0
            yMarkers[cidx]  = QDLDL_UNUSED;

        end #end for i

        #Maintain a count of the positive entries
        #in D.  If we hit a zero, we can't factor
        #this matrix, so abort
        if(D[k] == 0.0) return -1 end
        if(D[k]  > 0.0) positiveValuesInD += 1 end

        #compute the inverse of the diagonal
        Dinv[k]= 1/D[k]

    end #end for k

    return positiveValuesInD

end

# Solves (L+I)x = b, with x replacing b
function QDLDL_Lsolve!(n,Lp,Li,Lx,x)

    @inbounds for i = 1:n
        @inbounds for j = Lp[i]: (Lp[i+1]-1)
            x[Li[j]] -= Lx[j]*x[i];
        end
    end
    return nothing
end


# Solves (L+I)'x = b, with x replacing b
function QDLDL_Ltsolve!(n,Lp,Li,Lx,x)

    @inbounds for i = n:-1:1
        @inbounds for j = Lp[i]:(Lp[i+1]-1)
            x[i] -= Lx[j]*x[Li[j]]
        end
    end
    return nothing
end

# Solves Ax = b where A has given LDL factors,
# with x replacing b
function QDLDL_solve!(n,Lp,Li,Lx,Dinv,b)

    QDLDL_Lsolve!(n,Lp,Li,Lx,b)
    b .*= Dinv;
    QDLDL_Ltsolve!(n,Lp,Li,Lx,b)

end



# internal permutation and inverse permutation
# functions that require no memory allocations
function permute!(x,b,p)
  @inbounds for j = 1:length(x)
      x[j] = b[p[j]];
  end
  return x
end

function ipermute!(x,b,p)
 @inbounds for j = 1:length(x)
     x[p[j]] = b[j];
 end
 return x
end


"Given a sparse symmetric matrix `A` (with only upper triangular entries), return permuted sparse symmetric matrix `P` (only upper triangular) given the inverse permutation vector `iperm`."
function permute_symmetric(A::SparseMatrixCSC{Tv, Ti}, iperm::AbstractVector{Ti},
    Pr::AbstractVector{Ti} = zeros(Ti, nnz(A)),
    Pc::AbstractVector{Ti} = zeros(Ti, size(A, 1) + 1),
    Pv::AbstractVector{Tv} = zeros(Tv, nnz(A)) ) where {Tv <: AbstractFloat, Ti <: Integer}

    # perform a number of argument checks
    m, n = size(A)
    m != n && throw(DimensionMismatch("Matrix A must be sparse and square"))

    isperm(iperm) || throw(ArgumentError("pinv must be a permutation"))

    if n != length(iperm)
        throw(DimensionMismatch("Dimensions of sparse matrix A must equal the length of iperm, $((m,n)) != $(iperm)"))
    end
    return _permute_symmetric(A, iperm, Pr, Pc, Pv)
end

# the main function without extra argument checks
# following the book: Timothy Davis - Direct Methods for Sparse Linear Systems
function _permute_symmetric(A::SparseMatrixCSC{Tv, Ti}, iperm::AbstractVector{Ti}, Pr::AbstractVector{Ti}, Pc::AbstractVector{Ti}, Pv::AbstractVector{Tv}) where {Tv <: AbstractFloat, Ti <: Integer}
    # 1. count number of entries that each column of P will have
    n = size(A, 2)
    num_entries = zeros(Ti, n)
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


end #end module
