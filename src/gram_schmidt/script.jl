using Pkg
# Pkg.activate(joinpath(@__DIR__, "../src/mgs/"))

using Symbolics
using LinearAlgebra
using Test
using SparseArrays
using Random
using BenchmarkTools
using Plots
using StaticArrays

include("gram_schmidt.jl")
include("cgs.jl")
include("mgs.jl")

################################################################################
# CGS
################################################################################
# Generate a sparse random non-singular matrix
T = Float64
n = 20
Random.seed!(100)
A = sprand(n,n,0.16)
while rank(A) < n
    A = sprand(n,n,0.16)
end

# Generate the code to handle matrices with A's sparsity pattern
cgs_data = CGSData(A,n)

# Compute the Q and R matrices
@time cgs!(cgs_data, A)
@benchmark cgs!(cgs_data, A)
@test norm(A - hcat(cgs_data.qs...)*triangularize(cgs_data.rs,n), Inf) < 1e-10

@variables q[1:n]
d = [2.0 * q[1:2]; -1.0]
d1 = simplify.(d ./ sqrt.(sum(d.^2, dims = 1)))





d2 = simplify.(d ./ sqrt(transpose(d) * d))
# Test it on a different matrix with the same sparsity pattern
A1 = deepcopy(A)
A1.nzval .+= rand(nnz(A))
@time cgs!(cgs_data, A1)
@test norm(A1 - hcat(cgs_data.qs...)*triangularize(cgs_data.rs,n), Inf) < 1e-10


# Execute the simple backslash solve
b = rand(SizedVector{n,T})
c = zeros(SizedVector{n,T})
@benchmark qr_solve!(cgs_data, c, b)
x = A1\b
@test norm(c - x, Inf) < 1e-10


function profile(M)
    for k = 1:M
        # cgs!(cgs_data, A)
        qr_solve!(cgs_data, c, b)
    end
    return nothing
end
@profiler profile(10000)


################################################################################
# MGS
################################################################################
# Generate a sparse random non-singular matrix
T = Float64
n = 20
Random.seed!(100)
A = sprand(n,n,0.16)
while rank(A) < n
    A = sprand(n,n,0.16)
end

# Generate the code to handle matrices with A's sparsity pattern
mgs_data = MGSData!(A,n)

# Compute the Q and R matrices
@time mgs!(mgs_data, A)
@benchmark mgs!(mgs_data, A)
@test norm(A - hcat(mgs_data.qs...)*triangularize(mgs_data.rs,n), Inf) < 1e-10

# Test it on a different matrix with the same sparsity pattern
A1 = deepcopy(A)
A1.nzval .+= rand(nnz(A))
@time mgs!(mgs_data, A1)
@test norm(A1 - hcat(mgs_data.qs...)*triangularize(mgs_data.rs,n), Inf) < 1e-10

# Execute the simple backslash solve
b = rand(SizedVector{n,T})
c = zeros(SizedVector{n,T})
@benchmark qr_solve!(mgs_data, c, b)
x = A1\b
@test norm(c - x, Inf) < 1e-10


function profile(M)
    for k = 1:M
        # mgs!(mgs_data, A)
        qr_solve!(mgs_data, c, b)
    end
    return nothing
end
@profiler profile(10000)



################################################################################
# Test with large condition number
################################################################################

T = Float64
n = 43
Random.seed!(100)
A = sprand(n,n,0.16)
while rank(A) < n
    A = sprand(n,n,0.16)
end

# Generate the code to handle matrices with A's sparsity pattern
mgs_data = MGSData!(A,n)

# Compute the Q and R matrices
@time mgs!(mgs_data, A)
@benchmark mgs!(mgs_data, A)
@test norm(A - hcat(mgs_data.qs...)*triangularize(mgs_data.rs,n), Inf) < 1e-10

# Test it on a different matrix with the same sparsity pattern
A1 = deepcopy(A)
A1.nzval .= exp.(20*rand(nnz(A)))
@show cnd = LinearAlgebra.cond(Matrix(A1))
@time mgs!(mgs_data, A1)
@test norm(A1 - hcat(mgs_data.qs...)*triangularize(mgs_data.rs,n), Inf) < 1e-15*cnd

# Execute the simple backslash solve
b = rand(SizedVector{n,T})
c = zeros(SizedVector{n,T})
@benchmark qr_solve!(mgs_data, c, b)
@benchmark $x = $A1\$b
x = A1\b
@test norm(c - x, Inf) < 1e-15*cnd
norm(A1*c - b, 1) / norm(A1*x - b, 1)
