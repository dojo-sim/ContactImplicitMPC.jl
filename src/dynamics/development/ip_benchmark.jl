using BenchmarkTools, InteractiveUtils

"""
    minimize   x' P x + q' x
    subject to    x >= 0
"""

n = 1000

_P = rand(n)
P = Diagonal(_P)
q = rand(n)
θ = [_P; q]
z = ones(2 * n)
r = zeros(2 * n)
rz = zeros(2 * n, 2 * n)
rθ = zeros(2 * n, 2 * n)
κ = 1.0

idx_ineq = collect(1:2 * n)

# residual
function _r!(r, z, θ, κ)
    x = z[1:1000]
    y = z[1001:2000]
    P = θ[1:1000]
    q = θ[1001:2000]

    r[1:1000] = 2.0 * Diagonal(P) * x + q - y
    r[1001:2000] = Diagonal(x) * y .- κ
    nothing
end

@variables r_sym[1:2000]
@variables z_sym[1:2000]
@variables θ_sym[1:2000]
@variables κ_sym[1:1]

parallel = false
_r!(r_sym, z_sym, θ_sym, κ_sym)
r_sym = simplify.(r_sym)
rf! = eval(ModelingToolkit.build_function(r_sym, z_sym, θ_sym, κ_sym,
    parallel = parallel)[2])
rz_exp = ModelingToolkit.sparsejacobian(r_sym, z_sym, simplify = true)
rθ_exp = ModelingToolkit.sparsejacobian(r_sym, θ_sym, simplify = true)
rz_sp = similar(rz_exp, Float64)
rθ_sp = similar(rθ_exp, Float64)
rzf! = eval(ModelingToolkit.build_function(rz_exp, z_sym, θ_sym, κ_sym,
    parallel = parallel)[2])
rθf! = eval(ModelingToolkit.build_function(rθ_exp, z_sym, θ_sym, κ_sym,
    parallel = parallel)[2])

# solver
ip = ContactControl.interior_point(2 * n, 2 * n, idx_ineq,
    r! = rf!, rz! = rzf!, rθ! = rθf!,
    rz = rz_sp,
    rθ = rθ_sp)

# options
ip.z .= copy(z)
ip.θ .= copy(θ)

rf!(ip.r, ip.z, ip.θ, 1.0)
rzf!(ip.rz, ip.z, ip.θ, 1.0)
rank(ip.rz)
_rz = Array(copy(ip.rz))
_rz_sp = copy(ip.rz)
@time ldiv!(ip.Δ, lu!(_rz), ip.r)
@time ldiv!(ip.Δ, lu(_rz_sp), ip.r)

lu(_rz_sp)
@time info = LAPACK.getrf!(_rz)
@time LAPACK.getrs!('N', info[1], info[2], ip.Δ)

fact = lu(ip.rz)
ldiv!(ip.Δ, fact, r)
# solve
@time status = ContactControl.interior_point!(ip, z, θ, opts = opts)

# test
@test status
@test norm(ip.r, Inf) < opts.r_tol
@test !ContactControl.inequality_check(ip.z, ip.idx_ineq)
@test ip.κ[1] < opts.κ_tol



A = sparse(rand(n, n))
B = A' * A
b = rand(n)
B \ b
F = lu(B)
sol = zeros(n)
ldiv!(sol, F, b)
sol

using SuiteSparse
struct LAsolve6{T}
    mat::SparseMatrixCSC{T,Int}
    fact#::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64}
    sol::Vector{T}
    b::Vector{T}
end

mm = LAsolve6(B, lu(B), zeros(n), rand(n))
mm.fact .= mm.fact
function solve!(x::LAsolve3)
    # x.fact .= lu(x.mat)
    ldiv!(x.sol, x.fact, x.b)
end


solve!(mm)
