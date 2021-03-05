include("interior_point.jl")
using InteractiveUtils, BenchmarkTools, ModelingToolkit

"""
    minimize   x' P x + q' x
    subject to    x >= 0
"""

n = 1000

_P = ones(n)
P = Diagonal(_P)
q = ones(n)
θ = [_P; q]
z = ones(2 * n)
r = zeros(2 * n)
rz = zeros(2 * n, 2 * n)
rθ = zeros(2 * n, 2 * n)
κ = 1.0

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

# @code_warntype _r!(r, z, θ, κ)
# @benchmark _r!($r, $z, $θ, $κ)

@variables r_sym[1:2000]
@variables z_sym[1:2000]
@variables θ_sym[1:2000]
@variables κ_sym[1:1]

parallel = false # ModelingToolkit.MultithreadedForm()
_r!(r_sym, z_sym, θ_sym, κ_sym)
r_sym = simplify.(r_sym)
r! = eval(ModelingToolkit.build_function(r_sym, z_sym, θ_sym, κ_sym,
    parallel = parallel)[2])
rz_exp = ModelingToolkit.sparsejacobian(r_sym, z_sym, simplify = true)
rθ_exp = ModelingToolkit.sparsejacobian(r_sym, θ_sym, simplify = true)
rz_sp = similar(rz_exp, Float64)
rθ_sp = similar(rθ_exp, Float64)
rz! = eval(ModelingToolkit.build_function(rz_exp, z_sym, θ_sym, κ_sym,
    parallel = parallel)[2])
rθ! = eval(ModelingToolkit.build_function(rθ_exp, z_sym, θ_sym, κ_sym,
    parallel = parallel)[2])

# @code_warntype r!(r, z, θ, κ)
# @benchmark r!($r, $z, $θ, $κ)
#
# @code_warntype rz!(rz_sp, z, θ, κ)
# rz!(rz_sp, z, θ, κ)
# @benchmark rz!($rz_sp, $z, $θ, $κ)
#
# @code_warntype rθ!(r, z, θ, κ)
# rθ!(rθ_sp, z, θ, κ)
# @benchmark rθ!($rθ_sp, $z, $θ, $κ)

num_var() = 2 * n
num_data() = 2 * n
idx_ineq = collect(1:num_var())
ip = interior_point(num_var(), num_data(), idx_ineq,
    rz = rz_sp,
    rθ = rθ_sp)

@time interior_point!(ip, z, θ, opts = InteriorPointOptions(diff_sol = true))
