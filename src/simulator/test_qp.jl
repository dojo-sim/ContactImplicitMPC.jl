include("interior_point.jl")
using InteractiveUtils, BenchmarkTools, ModelingToolkit

"""
    minimize   x' P x + q' x
    subject to    x >= 0
"""

n = 1000

_P = @SVector ones(n)
P = Diagonal(_P)
q = @SVector ones(n)
θ = [_P; q]
z = ones(2 * n)
r = zeros(2 * n)
rz = zeros(2 * n, 2 * n)
rθ = zeros(2 * n, 2 * n)

# residual
function _r!(r, z, θ)
    x = z[1:1000]
    y = z[1:1000]
    P = θ[1001:2000]
    q = θ[1001:2000]
    κ = θ[end]

    r[1:1000] = 2.0 * Diagonal(P) * x + q - y
    r[1001:2000] = Diagonal(x) * y .- κ
    nothing
end

@code_warntype _r!(r, z, θ)
@benchmark _r!($r, $z, $θ)

@variables r_sym[1:2000]
@variables z_sym[1:2000]
@variables θ_sym[1:2000]

parallel = false# ModelingToolkit.MultithreadedForm()
_r!(r_sym, z_sym, θ_sym)
r_sym = simplify.(r_sym)
r! = eval(ModelingToolkit.build_function(r_sym, z_sym, θ_sym,
    parallel = parallel)[2])
rz_exp = ModelingToolkit.sparsejacobian(r_sym, z_sym, simplify = true)
rθ_exp = ModelingToolkit.sparsejacobian(r_sym, θ_sym, simplify = true)
rz_sp = similar(rz_exp, Float64)
rθ_sp = similar(rθ_exp, Float64)
rz! = eval(ModelingToolkit.build_function(rz_exp, z_sym, θ_sym,
    parallel = parallel)[2])
rθ! = eval(ModelingToolkit.build_function(rθ_exp, z_sym, θ_sym,
    parallel = parallel)[2])

function r!(r, z, θ::ResidualData)
    r!(z, z, θ.θ)
end
@code_warntype r!(r, z, θ)
@benchmark r!($r, $z, $θ)

@code_warntype rz!(rz_sp, z, θ)
rz!(rz_sp, z, θ)
@benchmark rz!($rz_sp, $z, $θ)

@code_warntype rθ!(r, z, θ)
rθ!(rθ_sp, z, θ)
@benchmark rθ!($rθ_sp, $z, $θ)

num_var = 2 * n
num_data = 2 * n
idx_ineq = collect(1:num_var)
ip_data = interior_point_data(num_var, num_data, idx_ineq;
    rz = rz_sp,
    rθ = rθ_sp)

interior_point!(ip_data)#, opts = InteriorPointOptions())
