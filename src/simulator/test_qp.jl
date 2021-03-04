include("interior_point.jl")
using InteractiveUtils, BenchmarkTools, ModelingToolkit

"""
    minimize   x' P x + q' x
    subject to    x >= 0
"""

n = 100

_P = @SVector ones(n)
P = Diagonal(_P)
q = @SVector ones(n)
θ = [_P; q]
z = ones(2 * n)
r = zeros(2 * n)

# residual
function r!(r, z, θ)
    x = view(z, 1:100)
    y = view(z, 1:100)
    P = view(θ, 101:200)
    q = view(θ, 101:200)
    κ = θ[end]

    r[1:100] = P .* x
    r[1:100] .*= 2.0
    r[1:100] += q
    r[1:100] -= y

    r[101:200] = x .* y
    r[101:200] .-= κ
    nothing
end

@code_warntype r!(r, z, θ)
@benchmark r!($r, $z, $θ)

@variables r_sym[1:200] z_sym[1:200] θ_sym[1:200]

r!(r_sym, z_sym, θ_sym)
r_sym = simplify.(r_sym)
r_fast! = eval(ModelingToolkit.build_function(r_sym, z_sym, θ_sym)[2])
ModelingToolkit.jacobian(r_sym, )
ModelingToolkit.build_function(ModelingToolkit.jacobian(r_sym, parallel=false)
ModelingToolkit.jacobian(r_sym, θ_sym)


rz_fast! = eval(ModelingToolkit.build_function(ModelingToolkit.jacobian(r_sym, z_sym, simplify=true))[2])
rθ = ModelingToolkit.jacobian(r, θ, simplify=true)


# Build function
expr = Dict{Symbol, Expr}()
expr[:r]  = build_function(r,  dt, z, θ, ρ)[2]
expr[:rz] = build_function(rz, dt, z, θ, ρ)[2]
expr[:rθ] = build_function(rθ, dt, z, θ, ρ)[2]
@code_warntype r_fast!(r, z, θ)
@benchmark r_fast!($r, $z, $θ)

r
# residual Jacobian wrt z
function rz!(rz, z, θ)
    @warn "residual Jacobian wrt z not defined"
    nothing
end

# # residual Jacobian wrt θ
# function rθ!(rθ, z, θ)
#     @warn "residual Jacobian wrt θ not defined"
#     nothing
# end
