# Reference trajectory
model = get_model("quadruped", surf = "flat")
@load joinpath(pwd(), "src/dynamics/quadruped/gaits/gait_1.jld2") z̄ x̄ ū q̄ τ̄ λ̄ b̄ h̄

# time
h = mean(h̄)
T = length(τ̄)

[norm(dynamics(model, h̄[t], q̄[t], q̄[t+1], h̄[t] * τ̄[t], zeros(model.dim.w), h̄[t] * λ̄[t], h̄[t] * b̄[t], q̄[t+2]), Inf) for t = 1:T]


# initial conditions
q0 = SVector{model.dim.q}(q̄[1])
q1 = SVector{model.dim.q}(q̄[2])

function z_initialize!(z, model::Quadruped, q1)
	nq = model.dim.q

    z .= 1.0e-1
    z[1:nq] = q1
end

# simulator
sim = ContactControl.simulator(model, q0, q1, h, T,
    u = [SVector{model.dim.u}(h * u) for u in τ̄],
    r! = model.res.r, rz! = model.res.rz, rθ! = model.res.rθ,
    rz = model.spa.rz_sp,
    rθ = model.spa.rθ_sp,
    ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-6, κ_tol = 1.0e-6, κ_init = 1.0e-3),
    sim_opts = ContactControl.SimulatorOptions(warmstart = false))

@time step!(sim, 1)


# simulate
status = ContactControl.simulate!(sim, verbose = true)

include(joinpath(pwd(), "src/dynamics/quadruped/visuals.jl"))
vis = Visualizer()
render(vis)
# visualize!(vis, model, q̄, Δt = h)
visualize!(vis, model, sim.q, Δt = h)

# benchmark
using BenchmarkTools

A = copy(sim.ip.rz)
Ad = Array(copy(A))
b = copy(sim.ip.r)
x = zero(b)
sol = Ad \ b

function solve1!(x, A, b)
	x .= A \ b
end
@benchmark solve1!($x, $A, $b)
@benchmark solve1!($x, $Ad, $b)

function solve2!(x, A, b)
	x .= lu(A) \ b
end
@benchmark solve2!($x, $A, $b)
@benchmark solve2!($x, $Ad, $b)

function solve3!(x, A, b)
	x .= factorize(A) \ b
end
@benchmark solve3!($x, $A, $b)
@benchmark solve3!($x, $Ad, $b)

function solve4!(x, A, b)
	ldiv!(x, lu(A), b)
end
@benchmark solve4!($x, $A, $b)
@benchmark solve4!($x, $Ad, $b)

function solve5!(x, A, b)
	ldiv!(x, factorize(A), b)
end
@benchmark solve5!($x, $A, $b)
@benchmark solve5!($x, $Ad, $b)

function solve6!(x, A, b)
	ldiv!(x, factorize!(A), b)
end
# @benchmark solve6!($x, $A, $b)
@benchmark solve6!($x, $Ad, $b)

function solve7!(x, A, b)
	ldiv!(x, lu!(A), b)
end
# @benchmark solve6!($x, $A, $b)
@benchmark solve7!($x, $copy(Ad), $b)
Ad

function stuff!(A, Ad)
	Ad .= A
end
@benchmark stuff!($A, $Ad)
Ad = Array(copy(A))
info = LAPACK.getrf!(Ad)
Ad = Array(copy(A))
LAPACK.getrf
x .= b
# @benchmark LAPACK.getrs!('N', $info_cache[1], $info_cache[2], $x)
@benchmark LAPACK.getrs!($('N'), $(info_cache[1]), $(info_cache[2]), $x)
lu!(Ad, check = false)
norm(s - sol)
Ad
LAPACK.require_one_based_indexing(Ad)
LAPACK.chkstride1(Ad)
m, n = size(Ad)
lda  = max(1,stride(Ad, 2))
ipiv = similar(Ad, Int, min(m,n))
info = Ref{Int}()
function getrf!(A::AbstractMatrix{LAPACK.elty})
	require_one_based_indexing(A)
	chkstride1(A)
	m, n = size(A)
	lda  = max(1,stride(A, 2))
	ipiv = similar(A, BlasInt, min(m,n))
	info = Ref{BlasInt}()
	ccall((LAPACK.@blasfunc(getrf), liblapack), Cvoid,
		  (Ref{BlasInt}, Ref{BlasInt}, Ptr{LAPACK.elty},
		   Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}),
		  m, n, A, lda, ipiv, info)
	chkargsok(info[])
	A, ipiv, info[] #Error code is stored in LU factorization type
end
