include(joinpath(pwd(), "src/trajectory_optimization/moi.jl"))
include(joinpath(pwd(), "src/trajectory_optimization/utils.jl"))

T = 10
h = 0.1
n = 2
m = 1

x1 = [0.0; 0.0]
xT = [π; 0.0]

nz = n * T + m * (T - 1)
np = n * (T - 1)

z0 = rand(nz)

function f(x1, u1, x2)

	function dyn(x, u)
		[x[2];
        (u[1] / 0.25
        - 9.81 * sin(x[1]) / 0.5
        - 0.1 * x[2] / 0.25)]
	end

	return x2 - (x1 + h * dyn(0.5 * (x1 + x2), u1))
end

function obj(xt, ut)
	Qt = Diagonal(ones(n))
	Rt = Diagonal(ones(m))
	return transpose(xt - xT) * Qt * (xt - xT) + transpose(ut) * Rt * ut
end

function obj(x)
	QT = Diagonal(ones(n))
	return transpose(x - xT) * QT * (x - xT)
end

@variables x1_sym[1:n]
@variables x2_sym[1:n]
@variables u1_sym[1:m]

f_sym = f(x1_sym, u1_sym, x2_sym)
f_sym = simplify.(f_sym)
fx1_sym = Symbolics.jacobian(f_sym, x1_sym, simplify = true)
fu1_sym = Symbolics.jacobian(f_sym, u1_sym, simplify = true)
fx2_sym = Symbolics.jacobian(f_sym, x2_sym, simplify = true)

f_func = eval(Symbolics.build_function(f_sym, x1_sym, u1_sym, x2_sym)[1])
fx1_func = eval(Symbolics.build_function(fx1_sym, x1_sym, u1_sym, x2_sym)[1])
fu1_func = eval(Symbolics.build_function(fu1_sym, x1_sym, u1_sym, x2_sym)[1])
fx2_func = eval(Symbolics.build_function(fx2_sym, x1_sym, u1_sym, x2_sym)[1])

f_func(rand(n), rand(n), rand(m))
fx1_func(rand(n), rand(n), rand(m))
fu1_func(rand(n), rand(n), rand(m))
fx2_func(rand(n), rand(n), rand(m))

obj_sym = obj(x1_sym, u1_sym)
obj_sym = simplify.(obj_sym)

objx1_sym = Symbolics.gradient(obj_sym, x1_sym, simplify = true)
obju1_sym = Symbolics.gradient(obj_sym, u1_sym, simplify = true)

obj_func = eval(Symbolics.build_function(obj_sym, x1_sym, u1_sym))
objx1_func = eval(Symbolics.build_function(objx1_sym, x1_sym, u1_sym)[1])
obju1_func = eval(Symbolics.build_function(obju1_sym, x1_sym, u1_sym)[1])

objT_sym = obj(x1_sym)
objT_sym = simplify.(objT_sym)

objTx_sym = Symbolics.gradient(objT_sym, x1_sym, simplify = true)

objT_func = eval(Symbolics.build_function(objT_sym, x1_sym))
objTx1_func = eval(Symbolics.build_function(objTx_sym, x1_sym)[1])

x_idx = [(t - 1) * (n + m) .+ (1:n) for t = 1:T]
u_idx = [(t - 1) * (n + m) + n .+ (1:m) for t = 1:T-1]

function unpack(z)
	x = [z[x_idx[t]] for t = 1:T]
	u = [z[u_idx[t]] for t = 1:T-1]

	return x, u
end

unpack(z0)

function moi_obj(z)
	J = 0.0

	for t = 1:T-1
		xt = view(z, x_idx[t])
		ut = view(z, u_idx[t])
		J += obj_func(xt, ut)
	end

	xT = view(z, x_idx[T])
	J += objT_func(xT)

	return J
end

moi_obj(z0)

function ∇moi_obj!(g, z)
	for t = 1:T-1
		xt = view(z, x_idx[t])
		ut = view(z, u_idx[t])
		g[x_idx[t]] .= objx1_func(xt, ut)
		g[u_idx[t]] .= obju1_func(xt, ut)
	end

	xT = view(z, x_idx[T])
	g[x_idx[T]] .= objT_func(xT)

	nothing
end

g0 = zeros(nz)
∇moi_obj!(g0, z0)

function moi_con!(c, z)
	for t = 1:T-1
		x1 = view(z, x_idx[t])
		u1 = view(z, u_idx[t])
		x2 = view(z, x_idx[t+1])
		c[(t - 1) * n .+ (1:n)] = f_func(x1, u1, x2)
	end
	nothing
end

c0 = zeros(np)
moi_con!(c0, z0)

function ∇moi_con!(j, z)
	shift = 0
	for t = 1:T-1
		x1 = view(z, x_idx[t])
		u1 = view(z, u_idx[t])
		x2 = view(z, x_idx[t+1])

		r_idx = (t - 1) * n .+ (1:n)

		c_idx = x_idx[t]
		s = length(r_idx) * length(c_idx)
		j[shift .+ (1:s)] .= vec(fx1_func(x1, u1, x2))
		shift += s

		c_idx = u_idx[t]
		s = length(r_idx) * length(c_idx)
		j[shift .+ (1:s)] .= vec(fu1_func(x1, u1, x2))
		shift += s

		c_idx = x_idx[t+1]
		s = length(r_idx) * length(c_idx)
		j[shift .+ (1:s)] .= vec(fx2_func(x1, u1, x2))
		shift += s
	end
	nothing
end

j0 = zeros((n * n + n * m + n * n) * (T - 1))
∇moi_con!(j0, z0)

# user can overwrite sparsity_jacobian and sparsity_hessian
function sparsity_jacobian(nz, np)

    row = []
    col = []

	for t = 1:T-1
		r_idx = (t - 1) * n .+ (1:n)

		c_idx = x_idx[t]
		row_col!(row,col,r_idx,c_idx)

		c_idx = u_idx[t]
		row_col!(row,col,r_idx,c_idx)

		c_idx = x_idx[t+1]
		row_col!(row,col,r_idx,c_idx)
	end

    return collect(zip(row, col))
end

sparsity_jacobian(nz, np)

xl, xu = primal_bounds(nz)
xl[1:n] = x1
xu[1:n] = x1

xl[end-n+1:end] = xT
xu[end-n+1:end] = xT

x0 = linear_interpolation(x1, xT, T)
u0 = [randn(m) for t = 1:T-1]
z0 = vcat([[x0[t]; u0[t]] for t = 1:T-1]..., x0[T]...)

prob = ProblemMOI(nz, np,
    sparsity_jac=sparsity_jacobian(nz, np),
    primal_bounds=(xl, xu),
    constraint_bounds=constraint_bounds(np),
    hessian_lagrangian=false)

z_sol = solve(z0, prob,
        nlp=:ipopt,
        tol=1.0e-2,
        c_tol=1.0e-2,
        max_iter=1000)

c0 = zeros(np)
moi_con!(c0, z_sol)
norm(c0)
