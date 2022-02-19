n = 3
m = 2
T = 5

nz = n * (T - 1) + m * (T - 1)
nd = n * (T - 1)

u_idx = [collect((t - 1) * (m + n) .+ (1:m)) for t = 1:T-1]
x_idx = [collect((t - 1) * (m + n) + m .+ (1:n)) for t = 1:T-1]
n_idx = [(t - 1) * n .+ (1:n) for t = 1:T-1]

A = [rand(n,n) for t = 1:T-1]
B = [rand(n,m) for t = 1:T-1]
Q = [Diagonal(rand(n)) for t = 1:T]
R = [Diagonal(rand(m)) for t = 1:T-1]

x_ref = [randn(n) for t = 1:T-1]
u_ref = [randn(m) for t = 1:T-1]

x_init = ones(n)
z0 = rand(nz + nd)

function constraints(z)
	c = zeros(eltype(z), nd)

	for t = 1:T-1
		x2 = z[x_idx[t]]
		x1 = (t == 1 ? x_init : z[x_idx[t-1]])
		u1 = z[u_idx[t]]
		c[n_idx[t]] = x2 - A[t] * x1 - B[t] * u1
	end

	return c
end

constraints(z0)

function constraints_jacobian(z)
	C = zeros(eltype(z), nd, nz)

	for t = 1:T-1
		C[n_idx[t], u_idx[t]] = -B[t]
		C[n_idx[t], x_idx[t]] = Diagonal(ones(n))
		t == 1 && continue
		C[n_idx[t], x_idx[t-1]] = -A[t]
	end

	return C
end

constraints_jacobian(z0)

function objective_gradient(z)
	j = zeros(eltype(z), nz)
	for t = 1:T-1
		u1 = z[u_idx[t]]
		x1 = z[x_idx[t]]
		j[u_idx[t]] = R[t] * (u1 - u_ref[t])
		j[x_idx[t]] = Q[t] * (x1 - x_ref[t])
	end
	return j
end

objective_gradient(z0)

function objective_hessian(z)
	J = zeros(eltype(z), nz, nz)
	for t = 1:T-1
		J[u_idx[t], u_idx[t]] = R[t]
		J[x_idx[t], x_idx[t]] = Q[t]
	end
	return J
end

objective_hessian(z0)

function hessian(z)
	J = objective_hessian(z)
	C = constraints_jacobian(z)

	return [J C'; C zeros(eltype(z), nd, nd)]
end

function gradient(z)
	y = z[nz .+ (1:nd)]
	return [objective_gradient(z) + constraints_jacobian(z)' * y; constraints(z)]
end

u = [randn(m) for t = 1:T-1]
x = [x_init]
for t = 1:T-1
	push!(x, A[t] * x[end], + B[t] * u[t])
end
z = [vcat([[u[t]; x[t]] for t = 1:T-1]...); zeros(nd)]

res_norm = norm(gradient(z), 1)

Δ = hessian(z) \ gradient(z)

α = 1.0
z_cand = z - α * Δ

res_cand_norm = norm(gradient(z_cand), 1)

z_cand[nz .+ (1:nd)]
