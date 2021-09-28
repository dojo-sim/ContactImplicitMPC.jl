@testset "Static Dense Modified Gram Schmidt" begin
	T = Float64
	n = 16
	A = rand(n,n)
	a = Vector{SVector{n,T}}([SVector{n,T}(A[:,i]) for i=1:n])
	gs_solver = ContactImplicitMPC.SDMGSSolver(n)
	ContactImplicitMPC.factorize!(gs_solver, A)
	ContactImplicitMPC.factorize!(gs_solver, a)
	ContactImplicitMPC.factorize!(gs_solver)

	b = rand(SVector{n,T})
	ContactImplicitMPC.qr_solve!(gs_solver, b)
	@test norm(A*gs_solver.xs - b, Inf) < 1e-10
	# @benchmark A\b

	n = 43
	m = 33
	A = rand(n,n)
	gs_solver = ContactImplicitMPC.SDMGSSolver(A)
	ContactImplicitMPC.factorize!(gs_solver, A)
	X = rand(n,m)
	B = rand(n,m)
	ContactImplicitMPC.qr_matrix_solve!(gs_solver, X, B)
	@test norm(A * X - B) < 1e-10
end
