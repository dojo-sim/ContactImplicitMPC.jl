
function process(x, u; w=zeros(nw))
	q0, q1, γ0, b0 = unpack_x(x)
	sim = simulator(model, SVector{nq}(q0), SVector{nq}(q1), h, 1,
		p=open_loop_policy([u], h),
		d=open_loop_disturbances([w], h),
		ip_opts=InteriorPointOptions(r_tol=1e-8, κ_tol=2e-7, κ_init=1e-3, diff_sol=true),
		sim_opts=SimulatorOptions())
	simulate!(sim)
	q2 = sim.traj.q[3]
	γ1 = sim.traj.γ[1]
	b1 = sim.traj.b[1]
	x̄ = [q1; q2; γ1; b1]

	F = zeros(nx,nx)
	c1 = (1:nq)
	c2 = nq .+ (1:nq)

	off = 0
	F[off .+ (1:nq), c2] += I(nq); off += nq

	F[off .+ (1:nq), c2] = sim.deriv_traj.dq2dq1[1]; off += nq
	F[off .+ (1:nq), c1] = sim.deriv_traj.dq2dq0[1]

	F[off .+ (1:nc), c1] = sim.deriv_traj.dγdq0[1]
	F[off .+ (1:nc), c2] = sim.deriv_traj.dγdq1[1]; off += nc

	F[off .+ (1:nb), c1] = sim.deriv_traj.dbdq0[1]
	F[off .+ (1:nb), c2] = sim.deriv_traj.dbdq1[1]; off += nb

	H = zeros(nq,nx)
	H[1:nq,0  .+ (1:nq)] = copy(sim.deriv_traj.dq2dq0[1])
	H[1:nq,nq .+ (1:nq)] = copy(sim.deriv_traj.dq2dq1[1])
	return x̄, F, H
end




#
# x = rand(nx)
# u = rand(nu)
# w = rand(nw)
#
# process(x, u; w=w)
