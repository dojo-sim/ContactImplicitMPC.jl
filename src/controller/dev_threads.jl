using BenchmarkTools
using StaticArrays
Threads.nthreads()

mutable struct LinSys{T,n,nn}
	A::SMatrix{n,n,T,nn}
	b::SVector{n,T}
	x::SVector{n,T}
end

linsys = [LinSys(SMatrix{10,10}(rand(10, 10)), rand(SVector{10}), zeros(SVector{10})) for t = 1:1000]

function solve_linsys!(s::LinSys)
	for t = 1:10000
		s.x = s.A \ s.b
	end
	nothing
end

function solve_all!(s::Vector{LinSys{T,n,nn}}) where {T,n,nn}
	for ss in s
		solve_linsys!(ss)
	end
	nothing
end

function solve_all_threads!(s::Vector{LinSys{T,n,nn}}) where {T,n,nn}
	Threads.@threads for ss in s
		solve_linsys!(ss)
	end
	nothing
end

@code_warntype solve_all!(linsys)
@code_warntype solve_all_threads!(linsys)

@time solve_all!(linsys)
@time solve_all_threads!(linsys)

@benchmark solve_all!($linsys)
@benchmark solve_all_threads!($linsys)

@profiler solve_all!(linsys)
@profiler solve_all_threads!(linsys)
