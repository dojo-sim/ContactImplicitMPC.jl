using BenchmarkTools
Threads.nthreads()

struct LinSys{T}
	A::Array{T, 2}
	b::Vector{T}
	x::Vector{T}
end

linsys = [LinSys(rand(1000, 1000), rand(1000), zeros(1000)) for t = 1:1000]

function solve_linsys!(s::LinSys)
	s.x .= s.A \ s.b
	nothing
end

function solve_all!(s::Vector{LinSys{T}}) where T
	for ss in s
		solve_linsys!(ss)
	end
	nothing
end

function solve_all_threads!(s::Vector{LinSys{T}}) where T
	Threads.@threads for ss in s
		solve_linsys!(ss)
	end
	nothing
end

@code_warntype solve_all!(linsys)
@code_warntype solve_all_threads!(linsys)

@benchmark solve_all!($linsys)
@benchmark solve_all_threads!($linsys)

@profile solve_all_threads!(linsys)
