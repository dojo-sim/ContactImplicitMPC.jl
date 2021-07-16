################################################################################
# Linear Program
################################################################################
using Random
using LinearAlgebra

"""
https://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf

min   c' * x
s.t.  G * x + s = h
A * x = b
s ≧ 0

Primal variables:
x ∈ [n]
s ∈ [m]

Problem data:
c ∈ [n]
h ∈ [m]
b ∈ [p]
G ∈ [m,n]
A ∈ [p,n]

"""

################################################################################
# Methods
################################################################################

mutable struct LinearProgram11{T}
      n::Int
      m::Int
      p::Int
      c::AbstractVector{T}
      h::AbstractVector{T}
      b::AbstractVector{T}
      G::AbstractMatrix{T}
      A::AbstractMatrix{T}
end

function linear_program(n::Int, m::Int, p::Int; seed::Int = 100)
      Random.seed!(seed)
      c = rand(n)
      h = rand(m)
      b = rand(p)
      G = rand(m, n)
      A = rand(p, n)
      @assert rank(A) == p
      LinearProgram11(n, m, p, c, h, b, G, A)
end

mutable struct Candidate11{T}
      n::Int
      m::Int
      p::Int
      s::AbstractVector{T}
      κ::T
      x::AbstractVector{T}
      y::AbstractVector{T}
      z::AbstractVector{T}
      τ::T
end

function candidate(n::Int, m::Int, p::Int;)
      κ0 = 1.0
      τ0 = 1.0
      s0 = ones(m)
      x0 = ones(n)
      y0 = ones(p)
      z0 = ones(m)
      @assert all(s0 .> 0.0)
      @assert all(z0 .> 0.0)
      Candidate11(n, m, p, s0, κ0, x0, y0, z0, τ0)
end
################################################################################
# Problem Setup
################################################################################

n = 10
m = 5
p = 3

lp = linear_program(n, m, p)
ca = candidate(n, m, p)
