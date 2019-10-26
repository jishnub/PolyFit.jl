module PolyFit

using SpecialMatrices,LinearAlgebra, Polynomials
import Polynomials: polyfit

export polyfit,
SquareSystem,OverDeterminedSystem,UnderDeterminedSystem

abstract type LinearSystemDimension end
struct SquareSystem <: LinearSystemDimension end
struct OverDeterminedSystem <: LinearSystemDimension end
struct UnderDeterminedSystem <: LinearSystemDimension end

function linear_system_type(M::Integer,N::Integer)
	# M is the number of equations, N is the number of parameters
	if M == N
		return SquareSystem()
	elseif M < N
		return OverDeterminedSystem()
	else
		return UnderDeterminedSystem()
	end
end

@inline function _polyfit(V::AbstractMatrix,y::AbstractVector{<:Real},
	::Union{UniformScaling,Real,Nothing})

	a = V \ y
end

@inline function _polyfit(V::AbstractMatrix,y::AbstractVector{<:Real},
	Σ::AbstractMatrix)

	V′ =  Σ * V
	y′ =  Σ * y
	a = V′ \ y′
end

function polyfit(::UnderDeterminedSystem,x::AbstractVector{<:Real},
	y::AbstractVector{<:Real},n::Integer,sym::Symbol=:x,scaling=nothing)

	if isnothing(scaling)
		return polyfit(x,y,n,sym)
	end

	# Returns the least-square solution
	V = Vandermonde(collect(x))[:,1:n+1]
	a = _polyfit(V,y,scaling)
	Poly(a,sym)
end

function polyfit(::SquareSystem,x::AbstractVector{<:Real},
	y::AbstractVector{<:Real},n::Integer,sym::Symbol=:x,scaling=nothing)

	V = Vandermonde(collect(x))
	a = _polyfit(V,y,scaling)
	Poly(a,sym)
end

function polyfit(::OverDeterminedSystem,x::AbstractVector{<:Real},
	y::AbstractVector{<:Real},n::Integer,sym::Symbol=:x,scaling=nothing)

	V = Vandermonde(collect(x))
	a = _polyfit(V,y,scaling)
	a = vcat(a,zeros(n - length(x)))
	Poly(a,sym)
end

function polyfit(x::AbstractVector{<:Real},y::AbstractVector{<:Real},
	n::Integer,σy::Union{Real,UniformScaling},sym::Symbol=:x)
	
	T = linear_system_type(length(x),n)
	a = polyfit(T,x,y,n,sym,σy)
end

function polyfit(x::AbstractVector{<:Real},y::AbstractVector{<:Real},
	n::Integer,σy::AbstractVector{<:Real},sym::Symbol=:x)

	T = linear_system_type(length(x),n)
	sqrtΣinv = Diagonal(@. 1/σy)
	polyfit(T,x,y,n,sym,sqrtΣinv)
end

function polyfit(x::AbstractVector{<:Real},y::AbstractVector{<:Real},
	n::Integer,Σ::Diagonal{<:Real},sym::Symbol=:x)

	T = linear_system_type(length(x),n)
	sqrtΣinv = sqrt(Σ^-1)
	polyfit(T,x,y,n,sym,sqrtΣinv)
end

function polyfit(x::AbstractVector{<:Real},y::AbstractVector{<:Real},
	n::Integer,C::AbstractMatrix{<:Real},sym::Symbol=:x)

	T = linear_system_type(length(x),n)
	L = cholesky(C).L
	Linv = L^-1
	polyfit(T,x,y,n,sym,Linv)
end

end # module
