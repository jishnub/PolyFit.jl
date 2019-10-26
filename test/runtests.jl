using PolyFit
using Test,Polynomials,LinearAlgebra

@testset "no errors" begin
    # Should get the same answer as Polynomials.polyfit
    x = 1:10

    for m = 1:5
	    y = x .^ m
	    for n = 1:m
		    P = polyfit(PolyFit.linear_system_type(length(x),n),x,y,n)
		    P2 = Polynomials.polyfit(x,y,n)
		    @test P ≈ P2
		end
	end
end

@testset "same errors" begin
    # Should get the same answer as Polynomials.polyfit
    x = 1:5
	y = x .^ 2
	    
    P = polyfit(x,y,2,1)
    P2 = Polynomials.polyfit(x,y,2)
    @test P ≈ P2

    P = polyfit(x,y,2,I)
	@test P ≈ P2
end

@testset "independent errors" begin
    # Should get the same answer as Polynomials.polyfit
    x = 1:5
	y = x .^ 2
	    
    P = polyfit(x,y,2,ones(length(x)))
    P2 = Polynomials.polyfit(x,y,2)
    @test P ≈ P2

	P = polyfit(x,y,2,Diagonal(ones(length(x))))
	@test P ≈ P2

	P = polyfit(x,y,2,Diagonal(2ones(length(x))))
	@test P ≈ P2
end

@testset "covariance matrix" begin
    x = 1:5
	y = x .^ 2

	# Should get the same answer as Polynomials.polyfit
	P = polyfit(x,y,2,Array(Diagonal(ones(length(x)))))
	P2 = Polynomials.polyfit(x,y,2)
    @test P ≈ P2
end