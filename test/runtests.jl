using StochasticDominance
using DataFrames
using Dates
using Test

@testset "Testing StochasticDominance.jl" begin

    # Load data
    data = DataFrame(
        Date = Date.([Date("2024-07-01"), Date("2024-07-02"), Date("2024-07-03"), Date("2024-07-05"), Date("2024-07-08"), Date("2024-07-09"), Date("2024-07-10"), Date("2024-07-11"), Date("2024-07-12"), Date("2024-07-15"), Date("2024-07-16"), Date("2024-07-17"), Date("2024-07-18"), Date("2024-07-19"), Date("2024-07-22"), Date("2024-07-23"), Date("2024-07-24"), Date("2024-07-25"), Date("2024-07-26"), Date("2024-07-29"), Date("2024-07-30"), Date("2024-07-31")]),
        Agric = [-1.01, -2.50, -0.38, -1.11, 0.44, 0.05, -0.34, 4.00, 1.76, 1.53, 2.77, 0.71, -2.24, 0.08, 0.35, 2.55, -2.21, 1.67, 0.17, -0.97, -0.11, 0.24],
        Food = [-0.72, -0.22, 0.27, 0.15, -0.22, -1.49, 0.79, 1.60, 1.04, 0.02, 1.28, 0.69, -1.24, -1.51, 0.29, -0.09, 2.88, -0.21, 1.26, -0.79, 0.79, 1.89],
        Soda = [1.10, -0.86, -0.21, -0.33, -0.64, 1.52, 1.31, 2.53, -0.06, -3.42, 2.59, -1.44, -1.74, -0.23, -0.54, 0.08, -1.61, 1.05, 2.29, -0.48, 0.98, 0.72],
        Beer = [-1.80, -0.09, 0.41, 2.10, 0.59, -2.59, 2.75, 0.73, 0.41, -2.50, 0.40, 1.12, -1.72, -2.67, -0.42, -0.39, -2.51, 0.00, 0.93, -0.48, -1.65, -1.14],
        Smoke = [-0.65, 0.64, 0.41, -1.61, 0.33, -0.26, 1.93, 3.72, 2.17, -0.63, 1.83, 0.94, -1.78, -1.81, -0.25, 1.14, -0.41, 2.21, 2.13, -1.13, 0.07, -1.23]
    )

    #Define Portfolio matrix ξ ∈ R^{d×n} d assets and n scenarios
    ξ =  Matrix(select(data, Not(:Date)))'
    d, n = size(ξ) 
    # Define Benchmark 
    τ = fill(1/d,d)  # equally weights
    ξ_0 = vec(τ'*ξ)

    # Probability vectors for ξ and ξ_0
    p_ξ = fill(1/n,n)  
    p_ξ_0 = fill(1/n,n)  
    
    
    # Stochastic Order
    SDorder = 5.0 

    # parameter(s)

    β = 0.5    

    
    @testset "Mean Return: simplex" begin
        x_opt, t_opt = StochasticDominanceMeanReturn(ξ, ξ_0,p_ξ, p_ξ_0,SDorder;)

        @test isapprox(sum(x_opt),1, atol=1e-4)
    end # end of Mean Return: simplex test

    @testset "Mean Return: objective" begin
        x_opt, t_opt = StochasticDominanceMeanReturn(ξ, ξ_0,p_ξ, p_ξ_0,SDorder;)

        @test MeanReturn(x_opt,ξ,p_ξ) ≥ p_ξ_0'*ξ_0
    end # end of Mean Return: objective test 

    @testset "Mean Return: Higher Order Stochastic Dominance" begin
        x_opt, t_opt = StochasticDominanceMeanReturn(ξ, ξ_0,p_ξ, p_ξ_0,SDorder;)

        @test  isapprox(g_bar(x_opt, ξ, ξ_0, SDorder-1, p_ξ, p_ξ_0),0, atol=1e-4) 
    end # end of Mean Return: Higher Order Stochastic Dominance test 

    @testset "Risk Function: simplex" begin
        x_opt, q_opt, t_opt = StochasticDominanceRiskMeasure(ξ, ξ_0,p_ξ, p_ξ_0,SDorder;β)

        @test isapprox(sum(x_opt),1, atol=1e-4)
    end # end of Risk Function test

    @testset "Risk Function: objective" begin
        x_opt, q_opt, t_opt = StochasticDominanceRiskMeasure(ξ, ξ_0,p_ξ, p_ξ_0,SDorder;β)

        @test RiskFunction(x_opt, q_opt, ξ, SDorder-1, p_ξ, β) ≤  BenchmarkRiskFunction(q_opt, ξ_0, SDorder-1, p_ξ_0, β)
    end # end of Risk Function test
    
    @testset "Mean Return: Higher Order Stochastic Dominance" begin
        x_opt, q_opt, t_opt = StochasticDominanceRiskMeasure(ξ, ξ_0,p_ξ, p_ξ_0,SDorder;β)

        @test  isapprox(g_bar(x_opt, ξ, ξ_0, SDorder-1, p_ξ, p_ξ_0),0, atol=1e-4) 
    end # end of Mean Return: Higher Order Stochastic Dominance test 

end # end of master test

