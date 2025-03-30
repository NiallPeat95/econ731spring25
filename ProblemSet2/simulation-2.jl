#
#   Using hat algebra to solve MSEK model with simulated data.
#
#   Nels Lind, 2/21/2025
#

cd("/Users/niallpeat/Documents/GitHub/econ731spring25/ProblemSet2")
using Pkg
Pkg.activate("."); Pkg.instantiate()
using Pkg; Pkg.add(["FileIO", "DataFrames", "Chain", "Plots", "Distributions", "LinearAlgebra"])
using FileIO, DataFrames, Chain, Plots, Distributions, LinearAlgebra, FileIO

include("MSEK.jl")

# Simulation
J = 56
N = 44

# Generate random matrices for Π and Π_l
Π = [I + rand(Pareto(1, 0.001), N, N) for j = 1:J]
Π = [Π[j] ./ sum(Π[j], dims=1) for j = 1:J]
Π = cat(addDim.(Π, 1)..., dims=1)

Π_l = [I + rand(Pareto(1, 0.001), N, N) for j = 1:J]
Π_l = [Π_l[j] ./ sum(Π_l[j], dims=1) for j = 1:J]
Π_l = cat(addDim.(Π_l, 1)..., dims=1)

# Generate other random values
Y = rand(LogNormal(0, 2.), N)
D = rand(Uniform(-0.02, 0.02), N) .* Y
D = D .- mean(D)
α = rand(Uniform(0, 1), J, N)
α = α ./ sum(α, dims=1)  # alpha is preference over sector output
θ = rand(Uniform(2, 8), J)
θ = repeat(reshape(θ, 1, J), N, 1)

# μ is a J*N matrix of labor supply preferences
μ = rand(Uniform(0.1, 0.5), J, N)  # Random values for μ that are preference for labor supply in sector j and country n

# v is a constant float64 value (0.5)
v = ones(J)*0.5

include("MSEK.jl")

# Create an MSEK instance with the updated struct
m = MSEK(Π, Y, D, α, θ, μ, Π_l, v)

T̂ = ones(J,N)
μ̂ = ones(J, N)

# Set the entry for China-Manufacturing (you can specify the exact row/column indices)
μ̂[10, 23] = 1.1  # or another value depending on your needs

D′ = copy(m.D)
tol=1e-16;maxit=1e4;report=true

sum(Π,dims=2)

include("MSEK.jl")

Ŵ = tâtonnment(m,T̂,D′,report=true)

P̂ = exp.(dsum(μ.*log.(prices(m,Ŵ,T̂,τ̂,t′)),dims=1))

selfShare = [ Π[:,n,n]' for n=1:N ]



# This is for welfare analysis, ignore for now
scatter(log.(Y),Ŵ./P̂,legend=false)
scatter(D./Y,Ŵ./P̂,legend=false)
scatter(selfShare,Ŵ./P̂,legend=false)

