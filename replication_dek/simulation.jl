#
#   Using hat algebra to solve DEK with simulated data
#
#   Nels Lind, 2/17/2025
#


cd("/Users/niallpeat/Documents/GitHub/econ731spring25/replication_dek/")
#using Pkg; Pkg.add(["FileIO", "DataFrames", "Chain", "Plots", "Distributions", "LinearAlgebra"])
#using FileIO, DataFrames, Chain, Plots
#using Distributions, LinearAlgebra

using Pkg; Pkg.activate("."); Pkg.instantiate()

include("DEK.jl")

# simulate data
N = 20
Y = rand(LogNormal(0,2.),N)
D = rand(Uniform(-.02,.02),N) .* Y
D = D .- mean(D)
X = Y + D
Π = I + rand(Pareto(1,.001),N,N)
Π = Π ./ sum(Π,dims=1)
extrema(diag(Π))

# values in DEK
α = .188
β = .312
θ = 3.6

# implied manufacturing deficits
Dm = (I - (1-β)*Π)\( α * (X - Π*X))
m = DEK(Π,Y,D,Dm,α,β,θ)
extrema(excessDemand(m,ones(N),ones(N),ones(N,N),D,Dm))

# DEK counterfactual of zeroing out current accounts assuming CA = - D
Dm′ = Dm - D
D′ = zeros(N)
T̂ = ones(N)
τ̂ = ones(N,N)
Ŵ = tâtonnment(m,T̂,τ̂,D′,Dm′, maxit = 100000,report=true,reportrate=1,λ=.001)
P̂ = prices(m,Ŵ,T̂,τ̂,maxit = 100000)
scatter(diag(Π),Ŵ./P̂,legend=false)

# reduce trade costs by 10%
D′ = copy(D)
Dm′ = copy(Dm)
T̂ = ones(N)
τ̂ = I + (ones(N,N) - I).*fill(.9,N,N)
Ŵ = tâtonnment(m,T̂,τ̂,D′,Dm′,report=true,reportrate=1,λ=.001)
P̂ = prices(m,Ŵ,T̂,τ̂)
scatter(diag(Π),Ŵ./P̂,legend=false)




