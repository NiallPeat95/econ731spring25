#
#   Using hat algebra to solve MSEK model with simulated data.
#
#   Nels Lind, 2/21/2025
#

cd("/Users/niallpeat/Documents/GitHub/econ731spring25/ProblemSet2")
#using Pkg
#Pkg.activate("."); Pkg.instantiate()
using Pkg; Pkg.add(["FileIO", "DataFrames", "Chain", "Plots", "Distributions", "LinearAlgebra"])
using FileIO, DataFrames, Chain, Plots, Distributions, LinearAlgebra, FileIO

include("MSEK.jl")

# simulation
J = 56
N = 44
t = rand(Uniform(.01,.2),J,N,N)
Π = [I + rand(Pareto(1,.001),N,N) for j=1:J]
Π = [Π[j] ./ sum(Π[j],dims=1) for j=1:J]
Π = cat(addDim.(Π,1)...,dims=1)
Y = rand(LogNormal(0,2.),N)
D = rand(Uniform(-.02,.02),N) .* Y
D = D .- mean(D)
γ = [ I + rand(Pareto(1,.01),J,J) for n=1:N ]
γ = cat([ rand(Uniform(.2,.3),1,J) .* γ[n] ./ sum(γ[n],dims=1) for n=1:N ]...,dims=3)
μ = rand(Uniform(0,1),J,N)
μ = μ ./ sum(μ,dims=1)
#μ' = μ .* 1.01
θ = rand(Uniform(2,8),J)

m = MSEK(t,Π,Y,D,γ,μ,θ)

T̂ = ones(J,N)
τ̂ = ones(J,N,N)
t′ = m.t 
D′ = copy(m.D)
tol=1e-16;maxit=1e4;report=true

sum(Π,dims=2)

Ŵ = tâtonnment(m,T̂,τ̂,t′,D′,report=true)
P̂ = exp.(dsum(μ.*log.(prices(m,Ŵ,T̂,τ̂,t′)),dims=1))
selfShare = [ Π[:,n,n]' * μ[:,n] for n=1:N ]
scatter(log.(Y),Ŵ./P̂,legend=false)
scatter(D./Y,Ŵ./P̂,legend=false)
scatter(selfShare,Ŵ./P̂,legend=false)

