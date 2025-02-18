#
#   Using hat algebra to solve DEK with simulated data
#
#   Nels Lind, 2/17/2025
#


username = "nelslind"
cd("/Users/$username/Dropbox/teaching/emory/2024-2025/Econ 731 Spring 2025/code/econ731spring25/")
using Pkg
Pkg.activate("."); Pkg.instantiate()
using FileIO, DataFrames, Chain, Plots
using Distributions, LinearAlgebra

include("DEK.jl")

# example
N = 100
Y = rand(LogNormal(0,2.),N)
D = rand(Uniform(-.1,.1),N) .* Y
D = D .- mean(D)
Dm = rand(Uniform(-.1,.1),N) .* Y
Dm = Dm .- mean(Dm)
Π = I + .1*rand(N,N)
Π = Π ./ sum(Π,dims=1)
m = DEK(Π,Y,D,Dm,.5,.5,1.)
D′ = zeros(N)
Dm′ = zeros(N)
@time Ŵ = tâtonnment(m,D′,Dm′,report=true,reportrate=1)

scatter(D,Ŵ.*Y,legend=false)
scatter(Dm,Ŵ.*Y,legend=false)

scatter(log.(Y),log.(Ŵ),legend=false)