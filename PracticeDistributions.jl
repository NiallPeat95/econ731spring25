Pkg.add("Distributions")

using Distributions, Plots

α = 3 
θ = 2 
s = 1.0  
T = -s^(α)
F = Frechet(α, s)

z = 0.1:0.1:10  

y = exp.(T.*z.^(-θ))
# Compute CDF
cdf_values = cdf.(F, z)

# Plot
plot(z, cdf_values, label="Fréchet CDF", lw=2, xlabel="z", ylabel="CDF", title="Fréchet Distribution CDF")
plot!(z, y, label="CDF from PDF", lw=2)
