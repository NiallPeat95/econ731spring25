#
#   Defines the MSEK type for representing a multisector EK model as in 
#       Caliendo Parro (2015).
#
#   Nels Lind, 2/20/2025
#
#   j = 1,…,J, sectors
#
#   P_jn = ( ∑_o T_jo ( κ_jon C_jo )^-θ_j )^(-1/θ_j)
#
#   where κ_jon combines physical trade costs and tariffs in sector j
#
#   κ_jod = τ_jod * (1 + t_jod)
#
#   Within sector j, fraction of goods sourced from o by d
#
#   π_jod = T_jo ( κ_jod C_jo )^-θ_j / P_jn )^(-1/θ_j)
#
#   Cobb-Douglas preferences for sectoral composites with shares μ_jn
#
#   P_n = ∏_j P_jn^μ_jn
#
#   π_od = ∑_j π_jod * μ_jn
#
#   Cobb-Douglas production with γ_jkn = cost share of sector k in sector j
#       within country n. Share of value added = 1 - ∑_k γ_jkn
#
#   C_jn = W_n^(1-1 - ∑_i γ_ijn) * ∏_i P_in^γ_ijn
#
#   Imports and exports (excluding tariff payments)
#
#   X_jod = π_jod * X_jd / (1 + t_jod)
#
#   Expenditure on sector i by n
#
#   X_in = ∑_j γ_ijn ∑_d X_jnd + μ_in * I_n
#   I_n = W_n*L_n + R_n + D_n
#
#   where tariff revenue and the trade deficit are
#
#   R_n = ∑_jo t_jon X_jon
#   D_n = ∑_jo X_jon - ∑_jd X_jnd
#
#   Hat Algebra:
#
#   For any given change in fundmanetals (T̂_n,τ̂_od) and counterfatual tariffs t′_jod
#       and deficits D′_n
#
#   κ̂_jod = τ̂_jod + (1+t′_jod)/(1+t_jod)
#   Ĉ_jn = Ŵ_n^(1-1 - ∑_i γ_ijn) * ∏_i P̂_in^γ_ijn
#   P̂_jn = ( ∑_o π_jon T̂_jo ( κ̂_jon Ĉ_jo )^-θ_j )^(-1/θ_j)
#   π̂_jod = T̂_jo ( κ̂_jod Ĉ_jo )^-θ_j / P̂_jn )^(-1/θ_j)
#   X′_jod = π′_jod * X′_jd / (1 + t_jod)
#   X′_in = ∑_j γ_ijn ∑_d X′_jnd + μ_in * I′_n
#   I′_n = Ŵ_n*W_n*L_n + R′_n + D′_n
#   R′n = ∑_jo t′_jon X′_jon
#   ∑_jo X′_jon = ∑_jd X′_jnd + D′_n
#
#   The data requirements are then:
#
#   π_jod = sectoral trade shares
#   Y_n = W_n*L_n = value added
#   t_jod = initial tariffs
#   D_n = initial trade deficits
#   γ_jkn = intermediate input shares
#   μ_jn = final expenditure shares
#   θ_j = sector j trade elasticities

include("utilities.jl")

struct MSEK{T}
    t::Array{T,3}
    Π::Array{T,3}
    Y::Vector{T}
    D::Vector{T}
    γ::Array{T,3}
    μ::Matrix{T}
    θ::Vector{T}
end

# for testing

# simulate data
J = 3
N = 10
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
θ = rand(Uniform(2,8),J)

m = MSEK(t,Π,Y,D,γ,μ,θ)

Ŵ = rand(Uniform(.99,1.01),N)
T̂ = rand(Uniform(.99,1.01),J,N)
τ̂ = rand(Uniform(.99,1.01),J,N,N)
for n = 1:N
    τ̂[:,n,n] .= 1.
end
t′ = m.t .+ .001
D′ = copy(m.D)
tol=1e-16;maxit=1e4;report=true

sum(Π,dims=2)

function prices(m::MSEK{T},Ŵ::Vector{T},T̂::Matrix{T},τ̂::Array{T,3},t′::Array{T,3};tol=1e-16,maxit=1e4,report=false) where {T <:Number}
    #   κ̂_jod = τ̂_jod + (1+t′_jod)/(1+t_jod)
    #   Ĉ_jn = Ŵ_n^(1-1 - ∑_k γ_jkn) * ∏_k P̂_kn^γ_jkn
    #   P̂_jn = ( ∑_o π_jon T̂_jo ( κ̂_jon Ĉ_jo )^-θ_j )^(-1/θ_j)
    P̂ = ones(size(m.μ)...)
    κ̂ = τ̂ .* (1 .+ t′)./(1 .+ m.t)
    done = false
    iter = 0
    while !done
        iter += 1
        P̂old = copy(P̂)
        Ĉ = Ŵ'.^(1 .- dsum(m.γ,dims=1)) .* exp.( dsum( m.γ .* addDim(log.(P̂old),2) ,dims=1) )
        P̂ = dsum( m.Π .* T̂ .* ( κ̂ .* Ĉ ).^(.-m.θ),dims=2).^(.- 1 ./ m.θ)
        err = maximum(abs.(P̂ .- P̂old))
        done = (err < tol) || (iter ≥ maxit)
        if report
            @show iter,err
            # @show P̂
        end
        (iter ≥ maxit) && @warn "Maximum iterations reached when solving for prices."
    end
    return P̂
end
function tradeShares(m::MSEK{T},P̂::Matrix{T},Ŵ::Vector{T},T̂::Matrix{T},τ̂::Array{T,3},t′::Array{T,3}) where {T <:Number}
    κ̂ = τ̂ .* (1 .+ t′)./(1 .+ m.t)
    Ĉ = Ŵ'.^(1 .- dsum(m.γ,dims=1)) .* exp.( dsum( m.γ .* addDim(log.(P̂),2) ,dims=1) )
    out = m.Π .* T̂ .* ( κ̂ .* Ĉ ./ addDim(P̂,2) ).^(.-m.θ)
    return out # ./ sum(out,dims=2)
end
function excessDemand(m::MSEK{T},Ŵ::Vector{T},T̂::Matrix{T},τ̂::Array{T,3},
                        D′::Vector{T},Dm′::Vector{T}) where {T<:Number}
    P̂ = prices(m,Ŵ,T̂,τ̂,t′)
    Π′ = tradeShares(m,P̂,Ŵ,T̂,τ̂,t′)

    Π̃′ = Π′ ./ (1 .+ t′) 


#   X′_in = ∑_j γ_ijn ∑_d π̃′_jnd * X′_jd  + μ_in * ( Ŵ_n*W_n*L_n + ∑_jo t′_jon π̃′_jon * X′_jn  + D′_n )

#   X′_in - ∑_j γ_ijn ∑_d π̃′_jnd * X′_jd - μ_in * ∑_jo t′_jon π̃′_jon * X′_jn 
#        =  μ_in * ( Ŵ_n*W_n*L_n + D′_n )

    Δ = vec(μ .* ( Ŵ.*m.Y + D′ )')

#   X′_jod = π̃′_jod * X′_jd 
#   excessDemand = ∑_jd X′_jnd + D′_n - ∑_jo X′_jon


    return 
end

# excessDemand(m,Ŵ,T̂,τ̂,zeros(N),zeros(N))


function tâtonnment(m::DEK{T},T̂::Vector{T},τ̂::Matrix{T},D′::Vector{T},Dm′::Vector{T};λ = T(.1),tol=1e-10,maxit=1e4,report=false,reportrate=1) where {T<:Number}
    Ŵ = ones(length(m.Y))
    done = false
    iter = 0
    t0 = time()
    while !done
        iter += 1
        Ŵold = copy(Ŵ)
        ed = excessDemand(m,Ŵ,T̂,τ̂,D′,Dm′)
        Ŵ = Ŵold .* (1 .+ λ .* ed ./ (Ŵold .* m.Y))
        err = maximum(abs.(Ŵ .- Ŵold))
        done = (err < tol) || (iter ≥ maxit)
        if report && (time()- t0) > reportrate
            @show iter,err
            # @show Ŵ
            t0 = time()
        end
        (iter ≥ maxit) && @warn "Maximum iterations reached when solving for wages."
    end
    return Ŵ
end
