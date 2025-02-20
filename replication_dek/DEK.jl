#
#   Defines the DEK type for representing the Dekle, Eaton, and Kortum (2007) model
#
#   Nels Lind, 2/17/2025
#
#   There are manufactured goods, and final goods. Value added share in manufacturing is β,
#       manufactured good share in final goods is α. θ is the CES trade elasticity. Only
#       manufactured goods can be traded.
#
#   W,L = price and quantitiy of primary factor (labor)
#   Y = income = W*L
#   X = final good expenditure
#   Ym = gross output of manufacturing
#   Xm = expenditure on manufactured goods
#
#   Price level of manufactured goods in d:
#   P_d = sum_o( T_o * (τ_od * W_o^β * P_o^(1-β))^-θ )^(-1/θ)                                   (1)  
#
#   Share of d expenditure sourced from o amongst manufactures:
#   π_od = T_o * ( τ_od * W_o^β * P_o^(1-β) / P_d )^-θ                                          (2)
#
#   For each n = 1,...,N we have:
#
#   X_n = Y_n + D_n                 (definition of trade deficit)                               (3)
#   Xm_n = Ym_n + Dm_n              (definition of manufacturing trade deficit)                 (4)
#   Xm_n = α*X_n + (1-β)*Ym_n       (accounting identity for expenditure on manufacturing)      (5)
#   Ym_n = sum_d π_nd * Xm_d        (market clearing for manufactures in n)                     (6)
#
#   Hat Algebra:
#
#   For any change in fundamentals to (T̂_n,τ̂_od):
#
#   P̂_d = P′_d/P_d
#       = ( ∑_o T′_o * (τ′_od * W′_o^β * P′_o^(1-β) / P_d )^-θ )^(-1/θ)
#       = ( ∑_o T_o * (τ_od * W_o^β * P_o^(1-β)  / P_d )^-θ * T̂_o * (τ̂_od * Ŵ_o^β * P̂_o^(1-β) )^-θ )^(-1/θ)
#       = ( ∑_o π_od * T̂_o * (τ̂_od * Ŵ_o^β * P̂_o^(1-β) )^-θ )^(-1/θ)
#   
#   Implied trade shares
#  π′_od = π_od * π̂_od
#       = π_od * T̂_o * (τ̂_od * Ŵ_o^β * P̂_o^(1-β) / P̂_d )^-θ
#       
#   Implied market clearing condition
#
#   Ym′_o = ∑_d π′_od * Xm′_d
#   Ym′_o = ∑_d π′_od * ( Ym′_d + Dm′_d )
#   Ym_o = ∑_d π′_od * ( Ym′_d + Dm′_d )

struct DEK{T}
    Π::Matrix{T}
    Y::Vector{T}
    D::Vector{T}
    Dm::Vector{T}
    α::T
    β::T
    θ::T
end

# # for testing
# m = DEK(Π,Y,D,Dm,.5,.5,5.)
# Ŵ = rand(Uniform(.9,1.1),N)
# T̂ = rand(Uniform(.9,1.1),N)
# τ̂ = rand(Uniform(.9,1.1),N,N)
# for n = 1:N
#     τ̂[n,n] = 1.
# end
# tol=1e-16;maxit=1e4;report=true

function prices(m::DEK{T},Ŵ::Vector{T},T̂::Vector{T},τ̂::Matrix{T},;tol=1e-16,maxit=1e4,report=false) where {T <:Number}
    P̂ = ones(length(m.Y),1)
    done = false
    iter = 0
    while !done
        iter += 1
        P̂old = copy(P̂)
        P̂ = (sum(m.Π .* T̂ .* ( τ̂ .* Ŵ.^m.β .* P̂old.^(1 - m.β) ).^-m.θ,dims=1)').^-(1/m.θ)
        err = maximum(abs.(P̂ .- P̂old))
        done = (err < tol) || (iter ≥ maxit)
        if report
            @show iter,err
            # @show P̂
        end
        (iter ≥ maxit) && @warn "Maximum iterations reached when solving for prices."
    end
    return P̂[:]
end

# P̂ = prices(m,Ŵ,T̂,τ̂,report=true)

function tradeShares(m::DEK{T},P̂::Vector{T},Ŵ::Vector{T},T̂::Vector{T},τ̂::Matrix{T}) where {T <:Number}
    out = m.Π .* T̂ .* ( τ̂ .*  Ŵ.^m.β .* P̂.^(1 - m.β) ./ P̂' ).^-m.θ
    return out # ./ sum(out,dims=1) technically don't need normalization
end

function excessDemand(m::DEK{T},Ŵ::Vector{T},T̂::Vector{T},τ̂::Matrix{T},
                        D′::Vector{T},Dm′::Vector{T}) where {T<:Number}
    P̂ = prices(m,Ŵ,T̂,τ̂)
    Π′ = tradeShares(m,P̂,Ŵ,T̂,τ̂)
    Y′ = m.Y .* Ŵ
    X′ = Y′ + D′

    # Xm′ = α*X′ +(1-β)*Ym′
    # Xm′ = Ym′ + Dm′
    # ⟹ α*X′ +(1-β)*Ym′ = Ym′ + Dm′ ⟺ α*X′ - Dm′ = β*Ym′
    Ym′ = (m.α*X′ - Dm′)/m.β

    # Ym′ = Π′*Xm′ = Π*( Ym′ + Dm′ )
    return Π′*(Ym′ + Dm′) - Ym′
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
