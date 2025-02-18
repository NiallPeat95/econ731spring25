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
#   Share of d expenditure sourced from o amongst traded goods (manufactures):
#   π_od = T_o * (τ_od * W_o^β * P_o^(1-β))^-θ / sum_o( T_o * (τ_od * W_o^β * P_o^(1-β))^-θ )   (1)
#
#   Price level of traded goods in d:
#   P_d = sum_o( T_o * (τ_od * W_o^β * P_o^(1-β))^-θ )^(-1/θ)                                   (2)
#
#   For each n = 1,...,N we have:
#
#   X_n = Y_n + D_n                 (definition of trade deficit)                               (3)
#   Xm_n = Ym_n + Dm_n              (definition of manufacturing trade deficit)                 (4)
#   Xm_n = α*X_n + (1-β)*Ym_n       (accounting identity for expenditure on manufacturing)      (5)
#   Ym_n = sum_d π_nd * Xm_d        (market clearing for manufactures in n)                     (6)
#
#   Note that:
#   
#   
#

struct DEK{T}
    Π::Matrix{T}
    Y::Vector{T}
    D::Vector{T}
    Dm::Vector{T}
    α::T
    β::T
    θ::T
end
function prices(m::DEK{T},Ŵ::Vector{T};tol=1e-16,maxit=1e4,report=false) where {T <:Number}
    P̂ = ones(length(m.Y),1)
    done = false
    iter = 0
    while !done
        iter += 1
        P̂old = copy(P̂)
        P̂ = (sum(m.Π .* ( Ŵ.^m.β .* P̂old.^(1 - m.β) ).^-m.θ,dims=1)').^-(1/m.θ)
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
function tradeShares(m::DEK{T},P̂::Vector{T},Ŵ::Vector{T}) where {T <:Number}
    out = m.Π .* ( Ŵ.^m.β .* P̂.^(1 - m.β) ).^-m.θ
    return out ./ sum(out,dims=1)
end
function excessDemand(m::DEK{T},Ŵ::Vector{T},D′::Vector{T},Dm′::Vector{T}) where {T<:Number}
    P̂ = prices(m,Ŵ)
    Π = tradeShares(m,P̂,Ŵ)
    Y′ = m.Y .* Ŵ
    X′ = Y′ + D′

    # Xm = α*X +(1-β)*Ym
    # Xm = Ym + Dm
    Ym′ = (m.α*X′ - Dm′)/m.β

    # Ym = Π*Xm
    return Π*(Ym′ + Dm′) - Ym′
end
function tâtonnment(m::DEK{T},D′::Vector{T},Dm′::Vector{T};λ = T(.1),tol=1e-10,maxit=1e4,report=false,reportrate=1) where {T<:Number}
    Ŵ = ones(length(m.Y))
    done = false
    iter = 0
    t0 = time()
    while !done
        iter += 1
        Ŵold = copy(Ŵ)
        ed = excessDemand(m,Ŵ,D′,Dm′)
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
