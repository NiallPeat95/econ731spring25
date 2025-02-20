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
#   π_jod = T_jo ( τ_jon C_jo )^-θ_j / P_jn )^(-1/θ_j)
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
#   C_jn = W_n^(1-1 - ∑_k γ_jkn) * ∏_k P_kn^γ_jkn
#
#   Expenditure on sector j by n
#
#   X_jn = ∑_k γ_jkn ∑_i X_in

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
