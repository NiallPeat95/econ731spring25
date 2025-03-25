#
#   Defines the MSEK type for representing a multisector EK model as in 
#       Caliendo Parro (2015).
#
#   Nels Lind, 2/21/2025
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

function prices(m::MSEK{T},Ŵ::Vector{T},T̂::Matrix{T},τ̂::Array{T,3},t′::Array{T,3};tol=1e-16,maxit=1e4,report=false) where {T <:Number}
    #   κ̂_jod = τ̂_jod + (1+t′_jod)/(1+t_jod)
    #   Ĉ_jn = Ŵ_n^(1-1 - ∑_k γ_jkn) * ∏_k P̂_kn^γ_jkn
    #   P̂_jn = ( ∑_o π_jon T̂_jo ( κ̂_jon Ĉ_jo )^-θ_j )^(-1/θ_j)
    P̂ = ones(size(m.μ)...)
    κ̂ = τ̂ 
    done = false
    iter = 0
    while !done
        iter += 1
        P̂old = copy(P̂)
        Ĉ = Ŵ'
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
    κ̂ = τ̂ 
    Ĉ = Ŵ'
    out = m.Π .* T̂ .* ( κ̂ .* Ĉ ./ addDim(P̂,2) ).^(.-m.θ)
    return out ./ sum(out,dims=2)
end

function excessDemand(m::MSEK{T},Ŵ::Vector{T},T̂::Matrix{T},τ̂::Array{T,3},
                        t′::Array{T,3},D′::Vector{T}) where {T<:Number}
    P̂ = prices(m,Ŵ,T̂,τ̂,t′)
    Π′ = tradeShares(m,P̂,Ŵ,T̂,τ̂,t′)
    Π̃′ = Π′

#   X′_in = ∑_j γ_ijn ∑_d π̃′_jnd * X′_jd  + μ_in * ( Ŵ_n*W_n*L_n + ∑_jo t′_jon π̃′_jon * X′_jn  + D′_n )

#   X′_in - ∑_j γ_ijn ∑_d π̃′_jnd * X′_jd - μ_in * ∑_jo t′_jon π̃′_jon * X′_jn 
#        =  μ_in * ( Ŵ_n*W_n*L_n + D′_n )

    # Caliendo Parro (2015) vectorize this equation to form an equation of the form
    #   A*vec(X′) = vec(μ .* ( Ŵ.*m.Y + D′ )')
    # where the jd entry of  X′ is X′_jd. Julia vecorizes column wise so the
    # jd entry of X′_jd is in the j + J*(d-1) position of vec(X′).

    # Consider ∑_j γ_ijn ∑_d π̃′_jnd * X′_jd. This is a linear operator on vec(X′)
    # whose matrix representation has in × dj element of γ_ijn*π̃_jnd. We can organize
    # this matrix into J × J submatrices where the (n,d)th submatrix is
    #
    # γ[:,:,n] .* addDim(Π̃′[:,n,d],1)
    #
    A1 = blockmatrix([ γ[:,:,n] .* addDim(Π̃′[:,n,d],1) for n=1:N, d=1:N]) # Feed-through in intermediate demand.

    
    # Next, consider μ_in * ∑_jo t′_jon π̃′_jon * X′_jn. The matrix representation
    # of this linear operator has in × jd element of μ_in * ∑_o t′_jon π̃′_jon. Note
    # that this value doesn't depend on d. The (n,d)th J × J submatrix is
    #
    # μ[:,n] * sum(t′[:,:,n].*Π̃′[:,:,n],dims=2)'
    A2 = blockmatrix([ μ[:,n] * sum(t′[:,:,n].*Π̃′[:,:,n],dims=2)' for n=1:N, d=1:N])

    vX′ = (I-A1-A2)\vec(μ .* ( Ŵ.*m.Y + D′ )')

#   X′_jod = π̃′_jod * X′_jd 
#   excessDemand = ∑_jd X′_jnd + D′_n - ∑_jo X′_jon
    X′ = Π̃′ .* addDim(reshape(vX′,J,N),2)
    return dsum(X′,dims=(1,2)) + D′ - dsum(X′,dims=(1,3))
end
function tâtonnment(m::MSEK{T},T̂::Matrix{T},τ̂::Array{T,3},t′::Array{T,3},D′::Vector{T};
        λ = T(.0001),decay=T(.1),inflate=T(.01),tol=1e-8,maxit=1e6,report=false,reportrate=1) where {T<:Number}
    Ŵ = ones(length(m.Y))
    done = false
    iter = 0
    err = Inf
    t0 = time()
    while !done
        iter += 1
        Ŵold = copy(Ŵ)
        errold = copy(err)
        ed = excessDemand(m,Ŵ,T̂,τ̂,t′,D′)
        Ŵ = Ŵold .* (1 .+ λ .* ed ./ (Ŵold .* m.Y))
        err = maximum(abs.(Ŵ .- Ŵold))
        done = (err < tol) || (iter ≥ maxit)
        if err < errold
            λ *= 1+inflate
        else
            λ *= 1-decay
        end
        if report && (time()- t0) > reportrate
            @show iter,err,λ
            # @show Ŵ
            t0 = time()
        end
        (iter ≥ maxit) && @warn "Maximum iterations reached when solving for wages."
    end
    return Ŵ
end
