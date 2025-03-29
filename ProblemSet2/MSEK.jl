#
#   Defines the MSEK type for representing a multisector EK model as in 
#       Caliendo Parro (2015).
#
#   Nels Lind, 2/21/2025
#
#   j = 1,…,J, sectors

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
    Π::Array{T,3}
    Y::Vector{T}
    D::Vector{T}
    α::Matrix{T,2}
    θ::Vector{T} # Trade elasticities for each sector j
    μ::Matrix{T,2} # Preference parameters for final demand in each sector j and country n
    Π_l::Array{T,3}
    v::Array{T}
end

function prices(m::MSEK{T},Ŵ::Vector{T},T̂::Matrix{T}, μ̂::Array{T,2}, ;tol=1e-16,maxit=1e4,report=false) where {T <:Number}
 
    P̂ = ones(size(m.μ)...)
    done = false
    iter = 0
    while !done
        iter += 1
        P̂old = copy(P̂)
        Ĉ = Ŵ'
        P̂ = dsum( m.Π .* T̂ .* ( Ĉ ).^(.-m.θ),dims=2).^(.- 1 ./ m.θ)
        Ŵ = dsum( m.Π_l .* exp(μ̂).* (Ŵ_sn).^m.v,dims=2).^(1/m.v)
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
    out = m.Π .* T̂ .* ( Ĉ ./ addDim(P̂,2) ).^(.-m.θ)
    return out ./ sum(out,dims=2)
end

function laborshares(m::MSEK{T},P̂::Matrix{T},Ŵ::Vector{T}) where {T <:Number}
    out = m.Π_l .* exp(m.μ̂).* (addDim(Ŵ,2)).^m.v
    return out ./ sum(out,dims=2)
end

function excessDemand(m::MSEK{T},Ŵ::Vector{T},T̂::Matrix{T},τ̂::Array{T,3},
                        t′::Array{T,3},D′::Vector{T}) where {T<:Number}
    P̂ = prices(m,Ŵ,T̂,τ̂,t′)
    Π′ = tradeShares(m,P̂,Ŵ,T̂,τ̂,t′)
    Π̃′ = Π′ 
    Π_l′ = laborshares(m,P̂,Ŵ)

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
    A1 = blockmatrix([ γ[:,:,n] .* addDim(Π̃′[:,n,d],1) for n=1:N, d=1:N])
    
    # Next, consider μ_in * ∑_jo t′_jon π̃′_jon * X′_jn. The matrix representation
    # of this linear operator has in × jd element of μ_in * ∑_o t′_jon π̃′_jon. Note
    # that this value doesn't depend on d. The (n,d)th J × J submatrix is
    #
    # μ[:,n] * sum(t′[:,:,n].*Π̃′[:,:,n],dims=2)'
    A2 = blockmatrix([ α[:,n] * sum(t′[:,:,n].*Π̃′[:,:,n],dims=2)' for n=1:N, d=1:N])

    vX′ = (I-A1-A2)\vec(α .* ( Ŵ.*m.Y + D′ )')

#   X′_jod = π̃′_jod * X′_jd 
#   excessDemand = ∑_jd X′_jnd + D′_n - ∑_jo X′_jon
X′ = Π′ .* (m.μ .* (Ŵ .* m.Y .+ D′)')
return sum(X′, dims=(1,2)) + D′ - sum(X′, dims=(1,3))
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
