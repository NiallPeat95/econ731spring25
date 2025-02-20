
offdiag(A::AbstractMatrix) = [A[m,n] for m=1:size(A,1),n=1:size(A,2) if m!=n] 
@inline dsum(a::AbstractArray;dims=:) = dropdims(sum(a,dims=dims),dims=dims)
function addDim(A::T,pos::Integer) where {T<:AbstractArray}
    sz = size(A)
    return reshape(A,(sz[1:pos-1]...,1,sz[pos:end]...))
end
function addDims(A::T,pos::Integer,count::Integer) where {T<:AbstractArray}
    sz = size(A)
    return reshape(A,(sz[1:pos-1]...,fill(1,count)...,sz[pos:end]...))
end