import Base: eltype, Pair, Tuple, show, ==, iterate, length
import Graphs: AbstractEdge, SimpleEdge, src, dst, reverse, weights

export AbstractWeightedMultipleEdge, WeightedMultipleEdge, mul, weights

"""
    AbstractWeightedMultipleEdge{T, U, W} <: AbstractMultipleEdge{T, U}

An abstract type representing weighted multiple edges. 
"""
abstract type AbstractWeightedMultipleEdge{T, U, W} <: AbstractMultipleEdge{T, U} end

"""
    WeightedMultipleEdge{T, U, W} <: AbstractWeightedMultipleEdge{T, U, W}

A struct representing weighted multiple edges.

## Examples
```jltestdoc
julia> using Graphs, Multigraphs

julia> me = WeightedMultipleEdge(1, 2, 3, [1.0, 2.0, 3.0])
Multiple edge 1 => 2 with multiplicity 3 with weights [1.0, 2.0, 3.0]

julia> for e in me println(e) end
Edge 1 => 2 with weight 1.0
Edge 1 => 2 with weight 2.0    
Edge 1 => 2 with weight 3.0

```
"""
struct WeightedMultipleEdge{T<:Integer, U<:Integer, W<:Any} <: AbstractWeightedMultipleEdge{T, U, W}
    src::T
    dst::T
    mul::U
    ws::Vector{W}
    function WeightedMultipleEdge(src::T, dst::T, mul::U, weights::Vector{W}) where {T<:Integer, U<:Integer, W<:Any}
        if mul > 0
            @assert length(weights) == mul || "weights should match edge multiplicity"
            return new{T, U, W}(src, dst, mul, weights)
        else
            error("a multiple edge should have positive multiplicity")
        end
    end
end

function WeightedMultipleEdge(src::T, dst::T, mul::U, w::Type{W}) where {T<:Integer, U<:Integer, W<:Any}
    if isstructtype(w) 
        return WeightedMultipleEdge(src, dst, mul, [w() for i = 1:mul])
    end 
    return WeightedMultipleEdge(src, dst, mul, [one(w) for i = 1:mul])
end

WeightedMultipleEdge(src::T, dst::T, ws::U) where {T<:Integer, U<:Any} = WeightedMultipleEdge(src, dst, one(Int), [ws])
WeightedMultipleEdge(src::T, dst::T) where {T<:Integer} = WeightedMultipleEdge(src, dst, one(Int))
eltype(e::T) where {T<:AbstractWeightedMultipleEdge} = eltype(src(e))

src(e::WeightedMultipleEdge) = e.src
dst(e::WeightedMultipleEdge) = e.dst


"""
    mul(e)

Return the multiplicity of the weighted multiple edge `e`.

## Examples
```jltestdoc
julia> using Graphs, Multigraphs

julia> me = WeightedMultipleEdge(1, 2, 3, [1.0, 2.0, 3.0])
Multiple edge 1 => 2 with multiplicity 3 with weights [1.0, 2.0, 3.0]

julia> mul(me)
3

"""
mul(e::WeightedMultipleEdge) = e.mul

"""
    weights(e)

Return the weights of the weighted multiple edge `e`.

## Examples
```jltestdoc
julia> using Graphs, Multigraphs

julia> me = WeightedMultipleEdge(1, 2, 3, [1.0, 2.0, 3.0])
Multiple edge 1 => 2 with multiplicity 3 with weights [1.0, 2.0, 3.0]

julia> weights(me)
3-element Vector{Float64}:
 1.0
 2.0
 3.0

"""
weights(e::AbstractWeightedMultipleEdge) = e.ws

show(io::IO, e::AbstractWeightedMultipleEdge) = print(io, "Multiple edge $(src(e)) => $(dst(e)) with multiplicity $(mul(e)) and weights $(weights(e))")

Tuple(e::AbstractWeightedMultipleEdge) = (src(e), dst(e), mul(e), weights(e))


reverse(e::T) where {T<:AbstractWeightedMultipleEdge} = WeightedMultipleEdge(dst(e), src(e), mul(e), reverse(weights(e)))


struct WeightedSimpleEdge{U<:Any}
    src::Int
    dst::Int
    w::U
end

src(e::WeightedSimpleEdge) = e.src
dst(e::WeightedSimpleEdge) = e.dst

WeightedSimpleEdge(e::AbstractWeightedMultipleEdge, state::U) where U<:Integer = WeightedSimpleEdge(src(e), dst(e), weights(e)[state - 1]) 

==(e1::AbstractWeightedMultipleEdge, e2::AbstractWeightedMultipleEdge) = (src(e1) == src(e2) && dst(e1) == dst(e2) && mul(e1) == mul(e2) && weights(e1) == weights(e2))
==(e1::AbstractWeightedMultipleEdge, e2::WeightedSimpleEdge) = (src(e1) == src(e2) && dst(e1) == dst(e2) && any(==(e2.w), weights(e1)))
==(e1::WeightedSimpleEdge, e2::AbstractWeightedMultipleEdge) = (src(e1) == src(e2) && dst(e1) == dst(e2) && any(==(e1.w), weights(e2)))

show(io::IO, e::WeightedSimpleEdge) = print(io, "Edge $(src(e)) => $(src(e)) with weight $(e.w)")

function iterate(e::WeightedMultipleEdge{T, U, W}, state::U=one(U)) where {T<:Integer, U<:Integer, W<:Any}
    if state > mul(e)
        return nothing
    else
        state += one(U)
        return (WeightedSimpleEdge(e, state), state)
    end
end

length(me::WeightedMultipleEdge) = mul(me)