using Graphs, SparseArrays, LinearAlgebra

import Base: copy
import Graphs: nv, has_edge, add_edge!, rem_edge!, rem_vertex!,
    rem_vertices!, add_vertex!, add_vertices!, outneighbors, inneighbors, neighbors,
    vertices, adjacency_matrix, ne, is_directed, degree, indegree, outdegree, edges,
    has_vertex, all_neighbors, edgetype

export WeightedDiMultigraph, AbstractWeightedMultigraph, add_vertex!, add_edge!, add_vertices!,
       inneighbors, outneighbors, rem_vertices!, rem_vertex!, vertices, nv, ne, has_edge, has_vertex,
       degree, indegree, outdegree, adjacency_matrix, is_directed, edges, all_neighbors, edgetype,
       rem_edge!


abstract type AbstractWeightedMultigraph{T, U} <: AbstractMultigraph{T} end

edgetype(mg::AbstractWeightedMultigraph) = WeightedMultipleEdge{eltype(mg), multype(mg)}

mutable struct WeightedDiMultigraph{T<:Integer, U<:Any} <: AbstractWeightedMultigraph{T, U}
    adjlist::Dict{T, Vector{T}}
    weights::Dict{T, Vector{U}}
    _idmax::T

    function WeightedDiMultigraph{T, U}(d::Dict{T, Vector{T}}, w::Dict{T, Vector{U}}, _idmax::T) where T where U
        adjlist = deepcopy(d)
        weights = deepcopy(w)
        vs = keys(adjlist)
        for (v, l) in adjlist
            if l ⊆ vs
                idx = sortperm(l)
                weights[v] = weights[v][idx]
                l = l[idx]
            else
                error("Some vertices connected to $v is not in the multigraph!")
            end
        end
        _idmax >= 0 || (_idmax = isempty(adjlist) ? 0 : maximum(vs))
        new{T, U}(adjlist, weights, _idmax)
    end
    
end

WeightedDiMultigraph(adjlist::Dict{T, Vector{T}}, weights::Dict{T, Vector{U}}) where {T<:Integer, U<:Any} = WeightedDiMultigraph{T, U}(adjlist, weights, isempty(adjlist) ? 0 : maximum(keys(adjlist)))

function WeightedDiMultigraph(adjmx::AbstractMatrix{Vector{U}}) where {U<:Any}
    m, n = size(adjmx)
    if m != n
        error("Adjacency matrices should be square!")
    end
    adjlist = Dict(zip((1:m), [Int[] for _ = 1:m]))
    weights = Dict(zip((1:m), [U[] for _ = 1:m]))
    for v1 = 1:m
        for v2 = 1:m
            ws = adjmx[v1, v2]
            for i = 1:length(ws)
                push!(adjlist[v1], v2)
                push!(weights[v1], ws[i])
            end
        end
    end
    WeightedDiMultigraph{Int, U}(adjlist, weights, m)
end

function WeightedDiMultigraph(n::T, ::Type{U}) where {T<:Integer, U<:Any} 
    n >= 0 || error("Number of vertices should be non-negative")
    adjlist = Dict{T, Vector{T}}()
    weights = Dict{T, Vector{U}}()

    for i = 1:n
        adjlist[i] = T[]
        weights[i] = U[]
    end
    return WeightedDiMultigraph(adjlist, weights)
end
WeightedDiMultigraph(g::SimpleDiGraph{T}, w::Type{U}) where {T<:Integer, U<:Any} = WeightedDiMultigraph(
                                                                                 Dict(zip(T(1):nv(g), Graphs.SimpleGraphs.fadj(g))), 
                                                                                 Dict(zip(T(1):nv(g), ones(w, length.(Graphs.SimpleGraphs.fadj(g)))))
                                                                                )

copy(mg::WeightedDiMultigraph{T}) where {T} = WeightedDiMultigraph{T}(deepcopy(mg.adjlist), deepcopy(mg.weights), mg._idmax)

findparam(::Dict{<:Integer, Vector{U}}) where U = U
weighttype(mg::WeightedDiMultigraph) = findparam(mg.weights)

nv(mg::WeightedDiMultigraph{T}) where {T<:Integer} = T(length(mg.adjlist))
vertices(mg::WeightedDiMultigraph) = collect(keys(mg.adjlist))
weights(mg::WeightedDiMultigraph) = mg.weights
has_vertex(mg::WeightedDiMultigraph, v::Integer) = haskey(mg.adjlist, v)

function adjacency_matrix(mg::WeightedDiMultigraph)
    wtype = weighttype(mg)
    adjmx = [wtype[] for _ in CartesianIndices((nv(mg), nv(mg)))]

    ids = sort!(vertices(mg))
    for id1 in ids
        v1 = searchsortedfirst(ids, id1)
        for (i, id2) in enumerate(mg.adjlist[id1])
            v2 = searchsortedfirst(ids, id2)
            @inbounds push!(adjmx[v1, v2], mg.weights[v1][i])
        end
    end
    return adjmx
end

function add_edge!(mg::WeightedDiMultigraph, me::AbstractWeightedMultipleEdge)
    s = src(me)
    d = dst(me)
    m = mul(me)
    ws = weights(me)
    if has_vertex(mg, s) && has_vertex(mg, d)
        for i = 1:m
            idx = searchsortedfirst(mg.adjlist[s], d)
            insert!(mg.adjlist[s], idx, d)
            insert!(mg.weights[s], idx, ws[i])
        end
        return true
    end
    return false
end

has_edge(mg::AbstractWeightedMultigraph, t::Union{NTuple{2}, Pair}) = has_edge(mg, WeightedMultipleEdge(t)) 
has_edge(mg::AbstractWeightedMultigraph, src, dst) = has_edge(mg, WeightedMultipleEdge(src, dst)) 

add_edge!(mg::AbstractWeightedMultigraph, x, y) = add_edge!(mg, WeightedMultipleEdge(x, y))
add_edge!(mg::AbstractWeightedMultigraph, x, y, z) = add_edge!(mg, WeightedMultipleEdge(x, y, z))
add_edge!(mg::AbstractWeightedMultigraph, x::Vector{<:Integer}) = add_edge!(mg, WeightedMultipleEdge(x...))

function rem_edge!(mg::WeightedDiMultigraph, me::AbstractMultipleEdge)
    if has_edge(mg, me)
        s = src(me)
        d = dst(me)
        m = mul(me)
        for i = 1:m
            idx = searchsortedfirst(mg.adjlist[s], d)
            deleteat!(mg.adjlist[s], idx)
            deleteat!(mg.weights[s], idx)
        end
        return true
    else
        return false
    end
end

function rem_vertices!(mg::WeightedDiMultigraph{T}, vs::Vector{T}) where {T<:Integer}
    if all(has_vertex(mg, v) for v in vs)
        for v in vs
            for u in neighbors(mg, v)
                l = mg.adjlist[u]
                idx = searchsorted(l, v)
                deleteat!(l, idx)
                deleteat!(mg.weights[u], idx)
            end
            delete!(mg.adjlist, v)
            delete!(mg.weights, v)
        end
        if mg._idmax in vs
            mg._idmax = maximum(keys(mg.adjlist))
        end
        return true
    end
    return false
end

function add_vertices!(mg::WeightedDiMultigraph{T, U}, n::Integer) where {T<:Integer, U<:Any}
    idmax = mg._idmax
    mg._idmax += n
    new_ids = collect((idmax+1):(idmax+n))
    for i in new_ids
        mg.adjlist[i] = T[]
        mg.weights[i] = U[]
    end
    return new_ids
end

add_vertex!(mg::WeightedDiMultigraph{T, U}) where {T<:Integer, U<:Any} = add_vertices!(mg, one(T))

function outneighbors(mg::WeightedDiMultigraph, v::Integer; count_mul::Bool = false)
    has_vertex(mg, v) || error("Vertex not found!")
    if count_mul
        return copy(mg.adjlist[v])
    else
        return sort!(collect(Set(mg.adjlist[v])))
    end
end
function inneighbors(mg::WeightedDiMultigraph{T}, v::Integer; count_mul::Bool = false) where T
    has_vertex(mg, v) || error("Vertex not found!")
    
    innb = T[]
    for u in vertices(mg)
        mul_u_v = length(searchsorted(outneighbors(mg, u), v))
        if mul_u_v > 0
            if count_mul
                for i = 1:mul_u_v
                    push!(innb, u)
                end
            else
                push!(innb, u)
            end
        end
    end
    sort!(innb)
    return innb
end
neighbors(mg::WeightedDiMultigraph, v::Integer; count_mul::Bool = false) = outneighbors(mg, v, count_mul = count_mul)
all_neighbors(mg::WeightedDiMultigraph, v::Integer) = sort!(union(outneighbors(mg, v), inneighbors(mg, v)))

function mul(mg::WeightedDiMultigraph, s::Integer, d::Integer)
    (has_vertex(mg, s) && has_vertex(mg, d)) || error("Vertices not found!")
    return length(searchsorted(mg.adjlist[s], d))
end

function weights(mg::WeightedDiMultigraph, s::Integer, d::Integer)
    (has_vertex(mg, s) && has_vertex(mg, d)) || error("Vertices not found!")
    return mg.weights[s][searchsorted(mg.adjlist[s], d)]
end

is_directed(mg::WeightedDiMultigraph) = true

function ne(mg::WeightedDiMultigraph; count_mul::Bool = false)
    if count_mul
        return sum([length(mg.adjlist[v]) for v in vertices(mg)])
    else
        return sum([length(Set(mg.adjlist[v])) for v in vertices(mg)])
    end
end

function outdegree(mg::WeightedDiMultigraph{T}) where T
    degs = Dict{T, Int}()
    for v in vertices(mg)
        degs[v] = length(mg.adjlist[v])
    end
    return degs
end
function indegree(mg::WeightedDiMultigraph{T}) where T
    degs = Dict{T, Int}()
    for v in vertices(mg)
        degs[v] = 0
    end
    for v in vertices(mg)
        for u in mg.adjlist[v]
            degs[u] += 1
        end
    end
    return degs
end
degree(mg::WeightedDiMultigraph) = outdegree(mg)
degree(mg::WeightedDiMultigraph, v::Integer) = outdegree(mg, v)
indegree(mg::WeightedDiMultigraph, v::Integer) = length(inneighbors(mg, v; count_mul = true))
outdegree(mg::WeightedDiMultigraph, v::Integer) = length(outneighbors(mg, v; count_mul = true))
