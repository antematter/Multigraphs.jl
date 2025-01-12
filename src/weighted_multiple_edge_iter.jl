using Graphs

import Base: eltype, iterate, length

export WeightedMultipleEdgeIter, edges

edges(mg::AbstractWeightedMultigraph) = WeightedMultipleEdgeIter(mg)

struct WeightedMultipleEdgeIter{G<:AbstractWeightedMultigraph} <:AbstractEdgeIter
    g::G
end

function iterate(eit::WeightedMultipleEdgeIter{G}, state=(one(eltype(eit.g)), one(eltype(eit.g)))) where {G <: AbstractWeightedMultigraph}
    g = eit.g
    n = nv(g)
    vs = vertices(g)
    sort!(vs)
    u, i = state

    @inbounds while u <= n
        list_u = outneighbors(g, vs[u])
        if i > length(list_u)
            u += 1
            i = one(u)
            continue
        end
        if is_directed(g)
            e = WeightedMultipleEdge(vs[u], list_u[i], mul(g, vs[u], list_u[i]), weights(g, vs[u], list_u[i]))
            state = (u, i + 1)
            return e, state
        else
            if list_u[i] >= vs[u]
                e = WeightedMultipleEdge(vs[u], list_u[i], mul(g, vs[u], list_u[i]), weights(g, vs[u], list_u[i]))
                state = (u, i + 1)
                return e, state
            else
                i += 1
            end
        end
    end

    if n == 0 || u > n
        return nothing
    end
end

length(eit::WeightedMultipleEdgeIter) = ne(eit.g)