"""
    Multigraphs

An multigraph extension for `Graphs`.
"""
module Multigraphs

include("multiple_edge.jl")
include("abstract_multigraph.jl")
include("multiple_edge_iter.jl")
include("multigraph_adjlist.jl")
include("di_multigraph_adjlist.jl")
include("weighted_multiple_edge.jl")
include("di_weighted_multigraph_adjlist.jl")

end # module
