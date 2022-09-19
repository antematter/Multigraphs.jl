# Multigraphs

[![CI](https://github.com/QuantumBFS/Multigraphs.jl/workflows/CI/badge.svg)](https://github.com/QuantumBFS/Multigraphs.jl/actions)
[![Codecov](https://codecov.io/gh/QuantumBFS/Multigraphs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/QuantumBFS/Multigraphs.jl)

Multigraphs extension for `Graphs.jl`.

## Note

This fork adds a weighted version of the Multigraphs. The usage is as follows:

```julia
julia> wmg = WeightedDiMultigraph(4, Float64)
{4, 0} directed Int64 multigraph

julia> add_edge!(wmg, 1, 2)
true

julia> add_edge!(wmg, 3, 4)
true

julia> add_edge!(wmg, 2, 3)
true

julia> for e in edges(wmg) println(e) end
Multiple edge 1 => 2 with multiplicity 1 and weights [1.0]
Multiple edge 2 => 3 with multiplicity 1 and weights [1.0]
Multiple edge 3 => 4 with multiplicity 1 and weights [1.0]

julia> we = WeightedMultipleEdge([1, 2, 3], [2, 3, 4])
Multiple edge 1 => 2 with multiplicity 3 and weights [2, 3, 4]

julia> add_edge!(wmg, we)
true

julia> for e in edges(wmg) println(e) end
Multiple edge 1 => 2 with multiplicity 4 and weights [4.0, 3.0, 2.0, 1.0]
Multiple edge 2 => 3 with multiplicity 1 and weights [1.0]
Multiple edge 3 => 4 with multiplicity 1 and weights [1.0]
```

## Installation

<p>
Multigraphs is a &nbsp;
    <a href="https://julialang.org">
        <img src="https://julialang.org/favicon.ico" width="16em">
        Julia Language
    </a>
    &nbsp; package. To install Multigraphs,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd> key in the REPL to use the package mode, then type the following command
</p>

```julia
pkg> add Multigraphs
```

## Examples

```julia
using Graphs, Multigraphs

# create a undirected multigraph with 3 vertices and 0 multiple edges
# use DiMultigraph for directed multigraphs
julia> mg = Multigraph(3)
{3, 0} undirected Int64 multigraph with Int64 multiplicities

# add a multiple edge from 1 to 2 with multiplicity 2
julia> add_edge!(mg, 1, 2, 2)
true

# add a simple edge (multiple edge with multiplicity 1) from 2 to 3
julia> add_edge!(mg, 2, 3)
true

# this will increase multiplicity of the edge from 2 to 3 by 2
julia> add_edge!(mg, 2, 3, 2) 
true

# this will decrease multiplicity of the edge from 2 to 3 by 1
julia> rem_edge!(mg, 2, 3, 2) 

# here me is a MultipleEdge
julia> mes = [me for me in edges(mg)]
2-element Array{MultipleEdge{Int64,Int64},1}:
Multiple edge 1 => 2 with multiplicity 2
Multiple edge 2 => 3 with multiplicity 1

# here e is a Graphs.SimpleEdge
julia> for e in mes[1] 
           println(e)
       end
Edge 1 => 2
Edge 1 => 2

```

## License

MIT License
