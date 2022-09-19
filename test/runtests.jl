using Multigraphs, Graphs, SparseArrays
using Test

@testset "multiple_edge.jl" begin
    include("multiple_edge.jl")
end

@testset "multigraph_adjlist.jl" begin
    include("multigraph_adjlist.jl")
end

@testset "multiple_edge_iter.jl" begin
    include("multiple_edge_iter.jl")
end

@testset "di_multigraph.jl" begin
    include("di_multigraph_adjlist.jl")
end

@testset "weighted_multiple_edge.jl" begin 
    include("weighted_multiple_edge.jl")
end

@testset "weighted_multiple_edge_iter.jl" begin 
    include("weighted_multiple_edge_iter.jl")
end

@testset "weighted_di_multigraph_adjlist.jl" begin 
    include("weighted_di_multigraph_adjlist.jl")
end