using Multigraphs, Graphs

samplezeros(t, x, y) = [t[] for _ in CartesianIndices((x, y))]

try
    m2 = samplezeros(AbstractFloat, 2, 3)
    dg = WeightedDiMultigraph(m2)
catch e
    @test e !== nothing
end
try
    m2 = samplezeros(Int, 2, 2)
    m2[1, 2] = [2, 2]
    dg = WeightedDiMultigraph(m2)
catch e
    @test e !== nothing
end

m = samplezeros(Int, 4, 4)
m
m[1,2] = [2]
m[2,1] = [2]
m[2,3] = [2]
m[3,2] = [2]

g = WeightedDiMultigraph(m)
g = WeightedDiMultigraph(Matrix(m))
g0 = WeightedDiMultigraph(2, Int)
@test !add_edge!(g0, 2, 3) && !rem_edge!(g0, 1, 2)
@test adjacency_matrix(g) == m

@test is_directed(g)
@test edgetype(g) == WeightedMultipleEdge{Int, Int}
@test size(adjacency_matrix(g), 1) == 4

@test nv(g) == 4 && ne(g, count_mul = true) == 4 && ne(g) == 4

add_vertices!(g, 3)
@test nv(g) == 7

@test has_edge(g, 1, 2, 1)
@test rem_vertices!(g, [7, 5, 4, 6])
add_edge!(g, [2, 3, 2])
rem_edge!(g, [2, 3, 2])
add_edge!(g, 2, 3)
rem_edge!(g, 2, 3)
add_edge!(g, 2, 3, 2)
rem_edge!(g, 2, 3, 1)

@test has_edge(g, 2, 3) && has_edge(g, (2, 3))
@test !has_edge(g, 3, 4) && !has_edge(g, 2, 5)
@test has_vertex(g, 1) && !has_vertex(g, 5)
for v in vertices(g)
    @test inneighbors(g, v) == outneighbors(g, v)
    @test degree(g, v) == indegree(g, v) && indegree(g, v) == outdegree(g, v)
end
add_vertex!(g)
@test indegree(g) == outdegree(g)

mg0 = WeightedDiMultigraph(0)
@test nv(mg0) == ne(mg0) == 0