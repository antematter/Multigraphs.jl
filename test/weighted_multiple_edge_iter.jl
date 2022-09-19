mg = WeightedDiMultigraph(3, AbstractFloat)
add_vertices!(mg, 3)
rem_vertices!(mg, [1, 3])
add_edge!(mg, 2, 5)
add_edge!(mg, 2, 4, 2)

@test outneighbors(mg, 2) == [4, 5]
eit = edges(mg)
iterate(eit)[1]
vertices(mg)
@test iterate(eit)[2] == (1, 2)
mes = [me for me in edges(mg)]
@test length(mes) == length(eit)
