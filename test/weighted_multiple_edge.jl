import Base.one 


struct ValidWeight end 
ValidWeight() = 1.0

wme = WeightedMultipleEdge(1, 2, 3, ValidWeight)

try 
    WeightedMultipleEdge(1, 2, 1, 1.0)
catch err 
    @test err !== nothing 
end

try
    WeightedMultipleEdge(1, 2, 0)
catch err
    @test err !== nothing
end

@test src(wme) == 1 && dst(wme) == 2 && mul(wme) == 3 && weights(wme) == [1.0, 1.0, 1.0]
e0 = WeightedMultipleEdge([1, 2, 1], [1.0])
@test WeightedMultipleEdge(1, 2, 1, Float64) == e0
@test e0 == WeightedMultipleEdge((1, 2))
@test e0 == WeightedMultipleEdge(1 => 2)
@test reverse(wme) == WeightedMultipleEdge(2, 1, 3, ValidWeight)
@test eltype(wme) == Int

@test iterate(wme)[2] == 2
@test [e0 == e for e in wme] == [true for i = 1:mul(wme)]
@test Tuple(wme) == (1,2,3, [1.0, 1.0, 1.0])
length(wme) == mul(wme)