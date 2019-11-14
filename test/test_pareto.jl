using EfficientGlobalOptimization
@testset "paretoset" begin
    p = hcat(([x,y] for x in 0:0.2:1 for y in 0:0.2:1 if (x-1)^2+(y-1)^2<=1)...)
    @test p[:,ParetoSet(p)] == [0.0 0.2 0.4 1.0; 1.0 0.4 0.2 0.0]
end
