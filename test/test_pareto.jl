@testset "paretoset" begin
    p = hcat(([x,y] for x in 0:0.2:1 for y in 0:0.2:1 if (x-1)^2+(y-1)^2<=1)...)
    @test p[:,EGO.ParetoSet(p)] == [0.0 0.2 0.4 1.0; 1.0 0.4 0.2 0.0]
    @test p[:,EGO.ParetoSet(p,[:Min, :Min])] == [0.0 0.2 0.4 1.0; 1.0 0.4 0.2 0.0]
    @test p[:,EGO.ParetoSet(p,[EGO.Min, EGO.Min])] == [0.0 0.2 0.4 1.0; 1.0 0.4 0.2 0.0]
    p = hcat(([x,y] for x in 0:0.2:1 for y in 0:0.2:1 if x^2+y^2<=1)...)
    @test p[:,EGO.ParetoSet(p, [:Max, :Max])] == [0.0 0.6 0.8 1.0; 1.0 0.8 0.6 0.0]
    @test p[:,EGO.ParetoSet(p, [EGO.Max, EGO.Max])] == [0.0 0.6 0.8 1.0; 1.0 0.8 0.6 0.0]
end
