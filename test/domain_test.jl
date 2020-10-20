using QuasiGeostrophy

@testset "Domain Tests" begin
    Ω = S¹(0,2π)
    @test dim(Ω)==1
    @test dim(Ω×Ω)==2
    @test dim(Ω×Ω×Ω)==3
end