using Test
include("/Users/matteopolimeno/Local/Negotium/Projects/2025_10/BirdLung/PostProcessingJuliaCodes/DataParsing.jl")
using .DataParsing


@testset "DataParsing" begin
    include("comsol_parsing_tests.jl")
end