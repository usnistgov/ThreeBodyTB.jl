ENV["GKSwstype"] = "100" #magic fix tests on github headless machines with GR as plotting backend


using ThreeBodyTB
using Test


include("crystal_testing.jl")
include("test_ewald.jl")
include("sym_test.jl")
include("test_laguerre.jl")
include("test_U.jl")
include("test_scf.jl")
include("test_forces.jl")
include("test_proto.jl")
include("test_makescf.jl")
include("test_atomicproj.jl")


include("test_examples.jl")

delete!(ENV, "GKSwstype") #undo magic to avoid side effects


Nothing
