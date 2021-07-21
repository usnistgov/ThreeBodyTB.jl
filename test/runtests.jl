
println("using plot")

using Plots

println("try plot")

plot([0,1,2,3], [0.1, 0.5, 0.5, 0.1])
println("did plot")
println("savefig")
savefig("sf.pdf")
println("did savefig ")

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


Nothing
