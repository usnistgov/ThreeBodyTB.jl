using ThreeBodyTB


#precompile all codes used in examples

println("precompile")

#include(joinpath( ThreeBodyTB.TESTDIR ,"test_examples.jl"))


for t in ["calculate_energy_plot_fcc_al.jl", "load_dft_plot.jl", "load_tbc_plot.jl", "plot_Sc_P.jl", "relax_fcc_al.jl"]
    include(joinpath( ThreeBodyTB.EXAMPLESDIR ,t))
end


