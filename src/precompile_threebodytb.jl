using ThreeBodyTB

#precompile all codes used in examples

println("Precompiling, this may take a while...")

#include(joinpath( ThreeBodyTB.TESTDIR ,"test_examples.jl"))

try  #use Suppressor to avoid excessive output

    using Suppressor
    @suppress for t in ["calculate_energy_plot_fcc_al.jl", "load_dft_plot.jl", "load_tbc_plot.jl", "plot_Sc_P.jl", "relax_fcc_al.jl"]
        include(joinpath( ThreeBodyTB.EXAMPLESDIR ,t))
    end
    
catch #Suppressor is missing
    
    for t in ["calculate_energy_plot_fcc_al.jl", "load_dft_plot.jl", "load_tbc_plot.jl", "plot_Sc_P.jl", "relax_fcc_al.jl"]
        include(joinpath( ThreeBodyTB.EXAMPLESDIR ,t))
    end
    
end    

