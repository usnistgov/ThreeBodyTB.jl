using ThreeBodyTB

println("Precompiling, this may take a while...")

ThreeBodyTB.my_precompile()

println("Done Precompiling.")
println()


##
###precompile all codes used in examples
##
##println("Precompiling, this may take a while...")
##
###include(joinpath( ThreeBodyTB.TESTDIR ,"test_examples.jl"))
##
##
##ThreeBodyTB.BandStruct.set_no_display(true) #do not display plots
##
##tryusing(pkgsym) = try
##    @eval using $pkgsym
##    return true
##catch e
##    return e
##end
##
##
##if tryusing(:Suppressor) #use Suppressor to avoid excessive output
##
##    using $Suppressor
##    @suppress for t in ["calculate_energy_plot_fcc_al.jl", "load_dft_plot.jl", "load_tbc_plot.jl", "plot_Sc_P.jl", "relax_fcc_al.jl"]
##        include(joinpath( ThreeBodyTB.EXAMPLESDIR ,t))
##    end
##    
##else #Suppressor is missing
##    
##    for t in ["calculate_energy_plot_fcc_al.jl", "load_dft_plot.jl", "load_tbc_plot.jl", "plot_Sc_P.jl", "relax_fcc_al.jl"]
##        include(joinpath( ThreeBodyTB.EXAMPLESDIR ,t))
##    end
##    
##end    
##
##ThreeBodyTB.BandStruct.set_no_display(false) #default behavior


