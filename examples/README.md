## Examples

1) `calculate_energy_plot_fcc_al.jl` - Create a crystal object of
face-centered-cubic Al, construct the tight-binding model, and
calculate the energy and band structure. Uses pre-fit coefficient database.

2) `relax_fcc_al.jl` - Now calculate forces/stresses and relax the fcc Al
structure.

3) `plot_Sc_P.jl` - Calculate the energy, band structure, and density
of states (DOS) for ScP in the rocksalt structure.

4) `load_tbc_plot.jl` - Load a tight binding dataset, and plot the
band structure.

5) `load_dft_plot.jl` - Load output from a DFT calculation and compare
band structure to a tight-binding model.

6) `fit_model.jl` - Fit your own tight-binding coefficients from your
own DFT calculations. Shows how to run Quantum Espresso using the
`ThreeBodyTB` code, construct tight-binding models for individual DFT
calculations, and then fit coefficients to those DFT/TB models. This
example is only for advanced users who want to fit their own coefficients,
rather than use the provided coefficients.
