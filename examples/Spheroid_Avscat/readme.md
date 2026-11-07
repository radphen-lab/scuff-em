## instructions
Example file of scuff-avescatter for oblate(prolate) spheroids

to run the each simulation, do:
scuff-avescatter --geometry oblate.scuffgeo --OmegaFile OmegaList.dat

Use "compare_results.py" to compare the simulation with the exact solution.
The script reads the current .SIFlux.AVSCAT output, whose third column contains
the surface label, and falls back to the supplied historical .AVSCAT.EMTPFT file.
