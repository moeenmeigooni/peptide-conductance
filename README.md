# peptide-conductance
Python scripts for MD simulations and conformational analysis of peptides in STM break junction.


Directories organized by amino acid composition. Run scripts (`run_hold.sh`) are configured for running on DGX-2 with 16xV100 GPUs.


`build` folder contains all files necessary for creating 16 replicates of the given peptide with randomized phi,psi angles. `build.sh` runs necessary scripts. Packages/softwares used include:
 - VMD 1.9.4a54 (for tcl scripts)
 - PeptideBuilder 1.1.0


`hold*-orient-efield` folder contains python script for OpenMM simulations with custom potentials for the given holding stage. Packages/softwares used include:
 - OpenMM 7.7.0
 - MDTraj 1.9.3


`analyze-extract.ipynb` notebook contains code for loading/analyzing trajectories and extracting relevant conformers for NEGF/DFT calculations, and for generating plots. Packages used include:
 - mdtraj 1.9.3 for loading trajectories
 - numpy 1.21.5 and custom functions (jit-accelerated with numba 0.53.1) for analysis e.g. computing distances and angles
 - scikit-learn 1.0.2 for PCA
 - pyemma 2.5.11 for free energy plots
 - matplotlib 3.5.2 for all other plots
