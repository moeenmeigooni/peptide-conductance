Python code for MD simulations and conformational analysis of peptides in STM break junction corresponding to the following paper:

R. Samajdar, M. Meigooni, H. Yang, J. Li, X. Liu, N. E. Jackson, M. A. Mosquera, E. Tajkhorshid, and C. M. Schroeder (2024). Secondary structure determines electron transport in peptides. *Proceedings of the National Academy of Sciences*, 121(32):e2403324121.


Directories are organized by amino acid composition. 


`build` folders contain all files necessary for creating 16 replicates of the given peptide with randomized phi,psi angles. `build.sh` runs necessary scripts. Packages/softwares used include:
 - VMD 1.9.4a54 (for tcl scripts)
 - PeptideBuilder 1.1.0


`hold*-orient-efield` folder contains python script for OpenMM simulations with custom potentials for the given holding stage. Run scripts provided (`run_hold.sh`) are configured for running on DGX-2 with 16xV100 GPUs. Packages/softwares used include:
 - OpenMM 7.7.0
 - MDTraj 1.9.3


`analyze-extract.ipynb` notebook contains code for loading/analyzing trajectories and extracting relevant conformers for NEGF/DFT calculations, and for generating plots. Packages used include:
 - MDTraj 1.9.3 for loading trajectories
 - NumPy 1.21.5 and custom functions (jit-accelerated with numba 0.53.1) for analysis e.g. computing distances and angles
 - scikit-learn 1.0.2 for PCA
 - PyEMMA 2.5.11 for free energy plots
 - Matplotlib 3.5.2 for all other plots


Solvent-stripped trajectory files are available at doi.org/10.5281/zenodo.7843691
