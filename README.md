Radial Distribution Function (RDF) Calculation for CP2K

This script calculates the Radial Distribution Function (RDF) for a given set of atomic positions from XYZ files. The script reads the atomic positions, applies periodic boundary conditions, calculates the distances between specified atom types, and generates the RDF. The results are saved as a plot and a CSV file.

Requirements:
Python 3.x
NumPy
Matplotlib
Glob
Multiprocessing
TQDM

You can install the required packages using pip:
pip install numpy matplotlib glob2 tqdm

Parameters:
file_pattern: Pattern to match XYZ files (e.g., 'water-traj.xyz-pos-1.xyz')
inp_file: Input file containing cell lengths (e.g., 'water.inp')
corrected_xyz_file: Output file for corrected atomic positions (e.g., 'corrected_positions.xyz')
rdf_plot_file: Output file for RDF plot (e.g., 'rdf_plot.png')
rdf_data_file: Output file for RDF data (e.g., 'rdf_data.csv')
atom_type1: First atom type for RDF calculation (e.g., 'Na')
atom_type2: Second atom type for RDF calculation (e.g., 'O')
dr: Bin width for RDF calculation (e.g., 0.1)
r_max: Maximum distance for RDF calculation (e.g., 10.0)

Usage:
Read XYZ Files: The script reads atomic positions from XYZ files matching the specified pattern.
Read Cell Lengths: The script reads cell lengths from the input file.
Apply Periodic Boundary Conditions: The script applies periodic boundary conditions to the atomic positions.
Calculate RDF: The script calculates the RDF for the specified atom types.
Save Results: The script saves the RDF plot as a PNG file and the RDF data as a CSV file.
