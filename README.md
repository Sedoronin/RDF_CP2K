
# Radial Distribution Function (RDF) Calculation

This script calculates the Radial Distribution Function (RDF) for a given set of atomic positions from XYZ files. The script reads the atomic positions, applies periodic boundary conditions, calculates the distances between specified atom types, and generates the RDF. The results are saved as a plot and a CSV file.

## Requirements

- Python 3.x
- NumPy
- Matplotlib
- Glob
- Multiprocessing
- TQDM

You can install the required packages using pip:

pip install numpy matplotlib glob2 tqdm


## Parameters

- `file_pattern`: Pattern to match XYZ files (e.g., 'water-traj.xyz-pos-1.xyz')
- `inp_file`: Input file containing cell lengths (e.g., 'water.inp')
- `corrected_xyz_file`: Output file for corrected atomic positions (e.g., 'corrected_positions.xyz')
- `rdf_plot_file`: Output file for RDF plot (e.g., 'rdf_plot.png')
- `rdf_data_file`: Output file for RDF data (e.g., 'rdf_data.csv')
- `atom_type1`: First atom type for RDF calculation (e.g., 'Na')
- `atom_type2`: Second atom type for RDF calculation (e.g., 'O')
- `dr`: Bin width for RDF calculation (e.g., 0.1)
- `r_max`: Maximum distance for RDF calculation (e.g., 10.0)

## Usage

1. **Read XYZ Files**: The script reads atomic positions from XYZ files matching the specified pattern.

2. **Read Cell Lengths**: The script reads cell lengths from the input file.

3. **Apply Periodic Boundary Conditions**: The script applies periodic boundary conditions to the atomic positions.

4. **Calculate RDF**: The script calculates the RDF for the specified atom types.

5. **Save Results**: The script saves the RDF plot as a PNG file and the RDF data as a CSV file.

## Example

```python
# Set parameters
file_pattern = 'water-traj.xyz-pos-1.xyz'
inp_file = 'water.inp'
corrected_xyz_file = 'corrected_positions.xyz'
rdf_plot_file = 'rdf_plot.png'
rdf_data_file = 'rdf_data.csv'
atom_type1 = 'Na'
atom_type2 = 'O'
dr = 0.1
r_max = 10.0

# Read XYZ files
atom_list, positions, num_frames = read_xyz_files(file_pattern)

# Read cell lengths from .inp file
cell_matrix = read_cell_lengths(inp_file)

# Calculate RDF for two types of atoms (e.g., 'Na' and 'O')
r, rdf = calculate_rdf(positions, atom_list, atom_type1, atom_type2, cell_matrix, dr=dr, r_max=r_max)

# Save RDF plot to PNG file
plt.figure(figsize=(8, 6))
plt.plot(r, rdf, label='RDF')
plt.xlabel('Distance (r)')
plt.ylabel('g(r)')
plt.title('Radial Distribution Function')
plt.legend()
plt.savefig(rdf_plot_file)
plt.close()

# Save RDF data to CSV file
rdf_data = np.column_stack((r, rdf))
np.savetxt(rdf_data_file, rdf_data, delimiter=',', header='Distance, g(r)', comments='')

print(f"RDF plot saved as '{rdf_plot_file}'")
print(f"RDF data saved as '{rdf_data_file}'")

