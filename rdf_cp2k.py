# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import glob
from multiprocessing import Pool
from tqdm import tqdm  # Import tqdm for progress bar

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

def read_xyz_files(file_pattern):
    files = glob.glob(file_pattern)
    positions = []
    atom_list = []
    num_frames = 0
    for file in tqdm(files, desc="Reading XYZ files"):
        with open(file, 'r') as f:
            lines = f.readlines()
            num_atoms = int(lines[0].strip())
            if not atom_list:
                atom_list = [line.split()[0] for line in lines[2:2+num_atoms]]
            for i in range(2, len(lines), num_atoms + 2):
                pos = np.array([list(map(float, lines[j].split()[1:4])) for j in range(i, i + num_atoms)])
                positions.append(pos)
                num_frames += 1
    print(f"Number of frames: {num_frames}")
    return atom_list, np.array(positions), num_frames

def read_cell_lengths(inp_file):
    with open(inp_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.strip().startswith('&CELL'):
                cell_lengths = []
                for cell_line in lines[lines.index(line)+1:]:
                    if cell_line.strip().startswith('&END CELL'):
                        break
                    # Extract only numerical values from the line
                    values = cell_line.split()
                    for value in values:
                        try:
                            cell_lengths.append(float(value))
                        except ValueError:
                            continue
                cell_matrix = np.array(cell_lengths[:9]).reshape(3, 3)  # Consider only the first 9 values and form a 3x3 matrix
                print("Cell matrix:")
                print(cell_matrix)
                return cell_matrix

def apply_periodic_boundary_conditions(positions, cell_matrix):
    inv_cell_matrix = np.linalg.inv(cell_matrix)
    for i in tqdm(range(positions.shape[0]), desc="Applying PBC"):
        for j in range(positions.shape[1]):
            pos = positions[i, j]
            pos = np.dot(pos, inv_cell_matrix)
            pos -= np.floor(pos + 0.5)  # Apply periodic boundary conditions
            positions[i, j] = np.dot(pos, cell_matrix)
    return positions

def write_corrected_xyz(file_name, atom_list, positions):
    with open(file_name, 'w') as f:
        num_atoms = len(atom_list)
        for i in tqdm(range(positions.shape[0]), desc="Writing corrected XYZ"):
            f.write(f"{num_atoms}\n")
            f.write(f"Frame {i+1}\n")
            for j in range(num_atoms):
                f.write(f"{atom_list[j]} {positions[i, j][0]} {positions[i, j][1]} {positions[i, j][2]}\n")

def calculate_distances(config, positions, atom_list, atom_type1, atom_type2, cell_matrix, r_max):
    num_atoms = positions.shape[1]
    distances = []
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            if atom_list[i] == atom_type1 and atom_list[j] == atom_type2:
                diff = positions[config, i] - positions[config, j]
                dist = np.linalg.norm(diff)
                if dist < r_max:
                    distances.append(dist)
    return distances

def calculate_rdf(positions, atom_list, atom_type1, atom_type2, cell_matrix, dr=0.1, r_max=10.0):
    num_configs = positions.shape[0]
    
    # Apply periodic boundary conditions to all positions
    positions = apply_periodic_boundary_conditions(positions, cell_matrix)
    
    # Write corrected positions to new XYZ file
    write_corrected_xyz(corrected_xyz_file, atom_list, positions)
    
    # Calculate distances with multiprocessing
    with Pool() as pool:
        distances = list(tqdm(pool.starmap(calculate_distances, [(config, positions, atom_list, atom_type1, atom_type2, cell_matrix, r_max) for config in range(num_configs)]), desc="Calculating distances"))
    
    # Flatten the list of distances
    distances = [dist for sublist in distances for dist in sublist]
    
    # Calculate RDF
    distances = np.array(distances)
    bins = np.arange(0, r_max + dr, dr)
    hist, bin_edges = np.histogram(distances, bins=bins, density=True)
    
    # Normalize RDF
    shell_volumes = 4/3 * np.pi * (bin_edges[1:]**3 - bin_edges[:-1]**3)
    number_density = len(distances) / (positions.shape[0] * np.prod(cell_matrix.diagonal()))
    rdf = hist / (shell_volumes * number_density)
    
    return bin_edges[1:], rdf

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
