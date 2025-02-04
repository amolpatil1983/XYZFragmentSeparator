import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
from ase import Atoms, io
from ase.io import read
from ase.neighborlist import neighbor_list
from ase.geometry.analysis import Analysis
import shutil

"""
Separates fragments from an xyz file into individual .xyz files, if any.
Uses numpy, ase, scipy and shutil.
Written By: Dr. Amol Patil, Chuo University, Japan 
"""

def get_ase_adjascency_matrix(inp_xyz_fl):
    ase_mol = read(inp_xyz_fl)
    symbols = ase_mol.get_chemical_symbols()
    natoms = len(symbols)
    ana = Analysis(ase_mol)
    bonds = ana.unique_bonds[0] #ana.all_bonds[0]
    adjascency_matrix = np.zeros((natoms, natoms), dtype=int)
    bonded_pairs = []
    for i in range(natoms):
        if len(bonds[i]) > 0:
            for j in range(len(bonds[i])):
                bonded_pairs.append([i, bonds[i][j]])
    for i in range(natoms):
        for j in range(i+1, natoms):
            initial_bond_order = adjascency_matrix[i, j]
            if [i, j] in bonded_pairs or (j, i) in bonded_pairs:
                new_bond_order = initial_bond_order + 1
                adjascency_matrix[i, j] = adjascency_matrix[j, i] = new_bond_order
    return adjascency_matrix, symbols, natoms
    
def group_fragments_by_connectivity(adjacency_matrix):
    sparse_matrix = csr_matrix(adjacency_matrix)
    num_components, labels = connected_components(sparse_matrix)

    fragments = [[] for _ in range(num_components)]
    for atom_index, component_label in enumerate(labels):
        fragments[component_label].append(atom_index)

    return fragments

def replace(xyzf):
    with open(xyzf,'r') as myxyz:
        with open('xyz_dummy.xyz','w') as mynewxyz:
            for line in myxyz.readlines():
                if '}' not in line:
                    mynewxyz.write(f'{line}')
    shutil.copy('xyz_dummy.xyz',xyzf)                
    return None

def extract_fragments_from_xyz(file_path, cutoff=1.2):
    basename = file_path.split('.')[0]
    replace(file_path)
    atoms = io.read(file_path)
    adjacency_matrix, symbols, natoms = get_ase_adjascency_matrix(file_path)
    fragments = group_fragments_by_connectivity(adjacency_matrix)

    if len(fragments) == 1:
        print("No separate fragments detected.")
    else:
        for idx, fragment_indices in enumerate(fragments):
            fragment_atoms = atoms[fragment_indices]
            io.write(f"{basename}_frag_{idx + 1}.xyz", fragment_atoms)
            print(f"Fragment {idx + 1} saved as {basename}_frag_{idx + 1}.xyz")

# Example usage.
extract_fragments_from_xyz("mymol.xyz")
