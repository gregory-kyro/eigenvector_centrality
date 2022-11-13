import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import MDAnalysis as mda
from scipy.spatial import distance_matrix

# create universe 
top = '\path\to\top.pdb'
traj = '\path\to\traj.dcd'
u_prot = mda.Universe(top,traj)

# select alpha carbon atoms
ca_prot = u_prot.select_atoms('name CA')

# contruct an array of alpha carbon coordinates
number_of_frames = 2002
number_of_atoms = 122

coordinates = np.ndarray((number_of_frames, number_of_atoms, 3))
counter = 0
for frame in u_prot.trajectory:
    coordinates[counter] = ca_prot.positions
    counter += 1
    
# define function to reshape and average coordinate array
def coordinate_reshape_and_average(coordinates_mat, nframes, natoms):
    coords_reshaped = np.zeros((nframes + 1, 3 * natoms))
    for i in range(0, nframes):
        counter = 0
        for j in range(0, natoms):
            coords_reshaped[i][counter] = coordinates_mat[i, j, 0]
            coords_reshaped[i][counter + 1] = coordinates_mat[i, j, 1]
            coords_reshaped[i][counter + 2] = coordinates_mat[i, j, 2]
            counter += 3

    average = np.mean(coords_reshaped, axis=0)
    coords_reshaped[nframes] = average

    return coords_reshaped

# reshape and average coordinate array
coords = coordinate_reshape_and_average(coordinates, number_of_frames, number_of_atoms)

# define function to construct displacement matrix
def compute_displacements(coords_reshaped, natoms, nframes):
    displacements = np.zeros((np.shape(coords_reshaped)[0]-1, np.shape(coords_reshaped)[1]))
    for i in range(0, 3 * natoms):
        for j in range(0, nframes):
            displacements[j][i] = coords_reshaped[j][i] - coords_reshaped[-1][i]
    return displacements

# construct displacement matrix
displacement_matrix = compute_displacements(coords, number_of_atoms, number_of_frames)

# define function to compute distance matrix
def compute_distance_matrix(coords_reshaped, natoms, nframes, id=None, out_path=None, DEBUG=False, RMSF=False):
    coords_cop = coords_reshaped[nframes].view()
    coords_cop = np.reshape(coords_cop, [natoms, 3])
    dist_matrix = distance_matrix(coords_cop, coords_cop)
    jj = 0
    ii = 0
    if DEBUG:
        file_name = 'distances_K810A' + str(id)
        dump_distance_matrix(dist_matrix, out_path, file_name)

        if RMSF:
            ii = 0
            with open(out_path + str('RMSF_') + str(id) + '.txt', 'w+') as f:
                for i in range(natoms):
                    rmsf = 0.0
                    rmtemp = 0.0
                    temp0 = coords_reshaped[nframes, ii + 0]
                    temp1 = coords_reshaped[nframes, ii + 1]
                    temp2 = coords_reshaped[nframes, ii + 2]

                    for j in range(nframes):
                        rmtemp = (coords_reshaped[j][ii + 0] - temp0) ** 2 + \
                                 (coords_reshaped[j][ii + 1] - temp1) ** 2 + \
                                 (coords_reshaped[j][ii + 2] - temp2) ** 2
                        rmsf = rmsf + rmtemp / nframes

                    ii = ii + 3
                    rmsf = np.sqrt(rmsf)
                    rmsf = np.asarray(rmsf)
                    f.write('%10.6f\n' % rmsf)
                    rmsf = 0.

                temp0 = 0
                temp1 = 0
                temp2 = 0
                ii = 0
            f.close()
    return dist_matrix

# define function to dump distance matrix
def dump_distance_matrix(distance_matrix, path, filename):
    n = len(distance_matrix)
    f = open(path + filename + ".txt", "w+")
    for i in range(0, n):
        for j in range(0, n):
            f.write("%10.6f" % distance_matrix[i, j])
        f.write("\n")
    f.close()

output_path='/path/to/save'
# compute distance matrix
distances = compute_distance_matrix(coords, number_of_atoms, number_of_frames, 'filename', output_path, RMSF=True)
np.savetxt(output_path + "distance_matrix.txt", distances)

# define function to compute covariance matrix
def evaluate_covariance_matrix(displacement_matrix, id=None, out_path=None, img_path=None, DEBUG=False, plot=False):
    covar_mat = np.cov(displacement_matrix, rowvar=False, bias=True)
    return covar_mat

plot_path='/path/to/save/plot'
# compute covariance matrix
covariance_matrix = evaluate_covariance_matrix(displacement_matrix, 'filename', output_path, plot_path, plot=True)

# print the shape of the WT covariance matrix (3N_residues x 3N_residues)
print('The shape of the WT covariance matrix is: ', covariance_matrix.shape)

# plot covariance matrix
plt.figure(figsize=(10,10))
plt.imshow(covariance_matrix, origin='lower')
plt.xlabel('Node-Coordinate Index')
plt.ylabel('Node-Coordinate Index')
plt.colorbar()

# define function to compute eigenvalues and eigenvectors of covariance matrix (diagonalize the matrix)
def diagonalize_square_matrix(matrix_to_diagonalize, id=None, out_path=None, file_name=None, DEBUG=False):
    eigenvalues, eigenvectors = LA.eig(matrix_to_diagonalize)
    return eigenvalues, eigenvectors

# define function to dump eigenvalues and eigenvectors (essential modes)
def dump_essentials(eigenvalues, eigenvectors, matrix_name, path):
    f1 = open(path + str(matrix_name) + "_eigenvec.txt", "w+")
    f2 = open(path + str(matrix_name) + "_eigenval.txt", "w+")
    n = len(eigenvalues)
    for i in range(0, n):
        f2.write('%10s %10.6f\n' % (i, eigenvalues[i].real))
        for j in range(0, len(eigenvectors)):
            f1.write('%10.6f\n' % eigenvectors[i, j].real)
    f1.close()
    f2.close()
    return

# compute eigenvalues and eigenvectors of covariance matrix
eigvals, eigvec = diagonalize_square_matrix(covariance_matrix, 'fileame', output_path)
