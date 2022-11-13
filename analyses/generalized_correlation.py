import matplotlib.pyplot as plt
imort numpy as np

# define function to compute generalized correlation coefficient matrix from mutual information matrix
def compute_generalized_correlation_coefficient(twobody_corr_matrix, natoms, id=None, out_path=None, img_path=None, save=True, plot=True):
    for j in range(0, natoms):
        for i in range(0, j):
            twobody_corr_matrix[j][i] = np.exp(-(2. / 3.) * twobody_corr_matrix[j][i])
            twobody_corr_matrix[j][i] = 1.0 - twobody_corr_matrix[j][i]

            twobody_corr_matrix[j][i] = np.emath.sqrt(twobody_corr_matrix[j][i])
            twobody_corr_matrix[i][j] = twobody_corr_matrix[j][i]
        twobody_corr_matrix[j][j] = 1.0

    if plot:
        fig, ax1 = plt.subplots(figsize=(5, 5))
        pos = ax1.imshow(twobody_corr_matrix, cmap='viridis', origin='lower', interpolation='none')
        fig.colorbar(pos, ax=ax1)
        plt.savefig(img_path + 'generalized_correlation_coefficient_' + str(id) + '.pdf')
        plt.close()

    if save:
        dump_matrix(twobody_corr_matrix, 'generalized_correlation_coefficient_' + str(id), out_path)
        dump_distance_matrix(twobody_corr_matrix, out_path, 'generalized_correlation_coefficient_matrix_'+ str(id))

    return twobody_corr_matrix

# define function to dump matrix
def dump_matrix(matrix, name_of_file, path):
    n = len(matrix)
    f = open(path + name_of_file + ".txt", "w+")
    for i in range(0, n):
        for j in range(0, n):
            f.write('%10s %10s %10.6f\n' % (i, j, matrix[i][j]))
    f.close()
    return

# compute generalized correlation coefficient matrix
generalized_correlation_coefficient = compute_generalized_correlation_coefficient(linearized_mutual_information, number_of_atoms, id='filename', out_path=output_path, img_path=plot_path, save=True, plot=True)

# compute eigenvalues and eigenvectors of generalized correlation coefficient matrix
gcc_eigvals, gcc_eigvec = np.linalg.eigh(generalized_correlation_coefficient)

# define function to check for symmetry
def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)

# check if generalized correlation coefficient matrix is symmetrical
print('Generalized correlation coefficient matrix is symmetrical: ',check_symmetric(generalized_correlation_coefficient))

# print ratio between largest and second-largest eigenvalues
print("lambda_1 / lambda_2 : ", gcc_eigvals[-1]/gcc_eigvals[-2])

# print minimum eigenvalue of matrix
print("Min. of the correlation Matrix: ", np.min(generalized_correlation_coefficient))

# print maximum eigenvalue of matrix
print("Max. of the correlation Matrix: ", np.max(generalized_correlation_coefficient))

# plot top ten eigenvalues
plt.figure(figsize=(12, 5))
plt.bar(np.linspace(1, 10, 10)-.3, gcc_eigvals[::-1][:10], width=0.2, label='WT', color='red')
plt.xticks(ticks=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
plt.xlabel('Eigenvalue Index', size=20)
plt.ylabel('Eigenvalue', size=20)
plt.xticks(size=20)
plt.yticks(size=20)
