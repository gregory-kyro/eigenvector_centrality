import numpy as np
import matplotlib.pyplot as plt

# define function to compute adjacency matrix from generalized correlation coefficient matrix and distance matrix
def gen_adjacency_matrix(gcc_matx, distance_matx, lambd):
    adj_matx = np.zeros(gcc_matx.shape)
    for i in range(0, len(gcc_matx)):
        adj_matx[i, i] = 0.0
        distance_matx[i, i] = 0.0
        for j in range(0, len(gcc_matx)):
            adj_matx[i, j] = gcc_matx[i, j] * np.exp(-1.0 * distance_matx[i, j] / lambd)
            adj_matx[j, i] = adj_matx[i, j]
    return adj_matx

# create list of lambda values (arbitrary)
distance_labels = [100.0, 20.0, 15.0, 10.0, 7.5, 5.0, 2.5]

# compute adjacency matrix from generalized correlation coefficient matrix and distance matrix
corr_matx = generalized_correlation_coefficient
dist_matx = distancesWT

adjacency_matrix = [gen_adjacency_matrix(corr_matx, dist_matx, distance_labels[i]) for i in range(len(distance_labels))]

# define function to compute eigenvector centrality and eigenvalues from adjacency matrix
def calc_eigvec_centrality(adj_matx, get_eigval):
    eigval, eigvec = np.linalg.eigh(adj_matx)
    evec_centrality = np.zeros(len(eigval))
    for i in range(len(eigval)):
        evec_centrality[i] = np.sum(adj_matx[i, :] * np.abs(eigvec[:, -1]))/eigval[-1]
    if get_eigval == True:
        return evec_centrality, eigval
    else:
        return evec_centrality

# compute eigenvector centrality and eigenvalues of adjacency matrix
cent_results = [(calc_eigvec_centrality(adjacency_matrix[i], True)) for i in range(0, len(adjacency_matrix))]

# define eigenvector centralities
centralities = [cent_results[i][0] for i in range(0, len(cent_results))]

# define maximum eigenvalues of adjacency matrix
max_eigvals = [cent_results[i][-1] for i in range(0, len(cent_results))]

# define function to plot top ten eigenvalues
def plot_topten_eighvals(eigvals1, label):
    plt.bar(np.linspace(1, 10,10),eigvals1[-10:][::-1], width=0.4, label=label)

# plot top ten eigenvalues as a function of lambda
plt.figure(figsize=(12, 5))
counter = 0
for i in range(0, len(max_eigvals)):
    plot_topten_eighvals(max_eigvals[i], distance_labels[counter])
    counter += 1
plt.legend(ncol=3,loc=1, prop={'size': 12}, title='$\lambda$')
plt.ylabel('Eigenvalue')
plt.xlabel('Eigenvalue Index')

# define function to normalize eigenvector centralities
def normalize_centralities(centrality_vals):
    max_value = np.max(centrality_vals)
    min_value = np.min(centrality_vals)
    norm_centralities = np.zeros(centrality_vals.shape)
    for i in range(len(centrality_vals)):
        norm_centralities[i] = 2 * ((centrality_vals[i] - min_value)/(max_value - min_value)) - 1.0
    return norm_centralities

# normalize centralities
normed_centralities = [normalize_centralities(centralities[i]) for i in range(0, len(centralities))]

# plot normalized eigenvector centralities at various lambda values 
plt.figure(figsize=(16, 8), )
for i in range(len(normed_centralities)):
    plt.plot(normed_centralities[i], label=distance_labels[i], linewidth=1.2, ls='-')
plt.legend(prop={'size': 10}, title='$\lambda$')
plt.xlabel('Amino Acid Index')
plt.title(r'Normalized Eigenvector Centrality at Various $\lambda$ Values')
plt.ylabel('Eigenvector Centrality')
plt.tight_layout()
plt.hlines(0,0, len(normed_centralities[i]), color='k', linewidth=1, ls='--')
