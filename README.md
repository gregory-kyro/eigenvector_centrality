# protein_network_analyses -- Conceptual Overview

## DISPLACEMENT VECTORS - OVERVIEW

We record the x, y, and z coordinates for each $\alpha-$carbon of each residue for every frame in the trajectory in the numpy.ndarray ```coordinates```. We then compute the average x, y, and z coordinates for each atom. These are used to construct the atomic position fluctuation (displacement) vectors: 

$$ \mathbf{ r_{i}} = \mathbf{q_{i}} - \left \langle \mathbf{q_{i}} \right \rangle $$

or more explicitly: 

$$ \mathbf{r_{i}} = \left(q_{x} - \left \langle q_{x} \right \rangle \right ) \hat i + \left(q_{y} - \left \langle q_{y} \right \rangle \right )  \hat j + \left(q_{z} - \left \langle q_{z} \right \rangle \right )  \hat k $$

##  COVARIANCE MATRIX - OVERVIEW

We compute the distance matrix, which contains the average distances between residues $i, j$ over the course of the trajectory. This can be used for the calculation of Eigenvector Centrality.

We then compute the covariance matrix defined in terms of the atomic displacement vectors as follows: 

$$ \large{\Sigma_{i, j}} = \frac{\left \langle \mathbf{r_{i}} \cdot \mathbf{r_{j}} \right \rangle }{\left \langle \mathbf{r_{i}}^{2} \right \rangle \left \langle \mathbf{r_{j}}^{2} \right \rangle} $$

Finally, we compute the eigenvalues and eigenvectors of the covariance matrix, or the "Essential Modes" of the covariance matrix. This can be used for other analysis such as Principal Components Analysis (PCA).

## MUTUAL INFORMATION - OVERVIEW

We estimate the degree of dynamical correlation in terms of the Mutual Information, defined as: 

$$ \mathbf{I} \left[ \mathbf{x}_{i}, \mathbf{x}_{j} \right] = \sum_{k, m} P \left( k, m \right) \log \left[ \frac{P \left(k, m \right) }{P_{\mathbf{x}_{i}} (k)\  P_{\mathbf{x}_{j}} ( m ) } \right]$$ 

which can be expressed in terms of marginal and joint Shannon Entropies: 

$$ I \left( \mathbf{x}_{i}, \mathbf{x}_{j} \right) = H \left( \mathbf{x}_{i} \right) + H \left( \mathbf{x}_{j} \right) - H \left( \mathbf{x}_{i}, \mathbf{x}_{j} \right) $$

where

$$ H \left[ \mathbf{x}_{i} \right] = - \int d\mathbf{x}_{i}\ p \left[ \mathbf{x}_{i} \right] \ln \left( p \left[ \mathbf{x}_{i} \right] \right) $$

and

$$ H \left[ \mathbf{x}_{i}, \mathbf{x}_{j} \right] = - \iint d\mathbf{x}_{i}\ d\mathbf{x}_{j}\ p\left[\mathbf{x}_{i}, \mathbf{x}_{j} \right]\ \ln \left( p \left[ \mathbf{x}_{i}, \mathbf{x}_{j} \right] \right) $$

are the marginal and joint Shannon Entropies obtained as ensemble averages over atomic fluctuations  with marginal and joint probability distributions computed over thermal fluctuations sampled by MD simulations of the system at equilibrium.

**The Mutual Information is a measure of information shared between residues, that is, the information about one residue that can be gained from knowledge of the other residue.**

## GENERALIZED CORRELATION - OVERVIEW

The generalized correlation coefficient is calculated from the Mutual Information according to the following equation: 

$$ \tilde{ \mathbf{r_{\mathbf{MI}}}} \left[ \mathbf{x}_{i}, \mathbf{x}_{j} \right] = \left[ 1 - \exp \left( - \frac{2 I \left( \mathbf{x}_{i}, \mathbf{x}_{j} \right) }{3} \right) \right]^{1/2} $$ 

The generalized correlation coefficients, $\mathbf{r_{MI}}$, are used to construct our adjacency matrix, as well as a metric for the edge weights in our protein graph, quantifying the relationship between amino acid residues $i, j$.

They will also be used to compute the so-called "direct communication coefficients", which allow us to trace the paths of direct information transfer within the protein network in a later step.

Once the generalized correlation coefficient matrix has been computed, we compute its eigenvalues and eigenvectors. These can be used for the computation of Eigenvector Centrality (in conjunction with the average distances computed above).

## CENTRALITY - OVERVIEW
One of the simplest concepts when computing node-level measures is that of **centrality**, i.e. how central a node or edge is in the graph. 

**DEGREE CENTRALITY**

The **degree centrality** counts the number of edges adjacent to a node. The degree of node *i* is the number of existing edges $c_{ij}$ in other nodes *j* in a network with *n* nodes:
$$d_{ij} =\sum\limits_{j=1}^{n} e_{ij} ~ where: ~ i \neq j$$

**EIGENVECTOR CENTRALITY**

**Eigenvector centrality** is the weighted sum of the centralities of all nodes that are connected to it by an edge, as defined in the **adjacency matrix** $A_{i, j}$. The eigenvector centrality can be calculated according to the following equation:
    
$$ c_{i} = \varepsilon^{-1} \sum_{j=1}^{n} \mathbf{A}_{ij} c_{j}, \ \ \ \ \ \textrm{where } \mathbf{c} \textrm{ is the eigenvector associated with the largest eigenvalue } \varepsilon \textrm{ of } \mathbf{A} $$ 

For this calculation, the adjacency matrix is defined as follows: 
    
$$ \mathbf{A}_{ij} = \begin{cases} \\ 0 & \ \ \ \ \ \textrm{ if } i=j \\ \mathbf{r}_{{MI}} \left[x_{i}, x_{j}\right] \exp \bigg( - \frac{d_{ij}}{\lambda}\bigg) & \ \ \ \ \ \textrm{ if } i \neq j \\ \end{cases} $$ 
        
where $\lambda$ has been introduced as a length/locality factor to augment the locality of the iterations considered in the exponential damping parameter and $d_{ij}$ is the matrix of average distance of $\alpha$-Carbons $i$ and $j$.

The idea of a so-called *locality factor* presents a way of probing the effect of physical distances between residues, consolidating our focus to only those interactions that occur within this distance.

The centrality values can be normalized such that the nodes with the largest centralities will have values closer to $1.0$, while those with smaller centralities will have values closer to $-1.0$. This is done according to the following equation: 

$$ c_{i}^{\textrm{norm.}} = 2 \frac{c_{i} - c_{min}}{\left( c_{max} - c_{min} \right)} - 1 $$ 

**BETWEENNESS CENTRALITY**

The **betweenness centrality** of a node or edge in a network measures the extent to which it lies on shortest paths. A higher betweenness indicates that the node or edge lies on more shortest paths and hence has more influence on the flow of information in the graph. The *shortest path* (geodesic path) between two nodes in a graph is the path(s) with the least number of edges.

The geodesic betweenness $B_{n}(i)$ of a **node** in a weighted, undirected network is
$$B_{n}(i) =  \sum_{s,t \in G} \frac{ \Psi_{s,t}(i) }{\Psi_{s,t}}$$
where nodes $s,t,i$ are all different from each other.

* $\Psi_{s,t}$ denotes the number of shortest paths (geodesics) between nodes $s$ and $t$
* $\Psi_{s,t}(i)$ denotes the number of shortest paths (geodesics) between nodes $s$ and $t$ that pass through node $i$.
* The geodesic betweenness $B_n$ of a network is the mean of $B_n(i)$ over all nodes $i$
