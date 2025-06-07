# Time-Series-Clustering

Time series clustering involves partitioning a set of time-dependent data series into groups (clusters) such that series in the same cluster are more similar to each other than to those in other clusters. It is an unsupervised learning task, meaning the clustering is done without predefined labels. The notion of "similarity" here depends on how we measure distance or similarity between time series and the objective of the analysis. Key components of a time-series clustering procedure include: (1) the choice of distance or similarity measure between time series, (2) the clustering algorithm or method used to form clusters, and (3) the evaluation of the clustering result. In literature, time-series clustering methods are often categorized by how they handle the data, for example as **shape-based**, **feature-based**, or **model-based** approaches (and we will also discuss some **structure-based** measures and other modern techniques). Below we outline the major aspects of time-series clustering:

---

## Distance and Similarity Measures

A fundamental step is defining a distance or dissimilarity measure between time series. This measure determines what it means for two series to be "similar." Different types of distance measures have been developed for time-series, which we group into shape-based (direct comparison of the time-series sequences), feature-based (comparison of derived attributes), and structure-based (comparison of underlying model or information content) categories.

### Shape-Based Measures

Shape-based distances compare time series *point-by-point* based on their raw shape (the sequence of values). These methods typically require the series to be aligned in time (or they employ mechanisms to align them) and focus on matching the observed patterns. We further distinguish **lock-step** distances (which compare values at the same time indices) and **elastic** distances (which allow some stretching or shifting of the time axis to find a better match).

#### Lock-step distances

| Distance          | Description                                                                    | R Packages/Functions                  |
| ----------------- | ------------------------------------------------------------------------------ | ------------------------------------- |
| Euclidean         | Standard L2 norm; measures direct, pointwise similarity (sensitive to shifts). | `dist()`, `TSdist::EuclideanDistance` |
| Manhattan         | L1 norm; sum of absolute differences (less sensitive to outliers).             | `dist()`, `TSdist::ManhattanDistance` |
| Minkowski         | Generalized Lp norm; allows tuning emphasis.                                   | `dist(method="minkowski")`            |
| Chebyshev         | Maximum absolute difference.                                                   | `TSdist::ChebyshevDistance`           |
| Canberra          | Weighted L1; emphasizes small values.                                          | `TSdist::CanberraDistance`            |
| Binary            | For binary/categorical series (Hamming distance).                              | `dist(method="binary")`               |
| Correlation       | Uses Pearson/Spearman/Kendall correlation (distance = 1 - correlation).        | `TSdist::correlationBasedDistance`    |
| Cross-correlation | Max correlation at any lag; captures phase-shifted similarity.                 | Custom, or `TSclust::diss.COR`        |
| Dissim            | Area between curves; handles uneven sampling.                                  | `TSdist::DissimDistance`              |

#### Elastic measures

| Measure                           | Description                                          | R Packages/Functions                            |
| --------------------------------- | ---------------------------------------------------- | ----------------------------------------------- |
| Dynamic Time Warping (DTW)        | Aligns series by warping time axis to minimize cost. | `dtw::dtw()`, `TSdist::DTWDistance`, `dtwclust` |
| Frechet distance                  | Minimum leash length between curves.                 | `TSdist::FrechetDistance`                       |
| Longest Common Subsequence (LCSS) | Finds longest subsequence allowing mismatches/gaps.  | `TSdist::LCSSDistance`                          |
| Shape-Based Distance (SBD)        | Max cross-correlation after shifting.                | `dtwclust::SBD`                                 |
| Shapelets/U-Shapelets             | Uses local subsequence patterns for clustering.      | Custom/Research code                            |
| Global Alignment Kernel (GAK)     | Kernel-based soft alignment.                         | `kertime` (less common in R)                    |

### Feature-Based Measures

Transform time series into feature vectors and apply standard clustering.

| Feature type        | Description                                       | R Packages/Functions                                   |
| ------------------- | ------------------------------------------------- | ------------------------------------------------------ |
| ACF/PACF            | Autocorrelation/Partial autocorrelation features. | `acf()`, `TSclust::diss.ACF`                           |
| Fourier/Periodogram | Frequency domain/spectral features.               | `TSA`, `TSdist::PerDistance`, `TSdist::IntPerDistance` |
| Wavelet             | Multi-scale, localized frequency info.            | `WaveletComp`, `waveslim`                              |
| SAX                 | Symbolic representation for coarse shapes.        | `jmotif`, `TSdist::MindistSAXDistance`                 |
| Trend/Stats         | Mean, variance, skewness, etc.                    | `feasts`, `tsfeatures`                                 |

### Structure-Based Measures

Compare underlying models or information content.

| Measure                    | Description                                              | R Packages/Functions                           |
| -------------------------- | -------------------------------------------------------- | ---------------------------------------------- |
| Model-based                | Distance between fitted models (e.g., AR coefficients).  | `TSclust::diss.AR.PIC`, `TSclust::diss.AR.MAH` |
| Cepstral                   | Compares cepstral/LPC coefficients.                      | `TSclust::diss.AR.LPC.CEPS`                    |
| Compression-based (NCD)    | Uses compression size (Kolmogorov complexity).           | `TSdist::NCDDistance`                          |
| Complexity-Invariant (CID) | Adjusts for series complexity (number of local changes). | `TSdist::CIDDistance`                          |
| Permutation Distribution   | Based on entropy of ordinal patterns.                    | `pdc`, `TSdist::PDCDistance`                   |

---

## Clustering Methods

### Partitioning

| Method                            | Description                                               | R Packages/Functions             |
| --------------------------------- | --------------------------------------------------------- | -------------------------------- |
| K-Means                           | Iterative assignment to minimize within-cluster variance. | `kmeans()`, `dtwclust::tsclust`  |
| K-Medoids (PAM/CLARA)             | Uses medoids as centers, works with arbitrary distances.  | `cluster::pam`, `cluster::clara` |
| Self-Organizing Map (SOM)         | Neural net for topological mapping.                       | `kohonen`                        |
| Expectation Maximization (EM/GMM) | Soft clustering via mixture models.                       | `mclust`, `mixtools`             |
| K-Shape                           | Specialized for time series, shape-based centroids.       | `dtwclust`                       |
| TADPole                           | Efficient DTW-based clustering.                           | `dtwclust`                       |
| Affinity Propagation (AP)         | Finds exemplars based on similarity matrix.               | `apcluster`                      |

### Hierarchical Clustering

| Method        | Description                   | R Packages/Functions            |
| ------------- | ----------------------------- | ------------------------------- |
| Agglomerative | Bottom-up merging by linkage. | `hclust()`, `dtwclust::tsclust` |
| Divisive      | Top-down splitting.           | `cluster::diana`                |

### Fuzzy Clustering

| Method              | Description                                                 | R Packages/Functions              |
| ------------------- | ----------------------------------------------------------- | --------------------------------- |
| Fuzzy C-Means (FCM) | Each series belongs to clusters with degrees of membership. | `e1071::cmeans`, `cluster::fanny` |

### Density-Based Clustering

| Method                        | Description                                         | R Packages/Functions |
| ----------------------------- | --------------------------------------------------- | -------------------- |
| DBSCAN                        | Density-based spatial clustering, detects outliers. | `dbscan::dbscan`     |
| Shared Nearest Neighbor (SNN) | Clusters based on shared neighbors.                 | `dbscan::sNNclust`   |

### Model-Based/Other

| Method                       | Description                               | R Packages/Functions   |
| ---------------------------- | ----------------------------------------- | ---------------------- |
| Finite Mixture Model         | Mixture of generative models.             | `mclust` (on features) |
| Deep Clustering              | Neural net feature learning + clustering. | `keras`, custom        |
| Robust Continuous Clustering | Graph-based continuous optimization.      | Custom/research code   |

---

## Cluster Validity Indices

### Internal Indices

| Index                                       | Description                                                          | R Packages/Functions                |
| ------------------------------------------- | -------------------------------------------------------------------- | ----------------------------------- |
| Silhouette                                  | Measures cohesion vs separation for each point; avg. near 1 is good. | `cluster::silhouette`, `factoextra` |
| Dunn                                        | Ratio of min inter-cluster distance to max intra-cluster diameter.   | `clusterCrit::intCriteria`          |
| Calinski-Harabasz                           | Ratio of between/within cluster variance.                            | `clusterCrit::intCriteria`          |
| Davies-Bouldin                              | Avg. similarity to closest cluster (lower is better).                | `clusterCrit::intCriteria`          |
| Ball-Hall, Ray-Turi, Xie-Beni, S\_Dbw, etc. | Other internal indices for specific criteria.                        | `clusterCrit`                       |

### External Indices

| Index                                   | Description                                     | R Packages/Functions       |
| --------------------------------------- | ----------------------------------------------- | -------------------------- |
| Rand, Adjusted Rand (ARI)               | Pair-counting index, corrected for chance.      | `clusterCrit::extCriteria` |
| Jaccard, Dice                           | Pairwise overlap between cluster/class labels.  | `clusterCrit::extCriteria` |
| Fowlkes-Mallows                         | Geometric mean of precision & recall for pairs. | `clusterCrit::extCriteria` |
| Hubert's Gamma, Phi, Sokal-Sneath, etc. | Other pair-counting/statistical indices.        | `clusterCrit::extCriteria` |

---

## Example: Clustering Synthetic Time Series

```r
# Load libraries
library(dtw)
library(proxy)
library(cluster)
library(clusterCrit)

# Create 4 example time series
set.seed(123)
time <- seq(0, 2*pi, length.out=100)
series1 <- sin(time)
series2 <- c(rep(0,10), sin(time[1:90]))
series3 <- rep(0,100); series3[50] <- 1
series4 <- rep(0,100); series4[70] <- 1
series_matrix <- rbind(series1, series2, series3, series4)
rownames(series_matrix) <- paste0("Series", 1:4)

# Compute distances
euclid_dist <- dist(series_matrix)
dtw_dist <- proxy::dist(series_matrix, method=function(x, y) dtw(x, y)$distance)

# Hierarchical clustering
hc_euclid <- hclust(euclid_dist, method="complete")
hc_dtw <- hclust(dtw_dist, method="complete")

# Cut to 2 clusters
cutree(hc_euclid, k=2)
cutree(hc_dtw, k=2)

# Silhouette index
sil_euclid <- silhouette(cutree(hc_euclid,2), euclid_dist)
sil_dtw <- silhouette(cutree(hc_dtw,2), dtw_dist)
summary(sil_euclid)$avg.width
summary(sil_dtw)$avg.width
```

---

## References & Further Reading

* Montero, P., & Vilar, J. A. (2014). TSclust: An R package for time series clustering. *Journal of Statistical Software*, 62(1), 1-43.
* Mori, U., Mendiburu, A., & Lozano, J. A. (2016). TSdist: A Package for Time Series Distance Measures in R. *The R Journal*, 8(2), 451-459.
* Paparrizos, J., & Gravano, L. (2015). k-Shape: Efficient and Accurate Clustering of Time Series. *SIGMOD*.
* dtwclust R package documentation and vignettes.
* Aghabozorgi, S., Shirkhorshidi, A. S., & Wah, T. Y. (2015). Time-series clustering–A decade review. *Information Systems*, 53, 16-38.
* clusterCrit R package documentation.

---

**Note:** For the latest algorithms and deep learning approaches, see recent research and Python packages such as `tslearn`, `sktime`, or deep clustering repositories.

General schema:

- Distance and similarity/dissimilarity
  - Shape-based
    - Lock-step
      - Euclidean
      - Manhattan
      - Minkowski
      - Mahalanobis
      - Maximum
      - Canberra
      - Binary
      - Correlation
        - Pearson
        - Spearman
        - Kendall tau
      - Cross-Correlation
      - Dissim
    - Elastic measure
      - Dynamic Time Warping (DTW)
      - Frechet distance
      - Longest Common Subsequence (LCSS)
      - Shape-based
      - U-Shapelets
      - Kernel
      - Global Alignement Kernel (GAK)
  - Features-based
    - Partial/Autocorrelation
    - Fourier
    - Wavelet
    - Periodogram
    - SAX
    - Spectral Density
  - Structure-based
    - Model
      - Cepstral
    - Compression
      - Compression
      - Complexity invariant
      - Permutation distribution based distance
      
    
- clustering methods
  - Partitioning
    - K-Means
    - Self-organizing map (SOM)
    - K-nearest Neighbour (KNN)
    - Expectation Maximization (EM)
    - K-Medoids (PAM)
    - Clustering Large Applications (CLARA)
    - K-Shape
    - TADPole clustering (TADP)
    - Affinity Propagation (AP)
  - Hierarchical clustering
    - Agglomerative
     - Single-linkeage
     - Average
     - Complete
     - Ward
     - McQuitty
     - Median
     - Centroid
    - Divisive
  - Fuzzy clustering
    - Fuzzy C-means (FCM)
  - Density-based clustering
    - Density-Based Spatial Clustering and Application with Noise (DBSCAN)
    - Shared Nearest Neighbor (SNN)
  - Model-based clustering
    - Finite Mixture Model
  - Grid-based clustering
    - STatistical INformation Grid-based method (STING)
    - Clustering In QUEst (CLIQUE)
  - Deep clustering
    - Deep Continuous Clustering (DOC)
  - continuous objective optimization
    - Robust Continuous Clustering
 
 - Cluster Validity Index
    - Internal
      - Ball-Hall index
      - Banfeld-Raftery index
      - C index
      - Calinski-Harabasz index
      - Davies-Bouldin index
      - Det Ratio index
      - Dunn index
      - Baker-Hubert Gamma index
      - GDI index
      - G plus index
      - Ksq DetW index
      - Log Det Ratio index
      - Log SS Ratio index
      - McClain-Rao index
      - PBM index
      - Point-Biserial
      - Ratkowsky-Lance index
      - Ray-Turi index
      - Scott-Symons index
      - SD index
      - S Dbw index
      - Silhouette index
      - Tau index
      - Trace W index
      - Trace WiB index
      - Wemmert-Gan¸carski index
      - Xie-Beni index
    - External
      - Czekanowski-Dice index
      - Folkes-Mallows index
      - Hubert Γ indexˆ
      - Jaccard index
      - Kulczynski index
      - McNemar index
      - Phi index
      - Rand index
      - Rogers-Tanimoto index
      - Russel-Rao index
      - Sokal-Sneath indices
    
### R packages for Time Series:
  
  - [TSclust](https://cran.r-project.org/web/packages/TSclust/index.html) A set of measures of dissimilarity between time series to perform time series clustering. Metrics based on raw data, on generating models and on the forecast behavior are implemented. Some additional utilities related to time series clustering are also provided, such as clustering algorithms and cluster evaluation metrics.
 
  - [dtwclust](https://cran.r-project.org/web/packages/dtwclust/index.html) Time series clustering along with optimized techniques related to the Dynamic Time Warping distance and its corresponding lower bounds. Implementations of partitional, hierarchical, fuzzy, k-Shape and TADPole clustering are available. Functionality can be easily extended with custom distance measures and centroid definitions. Implementations of DTW barycenter averaging, a distance based on global alignment kernels, and the soft-DTW distance and centroid routines are also provided. All included distance functions have custom loops optimized for the calculation of cross-distance matrices, including parallelization support. Several cluster validity indices are included.
  
  - [TSrepr](https://cran.r-project.org/web/packages/TSrepr/index.html) Methods for representations (i.e. dimensionality reduction, preprocessing, feature extraction) of time series to help more accurate and effective time series data mining. Non-data adaptive, data adaptive, model-based and data dictated (clipped) representation methods are implemented. Also min-max and z-score normalisations, and forecasting accuracy measures are implemented.
  
  - [LPWC](https://cran.r-project.org/web/packages/LPWC/index.html) Computes a time series distance measure for clustering based on weighted correlation and introduction of lags. The lags capture delayed responses in a time series dataset. The timepoints must be specified.
  
  - [clValid](https://cran.r-project.org/web/packages/clValid/index.html) Statistical and biological validation of clustering results.
  
  - [NbClust](https://cran.r-project.org/web/packages/NbClust/index.html) It provides 30 indexes for determining the optimal number of clusters in a data set and offers the best clustering scheme from different results to the user.
  
  - [mclust](https://cran.r-project.org/web/packages/mclust/index.html) Gaussian finite mixture models fitted via EM algorithm for model-based clustering, classification, and density estimation, including Bayesian regularization, dimension reduction for visualisation, and resampling-based inference.
  
  - [clusterCrit](https://cran.r-project.org/web/packages/clusterCrit/index.html) Compute clustering validation indices.


