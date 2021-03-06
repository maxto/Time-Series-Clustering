# Time-Series-Clustering

Time series clustering is to partition time series data into groups based on similarity or distance, so that time series in the same cluster are similar. The term "similar" is linked to the data type and the specific objective function we will apply.

Time series clustering belongs to the unsupervised learning methods and it can be divided in different parts:

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
