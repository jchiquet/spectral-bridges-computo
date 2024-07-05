# Spectral Bridges algorithm

Spectral Bridges is a clustering algorithm that  builds upon the traditional k-means and spectral clustering frameworks by subdividing data into small Voronoi regions, which are subsequently merged according to a connectivity measure. Drawing inspiration from Support Vector Machine's margin concept, a non-parametric clustering approach is proposed, building an affinity margin between each pair of Voronoi regions. This approach is characterized by minimal hyperparameters and delineation of intricate, non-convex cluster structures.


The Spectral Bridge algorithm is implemented both in 
- Python (<https://github.com/flheight/Spectral-Bridges>) and
-  R <https://github.com/cambroise/spectral-bridges>).
