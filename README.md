# MSc Astrophysics Dissertation: Reproducibility Repository

This repository contains the complete source materials, data, and analytical scripts used in my MSc Astrophysics dissertation: **"GaussCluster: A 3D Gaussian Mixture Model for Open Cluster Membership."**

The objective of this research was to develop and validate a probabilistic 3D Gaussian Mixture Model (GaussCluster) to astrometrically isolate open clusters in heavily obscured environments, specifically targeting the complex Cygnus X star-forming region (Berkeley 87).

## Repository Structure

The repository is organized to ensure full reproducibility of the analyses and figures presented in the final manuscript:

* **Writeup/**: Contains the complete LaTeX source code (.tex files) for the dissertation manuscript, including all chapters and appendices, as well as the final compiled PDF.
* **GaussCluster_Final/**: The standalone software package developed for this dissertation. Contains the core GaussCluster.py algorithm, the validation script compare_clusters.py, and the raw and processed astrometric catalogues for Berkeley 87 and NGC 957.
* **Figures/**: All final visual assets, plots, and dashboards generated during the analysis and embedded in the manuscript.
* **CMD Isochrone/**: Python scripts and data used to perform PARSEC isochrone fitting and age derivation for the kinematically extracted cluster members.
* **King profile/**: Scripts used to fit empirical King Profile models to the spatial distributions of the clusters to validate their gravitational boundedness.
* **H alpha/**: Code used to analyze LCO 0.4m photometric data to place physical constraints on active circumstellar decretion disks among the B-type populations.

## Associated Software

The core GaussCluster software pipeline has also been published as a standalone, open-source Python tool for astronomical research. You can view the dedicated software repository here:
[https://github.com/raffey05-sudo/GaussCluster](https://github.com/raffey05-sudo/GaussCluster)

## Citation

If you utilize this data, methodology, or source code in your own research, please cite this dissertation and the associated GaussCluster pipeline.

