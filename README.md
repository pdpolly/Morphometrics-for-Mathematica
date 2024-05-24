# Morphometrics-for-Mathematica

![Alt text](https://github.com/pdpolly/Morphometrics-for-Mathematica/blob/main/GMMIconGitHub.jpg)

This Mathematica© add-on package performs common geometric morphometric functions. Includes Procrustes superimposition, thin-plate spline graphics, Centroid Size calculation, Procrustes distance calculation, Mantel tests, Ancestral Node reconstruction, Euclidean Distance Matrix Analysis (EDMA), and import function for TPS files.

<b>Installation.</b> Download the ".m" format from <a href="https://github.com/pdpolly/Morphometrics-for-Mathematica/releases/latest">the latest release section</a>. Install using the "Install" item on the "File" menu. Once installed, you must load the package with the line <<PollyMorphometrics or whatever name you choose at installation.

<b>Changes in Version 13.0.</b> Added MapShapeChangeOntoTree[] function to plot a tree colored by the per-branch rates of shape change; TreeToMorphospace3D[] function to plot three dimensions of a morphospace with phylogenetic tree; added Rotate3DPointCloud[] function to align meshes with Procrustes superimposed landmarks from those meshes.  fixed bug in BreakOutlines[] that dropped the 3rd coordinate in 3D semilandmark data; fixed bug in Procrustes that allowed some 3D shapes to become mirrored as part of post-superimposition alignment; removed confidence ellipses from TreetoMorphospace[]; fixed a bug in CommonOrientation[] that caused reference to be rotated to target instead of the opposite in some cases.

This package is used in the course <a href="https://www.pollylab.org/courses/morphometrics">Geometric Morphometrics</a>.  Funding for development for this package has been provided by grants NSF EAR 1338298, the Robert R. Shrock fund at Indiana University, the Yale Institute for Biospheric Studies, and the Lilly Endowment through its support for the Indiana University Pervasive Technology Institute and the Indiana METACyt Initiative. 

<b>User manual:</b> https://github.com/pdpolly/Morphometrics-for-Mathematica/releases/download/v13.0/Guide.to.PollyMorphometrics.13.0.pdf 

<b>Cite as:</b> Polly, P.D. 2024. Geometric morphometrics for Mathematica. Version 13.0 https://github.com/pdpolly/Morphometrics-for-Mathematica [![DOI](https://zenodo.org/badge/513635593.svg)](https://zenodo.org/doi/10.5281/zenodo.11288554)
