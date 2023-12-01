# Morphometrics-for-Mathematica

![Alt text](https://github.com/pdpolly/Morphometrics-for-Mathematica/blob/main/GMMIconGitHub.jpg)

This Mathematica© add-on package performs common geometric morphometric functions. Includes Procrustes superimposition, thin-plate spline graphics, Centroid Size calculation, Procrustes distance calculation, Mantel tests, Ancestral Node reconstruction, Euclidean Distance Matrix Analysis (EDMA), and import function for TPS files.

This Mathematica© add-on package performs common geometric morphometric functions. Includes Procrustes superimposition, thin-plate spline graphics, Centroid Size calculation, Procrustes distance calculation, Mantel tests, Ancestral Node reconstruction, Euclidean Distance Matrix Analysis (EDMA), and import function for TPS files.

<b>Installation.</b> Download the ".m" format from <a href="https://github.com/pdpolly/Morphometrics-for-Mathematica/releases/latest">the latest release section</a>. Install using the "Install" item on the "File" menu. Once installed, you must load the package with the line <<PollyMorphometrics or whatever name you choose at installation.

<b>Changes in Version 12.x.</b> Upgrade of PrincipalCoordinates[] function to return principal coordinates scores from a data matrix based on Euclidean distance with option for Gower distance. Update in BreakOutlines[] function to space points at equal distances across the breakpoints in addition to within the outline segments. Fixed bugs that affected ShapeMANOVA[], ReconstructAncestorShapes[], TreeToMorphospace[], and PhylogeneticPrincipalComponentsOfShape[] in data sets with many taxa and landmarks. Added Mantel[] and MantelForMorphometrics[] functions to test similarity in covariance structure.  Added RVCoefficient[], CRCoefficient[], and AdamsCRTest[] functions for evaluating modularity and integration.  Also converted all "Modules" to "Blocks" for improved efficiency. 

This package is used in the course <a href="https://www.pollylab.org/courses/morphometrics">Geometric Morphometrics</a>.  Funding for development for this package has been provided by grants NSF EAR 1338298, the Robert R. Shrock fund at Indiana University, the Yale Institute for Biospheric Studies, and the Lilly Endowment through its support for the Indiana University Pervasive Technology Institute and the Indiana METACyt Initiative. 

<b>User manual:</b> https://github.com/pdpolly/Morphometrics-for-Mathematica/releases/download/v12.5/Guide.to.PollyMorphometrics.12.5.pdf 

<b>Cite as:</b> Polly, P.D. 2023. Geometric morphometrics for Mathematica. Version 12.5. https://github.com/pdpolly/Morphometrics-for-Mathematica [![DOI](https://zenodo.org/badge/513635593.svg)](https://zenodo.org/doi/10.5281/zenodo.10247111)
