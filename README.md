# robust-surface-curvature
An unsupervised learning and robust surface curvature measure for surface complementarity - Abhijit Gupta, Arnab Mukherjee

## Description
1) The src directory contains jupyter notebook files for running model comparison - where we compare our model with coleman's model using analytical dataset comprising of perturbed points on spheres with known radius and curvature (curvature = 1/R)
2) The other jupyter notebook contains demo for Protein-inhibitor system, which showcases our curvature based surface complementarity function.
3) The utils directory has two programs - hypersphere.py, which has our core surface curvature algorithm, and read_msms.py, which is used for reading MSMS files.

## Requirements
1) For generating SES surface (solvent excluded molecular surface), please install MSMS program or use Chimera, which has built in utility of generating SES surface, saved as dot molecular surface.
2) Python libraries - 
   (i) NumPy
   (ii) Scipy
   (iii) matplotlib
   (iv) regex
   (v) unittest
   (vi) pandas
   (vii) Biopython


