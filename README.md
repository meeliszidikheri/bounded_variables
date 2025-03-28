The aim of this project is to implement bounded variable versions of the LETKF and GETKF algorithms. There are few steps to this process:
- Create a new aerosol model driven by L95 winds to evaluate the skill improvement enabled by the new algorithms and make necessary adjustments prior to testing with more realistic models.
- Create the capability of generating Inverse-Gamma statistical distributions for observation errors
- Refactor exsiting pertubed Gaussian LETKF and GETKF algorithms
- Implement new Gaussian-Inverse-Gamma (GIG) LETKF and GETKF algorithms

The GIG analysis update equations for EnKF/LETKF/GETKF are isomorphic to the corresponding Gaussian equations currently implemented in JEDI. Hence "only" need to change the calculation of the error covariance matrix, R, and the perturbed observations, yi (see attached slide for clarification).

[Bounded-variable EnKF.pptx](https://github.com/user-attachments/files/19496991/Bounded-variable.EnKF.pptx)
