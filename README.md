# R codes for MVCM
Converted MATLAB to R to apply estimation procedure for varying coefficient functions and smoothing individual functions.


## R functions
`MVCM_test` is a test code to implement Zhu's (2012) method of a multivariate varying coefficient model (MVCM).
- `MVCM_read`: read raw data and generate, arc length, standardized design and response matrices, and related dimension parameters.
- `MVCM_lpks_wob`: read arc length, design and response matrices and generate optimal bandwidth for weighted least squares estimation.
- `MVCM_lpks_wb1`: read arc length, design and response matrices, and optimal bandwidth and generate the estimated coefficient functions, their first derivatives, and fitted responses using weighted least squares estimation.
- `MVCM_sif`: read arc length and the residuals of response and generate the estimations of η, ε, and Ση(s,t).


## References
1. Zhu et al.,2012. https://doi.org/10.1214/12-AOS1045
2. MATLAB codes for MVCM. https://www.nitrc.org/projects/fadtts/
3. MATLAB/R Reference. https://umaine.edu/mathematics/david-hiebeler/computing-software/matlab-r-reference/
