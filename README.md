**Code for computing lag analyses following procedure described in Raut et al. (2019) Neuroimage** 

**tdmx_template.m** -- The main script for performing lag analysis. Takes an input time series matrix and creates a pairwise time delay matrix. The script calls the following supporting functions (which should not require customization): \
**create_blocks.m** -- uses a temporal mask (e.g., indicating low-motion frames that should be included) of the time series to generate blocks of contiguous frames that meet the desired minimum block duration \
**lagged_cov.m** -- creates multivariate empirical cross-covariance functions (i.e., in units of the sampling interval) \
**parabolic_interp.m** -- performs parabolic interpolation to estimate the true peak of cross-covariance functions, and the corresponding time shift in seconds

**surrogate_TDE.m** is a function that can be used to generate surrogate pairs of time series that match the spectral characteristics of fMRI series, with zero-lag correlation and time delay, among other parameters, specified by the user. This function is useful for examining the dependency of time delay estimation error on data quantity and correlation magnitude (and likewise, estimating the error associated with an empirical time delay given a certain data quantity and correlation magnitude.
- f_alpha_gaussian.m pulled from (https://people.sc.fsu.edu/~jburkardt/m_src/cnoise/f_alpha_gaussian.m)

**Relevant citations:** \
Mitra, Anish, et al. "Lag structure in resting-state fMRI." Journal of Neurophysiology 111.11 (2014): 2374-2391. \
Raut, Ryan V., et al. "On time delay estimation and sampling error in resting-state fMRI." Neuroimage 194 (2019): 211-227.

=================

**Code for computing spectral analyses as in Raut et al. (2020) PNAS** 

**spectral_template.m** -- The main script for performing spectral analysis. Takes an input time series matrix and creates a Nx1 vector of "intrinsic timescales" for the N time series. The script calls the above functions "create_blocks.m" and "lagged_cov.m", as well as:      
**acf_hwhm.m** -- estimates intrinsic timescale for each time series, computed as half of the full-width-at-half-maximum of the autocorrelation function after spline fitting

**Relevant citation:** \
Raut, Ryan V., et al. "Hierarchical dynamics as a macroscopic organizing principle of the human brain." Proceedings of the National Academy of Sciences (in press).


**Please direct questions to ryanvraut [at] gmail.com**
