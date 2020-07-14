**The main script for lag analysis is tdmx_template.m**, which uses a time series matrix to create a corresponding time delay matrix. The script calls the following supporting functions (which should not require customization):  
**create_blocks.m** -- uses a temporal mask (e.g., indicating low-motion frames that should be included) of the time series to generate blocks of contiguous frames that meet the desired minimum block duration  
**lagged_cov.m** -- creates multivariate empirical cross-covariance functions (i.e., in units of the sampling interval)
**parabolic_interp.m** -- performs parabolic interpolation to estimate the true peak of cross-covariance functions, and the corresponding time shift in seconds

**surrogate_TDE.m** is a function that can be used to generate surrogate pairs of time series that match the spectral characteristics of fMRI series, with zero-lag correlation and time delay, among other parameters, specified by the user. This function is useful for examining the dependency of time delay estimation error on data quantity and correlation magnitude (and likewise, estimating the error associated with an empirical time delay given a certain data quantity and correlation magnitude.
- surrogate_TDE.m uses an external function, f_alpha_gaussian.m (https://people.sc.fsu.edu/~jburkardt/m_src/cnoise/f_alpha_gaussian.m) to simulate BOLD spectral content

**The main script for spectral analysis is spectral_template.m**, which takes in a time series matrix and creates a Nx1 vector of "intrinsic timescales" for the N time series. The script calls the above functions "create_blocks.m" and "lagged_cov.m", as well as:      
**acf_hwhm.m** -- estimates intrinsic timescale for each time series, computed as half of the full-width-at-half-maximum of the autocorrelation function after spline fitting

**Questions may be directed to ryanvraut@gmail.com**
