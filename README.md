Copyright (c) 2019 Washington University  
Code written by: Anish Mitra, Aaron Tanenbaum & Ryan Raut

DISCLAIMER:
Washington University hereby grants to you a non-transferable, non-exclusive, royalty-free,
non-commercial, research license to use and copy the computer code that may be downloaded within 
this site (the "Software").  You agree to include this license and the above copyright notice in 
all copies of the Software.  The Software may not be distributed, shared, or transferred to any third party.  
This license does not grant any rights or licenses to any other patents, copyrights, or other forms of 
intellectual property owned or controlled by Washington University.  
If interested in obtaining a commercial license, please contact Washington University's Office of Technology 
Management (otm@dom.wustl.edu).

YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS PROVIDED "AS IS", 
WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES 
OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY THIRD-PARTY PATENT, 
COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE 
OR WASHINGTON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR CONSEQUENTIAL DAMAGES 
ARISING OUT OF OR IN ANY WAY CONNECTED WITH THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, 
WHETHER IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS ADVISED OF THE POSSIBILITY OF SUCH 
DAMAGES. 


NOTES:
The main script for lag analysis is tdmx_template.m, which uses a time series matrix to create a corresponding time delay matrix. The script calls the following supporting functions (which should not require customization):  
create_blocks.m -- uses a temporal mask (e.g., indicating low-motion frames that should be included) of the time series to generate blocks of contiguous frames that meet the desired minimum block duration  
lagged_cov.m -- creates multivariate empirical cross-covariance functions (i.e., in units of the sampling interval)
parabolic_interp.m -- performs parabolic interpolation to estimate the true peak of cross-covariance functions, and the corresponding time shift in seconds  

surrogate_TDE.m is a script that can be used to generate surrogate pairs of time series that match the spectral characteristics of fMRI series, with zero-lag correlation and time delay, among other parameters, specified by the user. This function is useful for examining the dependency of time delay estimation error on data quantity and correlation magnitude (and likewise, estimating the error associated with an empirical time delay given a certain data quantity and correlation magnitude.
- surrogate_TDE.m uses an external function, f_alpha_gaussian.m (https://people.sc.fsu.edu/~jburkardt/m_src/cnoise/f_alpha_gaussian.m) to simulate BOLD spectral content
