%% surrogate_TDE.m
%   This function 1) creates surrogate pairs of BOLD time series 
%   with modeled time delays, and 2) uses lagged cross-covariance and
%   parabolic interpolation to recover the modeled lag, outputting the bias,
%   variance, and root mean squared error associated with this estimate
%
%   To specify spectral content, this function calls another function, 
%   f_alpha_gaussian.m, which is freely available at:
%   https://people.sc.fsu.edu/~jburkardt/m_src/cnoise/f_alpha_gaussian.m
%
% INPUTS
%   alpha = exponent defining spectral content (i.e., 1/f^alpha)
%           use [] to specify default value of .7
%   r = built-in correlation between time series (prior to time shift)
%   tau_in = built-in time delay between time series (in seconds)
%           use [] to pull phi, the modeled time delay, from a uniform 
%           distribution [-max lag:max lag], where max lag = maxshift*tr)
%   duration = duration of scan (in minutes)
%   maxshift = number of TR shifts
%           use [] to specify default value, computed as round(4/TR + 1)
%   `           (i.e., enough shifts to resolve a lag of 4 seconds)
%   TR = stimulus interval (repetition time)
%   trials_per_sim = number of independent simulations of time series
%           (only useful when tau_in is set to [])
%           use [] to specify default value of 1
%   num_sims = number of independent simulations of the same time series
% 
% OUTPUTS
%   taus = true time delays
%   lags = time delay estimates
%   corrs = zero-lag correlation estimates
%   bias = bias of the measured delays compared to true delays
%   var = variance of time delay estimates
%   RMSE = root mean squared error of time delay estimates
%
% Example usage:
%   For 2000 simulations of 15-minute BOLD time series pairs with sampling
%   interval (repetition time, TR) = 2 seconds, zero-lag correlation
%   magnitude = .4 and time delay = .5 s, and time delay estimation using
%   cross-covariance and parabolic interpolation: 
% 
%   [taus,lags,bias,var,RMSE] = surrogate_TD([],.2,.5,15,[],2,[],2000)
%
% =============================
% Author: Ryan Raut
% Last modified: 03/14/2019
%
% Reference: Raut, R.V., Mitra, A., Snyder, A.Z., & Raichle, M.E. (2019)
% On time delay estimation and sampling error in resting-state fMRI.
% Neuroimage.
% =============================
%%

function [taus,lags,bias,var,RMSE] = surrogate_TDE(alpha,r,tau_in,duration,maxshift,TR,trials_per_sim,num_sims)

tic

taus = zeros(num_sims*trials_per_sim,1); % modeled time delays
lags = taus;    % time delay estimates

% loop on creating time series
for N = 1:num_sims    
    %% Create 1/f, band-passed time series
    
    % create 1/f^alpha time series
    nROIs = 2;    % if greater than 2, cannot specify lag
    Fs = 1/TR;
    nframes = round(duration*60*Fs);
    if isempty(alpha)
        alpha = .7;
    end
    X = zeros(nframes,nROIs);
    for i = 1:nROIs
        X(:,i) = f_alpha_gaussian(nframes,1,alpha);        
    end
    
    % bandpass (2nd order butterworth)
    hi_cut = .1;
    lo_cut = .005;
    [b,a] = butter(1,[lo_cut hi_cut]/(Fs/2));
    X = filtfilt(b,a,X);
    
    %% Set correlation structure
    
    % create correlation matrix
    C = repmat(r,nROIs);
    C(logical(eye(nROIs))) = ones;
    
    % orthogonalize time series
    X = zscore(X);
    [~,~,V] = svd(C);
    X = zscore(X*V);
     
    % apply correlations
    if r == 1
        X(:,2) = X(:,1);
    else
        % Cholesky factorization of C (C = F'F)
        F = chol(C);
        X = X*F;
    end
    
    nframes = round(duration*60*Fs);
    X = X(1:nframes,:);
    
    
    if isempty(trials_per_sim)
        trials_per_sim = 1;
    end
    % loop on lag induction/estimation
    for n = 1:trials_per_sim
        
        %% Create lagged time series
        
        if isempty(maxshift)
            maxshift = round(4/TR + 1);
        end
         
        % set time shift
        if isempty(tau_in)
            tau = randn(1)*maxshift;
        else
            tau = tau_in;
        end
                
        % apply shift in frequency domain
        y = fft(X(:,2));
        freqs = Fs * (0:1/length(X):1-1/length(X))';
        y = y.*exp(1i*-2*pi*tau*freqs);
        x = X(:,1);
        y = ifft(y,length(x),'symmetric');
        
        %% Estimate lags
        
        deltas = -maxshift:maxshift;
        Cov = zeros(1,1,numel(deltas));

        % compute cross-covariance function
        Cov = Cov + lagged_cov(x,y,maxshift);

        % normalize CCF
        for k = 1:numel(deltas)
            Cov(:,:,k) = Cov(:,:,k)/(length(x) - abs(deltas(k)));
        end
       
        % parabolic interpolation of CCF
        try
            [lags(N*n),~] = parabolic_interp(Cov,TR);
        catch %#ok<CTCH>
            lags(N*n) = nan;
        end                                

        taus(N*n) = tau;
    end
end

%% Compute statistics
bias = nanmean(lags - taus);
var = nanmean((lags-nanmean(lags)).^2);
RMSE = sqrt(nanmean((lags-taus).^2));

disp(['Bias = ' num2str(bias)]);
disp(['Variance = ' num2str(var)]);
disp(['RMSE = ' num2str(RMSE)]);

toc
end