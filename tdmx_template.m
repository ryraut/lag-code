%% Main lag script used to compute time delay matrix and lag projection
%
%% General information for a few variables and their expected file types
%
%  subjects:  various compatible file formats; read as string array
%  BOLD:    .mat file; Time x Space (where the spatial dimension corresponds to a vector of spatial units...in other words, BOLD is a 2D matrix, where spatial units are vectorized)
%  good:    .mat file; Vector with length = number of spatial units; used to mask out bad spatial units
%  format:  various compatible file formats (e.g., text file); vector with length = number of time points; used to mask out bad time points
%%%%%%%%
%
%% Setup
% Set parameters
num_nodes = 264;    % number of time series (i.e., spatial dimension)
outdir = '';    % set directory for saving out images here
mkdir(outdir)
lag_lim = 4;    % lag limit (in seconds)
lags = -3:3;    % range of TR shifts; max(lags) = round(lag_lim/tr + 1)

% Specify data parameters
subjects = importdata('');
tr = 2; % sampling interval in seconds
motion_thresh = .2;    % important: must match motion criteria used during preproc

min_block_durn = (max(lags)+1)*tr;   % min. block duration (in seconds)

%% Loop over subjects
% initialize group matrices for running sum
grp_lags = single(nan(num_nodes));    % peak lags
grp_ZL = grp_lags;      % zero-lag correlation
grp_peak = grp_lags;    % peak correlation

grp_lags_nans = single(zeros(num_nodes));
grp_ZL_nans = grp_lags_nans;
grp_peak_nans = grp_lags_nans;

for s = 1:numel(subjects)
    tic
    subj = subjects{s};
    disp(['Processing ' subj]);
    
    % initialize subject matrices
    subj_lags = single(nan(num_nodes)); % peak lags
    subj_ZL = subj_lags;   % zero-lag correlation
    subj_peak = subj_lags; % peak correlation (correlation at optimal lag)
    
    
    BOLD = importdata(''); % read in time series matrix
    good = importdata(''); % read in spatial mask if desired; MAKE SURE THIS IS A LOGICAL!
    
    % read in temporal mask/motion time series (e.g., FD or DVARS); MAKE SURE THIS IS A LOGICAL!
    format = dlmread('') <= motion_thresh;
    
    % ignore pre-steady-state frames
    format(1:2) = false; % ignore first X frames
    
    FORMAT = create_blocks(format,min_block_durn,tr);
    
    %% Do the lagged correlation/covariance computation of TD matrices
    Cov = zeros([sum(good) sum(good) numel(lags)]);
    nblocks = numel(FORMAT);
    nframes = 0;
    
    % De-mean time series
    run_mean = nanmean(BOLD(format,:),1);
    BOLD = bsxfun(@minus,BOLD,run_mean);
    
    % Loop over blocks of contiguous frames
    for j = 1:numel(FORMAT)
        nframes = nframes + numel(FORMAT{j});
        FHCR = false(1,numel(format));
        FHCR(FORMAT{j}) = true;
        Cov = Cov + lagged_cov(BOLD(FHCR,good),BOLD(FHCR,good),max(lags));
    end
    
    % Normalize pairwise cross-covariance functions based on entire run
    for k = 1:numel(lags)
        Cov(:,:,k) = Cov(:,:,k)/(nframes - abs(lags(k))*nblocks);
    end
    
    % Parabolic interpolation to get peak lag/correlation
    [pl,pc] = parabolic_interp(Cov,tr);
    pl(abs(pl) > lag_lim) = nan; % Exclude long lags (generally occur when CCF is flat)
    
    % Get zero-lag correlation
    temp = Cov(:,:,lags==0);  % zero-lag correlation
    d = zeros(size(temp));
    d(logical(eye(length(temp)))) = sqrt(diag(temp));
    temp = d^(-1)*temp/d;
    temp = atanh(temp); % Fisher z transform
    temp(isnan(pl)) = nan;
    
    % Add to group running sum
    subj_lags(good,good) = pl;
    subj_ZL(good,good) = temp;
    subj_peak(good,good) = pc;
    
    grp_lags = cat(3,grp_lags,subj_lags);
    grp_lags = nansum(grp_lags,3);
    grp_ZL = cat(3,grp_ZL,subj_ZL);
    grp_ZL = nansum(grp_ZL,3);
    grp_peak = cat(3,grp_peak,subj_peak);
    grp_peak = nansum(grp_peak,3);
    
    % running sum of nans
    grp_lags_nans = grp_lags_nans + isnan(subj_lags);
    grp_ZL_nans = grp_ZL_nans + isnan(subj_ZL);
    grp_peak_nans = grp_peak_nans + isnan(subj_peak);

    toc
    
end

% Compute group averages
grp_lags_mean = grp_lags ./ (numel(subjects) - grp_lags_nans);
grp_peak_mean = grp_peak ./ (numel(subjects) - grp_peak_nans);
grp_ZL_mean = grp_ZL ./ (numel(subjects) - grp_ZL_nans);
grp_ZL_mean = tanh(grp_ZL_mean); % un-fisher z transform

%% Sort group matrices

assns = importdata(''); % import ROI network assignments for sorting TD matrix

% Sort by matrices by lag
[M,sorted_inds1] = sort(nanmean(grp_lags_mean));
assns_sort = assns(sorted_inds1);

grp_lags_temp = grp_lags_mean(sorted_inds1,sorted_inds1);
grp_peak_temp = grp_peak_mean(sorted_inds1,sorted_inds1);
grp_ZL_temp = grp_ZL_mean(sorted_inds1,sorted_inds1);

% Sort by network
[N,sorted_inds2] = sort(assns_sort);
sorted_inds2 = sorted_inds2(find(N,1):end);

grp_lags_mat = grp_lags_temp(sorted_inds2,sorted_inds2);
grp_peak_mat = grp_peak_temp(sorted_inds2,sorted_inds2);
grp_ZL_mat = grp_ZL_temp(sorted_inds2,sorted_inds2);

figure;imagesc(grp_lags_mat,[-.75,.75]);colorbar
figure;imagesc(grp_peak_mat,[-.7,.7]);colorbar
figure;imagesc(grp_ZL_mat,[-.7,.7]);colorbar

%% Make group-level lag projection maps
% Unweighted lag projection
grp_lags_proj_unweighted = nanmean(grp_lags_mean);

% Weighted lag projection (inversely weight lags by correlation magnitude
% to reduce sampling error)
lag_weights = tan((pi/2)*(1-abs(grp_ZL_mean))).^(-2);    % weighted by 1/f^2(r); f(r) = tan[(pi/2)(1-|r|)]
lag_weights(logical(eye(size(lag_weights)))) = 0;   % zero-out diagonal weights
lag_weights(isnan(grp_lags_mean)) = nan;
grp_lags_mean_wghtd = grp_lags_mean.*lag_weights;
grp_lags_proj = nansum(grp_lags_mean_wghtd)./nansum(lag_weights);
