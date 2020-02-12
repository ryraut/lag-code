%% Main script used to compute spectral content (intrinsic timescales)
%
%
%% Setup
% Set parameters
num_nodes = 264;    % number of time series 
outdir = '';    % set directory for saving out images here
mkdir(outdir)
lags = -6:6;    % range of TR shifts; should be sufficient to allow for all autocovariance functions (ACFs) to decay below .5

% Specify data parameters
subjects = importdata('');
tr = 2;
motion_thresh = .2;    % important: must match motion criteria used during preproc

min_block_durn = (max(lags)+1)*tr;   % min. block duration (in seconds)

%% Loop over subjects

% initialize group matrices
grp_acfs = single(nan(num_nodes,numel(lags),numel(subjects))); % keep all ACFs for all subjects for stats

for s = 1:numel(subjects)
    tic
    subj = subjects{s};
    disp(['Processing ' subj]);
    
    % initialize subject matrices
    subj_lags = single(nan(num_nodes)); % peak lags
    subj_ZL = subj_lags;   % zero-lag correlation
    subj_peak = subj_lags; % peak correlation (correlation at optimal lag)
    
    
    BOLD = importdata(''); % read in time series matrix
    good = importdata(''); % read in spatial mask if desired
    
    % read in temporal mask/motion time series (e.g., FD or DVARS)
    format = dlmread('') <= motion_thresh;
    
    % ignore pre-steady-state frames
    format(1:2) = false; % ignore first X frames
    
    FORMAT = create_blocks(format,min_block_durn,tr);
    
    %% Construct ACF
    ACFs = single(zeros(num_nodes,numel(lags)));
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
        for i = 1:sum(good)
            ACFs(i,:) = ACFs(i,:) + squeeze(lagged_cov(BOLD(FHCR,i),BOLD(FHCR,i),max(lags)))';
        end
    end
    
    % Normalize ACFs based on entire run
    for k = 1:numel(lags)
        ACFs(:,k) = ACFs(:,k)/(nframes - abs(lags(k))*nblocks);
    end
    ACFs = bsxfun(@rdivide,ACFs,ACFs(:,lags==0));
        
    % Store
    grp_acfs(good,:,s) = ACFs;

    toc
end

acf_mean = nanmean(grp_acfs,3);
acf_mean = bsxfun(@rdivide,acf_mean,acf_mean(:,lags==0));

% fit
hwhm = acf_hwhm(acf_mean',tr); % group-wise

% subject-wise
hwhms = zeros(num_nodes,numel(subjects));
for s = 1:numel(subjects)
    tic
    hwhms(:,s) = acf_hwhm(grp_acfs(:,:,s)',tr);
    toc
end
