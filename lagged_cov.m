% Computes (unnormalized) cross-covariance function out to +/- L lags 
% (TR shifts) between each column of Avg1 and Avg2

% Avg1 and Avg2 are time x region matrices/vectors

% See tdmx_template.m for appropriate normalization of lagged_cov output

function  r = lagged_cov(Avg1,Avg2,L)
    
	L1 = size(Avg1,2);
	L2 = size(Avg2,2);
	r = single(zeros(L1,L2,2*L+1));

	k = 1;
	for i = -L:L
		tau = abs(i);
         
        if i >=0
            Avg1_lagged = Avg1(1:end-tau,:);
            Avg2_lagged = Avg2(1+tau:end,:);
        else
            Avg1_lagged = Avg1(1+tau:end,:);
            Avg2_lagged = Avg2(1:end-tau,:);
        end    
        
		r(:,:,k) = Avg1_lagged'*Avg2_lagged;
		k = k+1;
	end

end
