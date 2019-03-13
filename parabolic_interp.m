%This function uses parabolic interpolation to find the lag using the
%extremum of the lagged cross correlation.

%The function first uses the
%sign of the cross correlation at zero to decide whether to find a maximum
%or minmum. Next, we look for the global max/min.

%lcc is the empirical lagged covariance curve, lags is a vector with the timepoints
%in each temporal direction (e.g. -8:2:8 for +/- 8 seconds with a 2 second TR). 

function [peak_lag,peak_cov] = parabolic_interp(lcc,tr)
	s = size(lcc);
	peak_lag = nan([1,s(1)*s(2)]);
	peak_cov = peak_lag;
	
    % linearize
	lcc = reshape(lcc,[s(1)*s(2),s(3)])';
    
    % find index of extremum (max or min determined by sign at zero-lag)
	[~,I]= max(bsxfun(@times,lcc,sign(lcc((s(3)+1)/2,:))),[],1);
    
    % ensure extremum is not at an endpoint (this would preclude parabolic interpolation)
	use = I>1 & I<s(3);
	lcc = lcc(:,use);
    
	% place peaks at center
	x0 = I(use) - (s(3)+1)/2;

	% set up three-point ccf for interpolation (y1,y2,y3)
	i = sub2ind([size(lcc),sum(use)],I(use),1:sum(use));
	lcc = [lcc(i-1);lcc(i);lcc(i+1)];

    % fit parabola: tau = TR * (y1-y3) / (2*(y1-2y2+y3))
	b = (lcc(3,:) - lcc(1,:))/2;
    a = (lcc(1,:) + lcc(3,:) - 2*lcc(2,:))/2;
	peak_lag(use) =  (-b./(2*a));
    
    % construct parabola to get covariance (y = ax^2 + bx + c)
	peak_cov(use) = a.*(peak_lag(use).^2) + b.*peak_lag(use) + lcc(2,:);
    
    % put back TR information
	peak_lag(use) = (peak_lag(use) + x0)*tr;

	peak_lag = reshape(peak_lag,[s(1) s(2)]);
	peak_cov = reshape(peak_cov,[s(1) s(2)]);
	
end


