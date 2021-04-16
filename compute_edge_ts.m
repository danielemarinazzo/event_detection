function [idxpeak,rss,rssr,corr_mats] = compute_edge_ts(ts,numshift)
%numshift counts how many times we actually circularly shift the time
%series, to control for the level of randomness
if nargin<2
    numshift=1;
end

% modified version of example_script.m, taking ts as argument

% zscore time series
z = zscore(ts);

% number of time points/nodes
[t,n] = size(z);

% upper triangle indices (node pairs = edges)
[u,v] = find(triu(ones(n),1));

% edge time series
ets = z(:,u).*z(:,v);

% calculate rss
rss = sum(ets.^2,2).^0.5;

% repeat with randomized time series
numrand = 100;

% initialize array for null rss
rssr = zeros(t,numrand);

% perform numrand randomizations
for irand = 1:numrand
    
    % create circularly shifted time series
    zr = z;
    if (irand>numshift || numshift==1) 
        for i = 1:n
            zr(:,i) = circshift(zr(:,i),randi(t));
        end
    end
    
    % edge time series with circshift data
    etsr = zr(:,u).*zr(:,v);
    
    % calcuate rss
    rssr(:,irand) = sum(etsr.^2,2).^0.5;
    
end

% calculate p-value
p = zeros(t,1);
for i = 1:t
    p(i) = mean(rssr(:) >= rss(i));
end

% apply statistical cutoff
pcrit = 0.001;

% find frames that pass statistical test
idx = find(p < pcrit);

% identify contiguous segments of frames that pass statistical test
dff = idx' - (1:length(idx));
unq = unique(dff);
nevents = length(unq);

% find the peak rss within each segment
idxpeak = zeros(nevents,1);
for ievent = 1:nevents
    idxevent = idx(dff == unq(ievent));
    rssevent = rss(idxevent);
    [~,idxmax] = max(rssevent);
    idxpeak(ievent) = idxevent(idxmax);
end

% get activity at peak
tspeaks = z(idxpeak,:);

% get co-fluctuation at peak
etspeaks = tspeaks(:,u).*tspeaks(:,v);



%% calculate mean co-fluctuation (edge time series) across all peaks
mu = nanmean(etspeaks,1);

% represent in matrix form
mat = zeros(n);
mat(triu(ones(n),1) > 0) = mu;
mat = mat + mat';

mat_full=corrcoef(ts);
corr_mats=corr(mat(triu(ones(n),1) > 0),mat_full(triu(ones(n),1) > 0),'Type','Spearman');

