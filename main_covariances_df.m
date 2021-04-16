%%
% This script generates covariance matrices with different degrees of freedom,
% and from there time series with an autoregressive model
%%
clear;clc;close all
% change the directory below, use https://github.com/tapios/arfit
addpath('/home/daniele.marinazzo@rc.fhsc.net/Dropbox/code/toolboxes_in_use/arfit');
warning off all
% number of time series to mimic the real data used in this example
N = 333;
t = 818;
AR_order=0;
nD=4;
Dfmax=N*3;
Dtot=round(linspace(N+1,Dfmax,nD));
figure;
for iDF=1:length(Dtot)
    D=Dtot(iDF);
    A = .085*randn(N,AR_order*N);
    w  = .01*randn(N,1);
    C  = generate_cov_wishart_df(N,D);
    ts = arsim(w,A,C,t,1000);
    zr = zscore(ts);
    [idxpeak,rss,rssr,corr_mats] = compute_edge_ts(ts);
    
    %% plot rss time series, null, and significant peaks
    
    subplot(4,1,iDF);
    ph = plot(1:t,rssr,'color',ones(1,3)*0.65);
    hold on;
    qh = plot(idxpeak,rss(idxpeak),'r*',1:t,rss,'k','linewidth',2);
    xlim([1,t])
    xlabel('frame'); ylabel('rss');
    legend([ph(1); qh],'null','significant','orig');
    %     title(['Df = ' num2str(D) ', corr edge-full FC = ' num2str(corr_mats)])
    title(['Df = ' num2str(D) ', N = ' num2str(N)])
end