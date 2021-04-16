%%
% This script generates covariance matrices with different degrees of freedom,
% and from there time series with an autoregressive model, then it computes
% the Spearman correlation between the edge FC and the full time series FC
%%
clear;clc;close all
% change the directory below, use https://github.com/tapios/arfit
addpath('/home/daniele.marinazzo@rc.fhsc.net/Dropbox/code/toolboxes_in_use/arfit');
warning off all
% number of time series to mimic the real data used in this example
N = 333;
t = 818;
AR_order=0;
nD=10;
Dfmax=850;
n_avg=100;
Dtot=round(linspace(N+1,Dfmax,nD));
randtot=[0 15 30];
CM=zeros(nD,length(randtot));NS=CM;
figure;iplot=0;
for irand=1:length(randtot)
    numshift=randtot(irand);
    for D = Dtot
        corr_mats=0;nsig=0;
        parfor i_avg=1:n_avg
            warning off all
            A = .085*randn(N,AR_order*N);
            w  = .01*randn(N,1);
            C  = generate_cov_wishart_df(N,D);
            ts = arsim(w,A,C,t,1000);
            zr = zscore(ts);
            [idxpeak,~,~,c_trial] = compute_edge_ts(zr,numshift);
            corr_mats=corr_mats+c_trial;
            nsig=nsig+length(idxpeak)/t;
        end
        CM(Dtot==D,irand)=corr_mats/n_avg;
        NS(Dtot==D,irand)=nsig/n_avg;
    end
end
%%
figure;
yyaxis left;plot(Dtot(1:nD),CM(1:nD,1:3),'LineWidth',2);
yyaxis right;plot(Dtot(1:nD),NS(1:nD,1:3),'LineWidth',2);
yyaxis left
title('Covariance matrix, N = 333, 100 realizations')
xlabel('Df for Wishart distribution')
ylabel('Spearman correlation between edge FC and Pearson FC')
yyaxis right
ylabel('percentage of significant RSS')
%xlim([300 900]);