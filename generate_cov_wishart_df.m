function C=generate_cov_wishart_df(N,Df,Sig)
% inspired by Alexandre Barachant's code
% https://github.com/alexandrebarachant/covariancetoolbox/blob/master/lib/simulation/generate_wishart_set.m

if nargin<3
    [Q,~] = qr(randn(N));
    Sig = Q'*diag(5*rand(N,1))*Q;
end

C(:,:) = wishrnd(Sig,Df);
