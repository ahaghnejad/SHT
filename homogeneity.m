%% Slope Homogeneity Tests:
%
% Swamy test               >>> Swamy, P. A. (1970). Efficient inference in a random coefficient regression model. Econometrica: Journal of the Econometric Society, 38(2), 311-323.
% Delta tests              >>> Pesaran, M. H., & Yamagata, T. (2008). Testing slope homogeneity in large panels. Journal of Econometrics, 142(1), 50-93.
%--------------------------------------------------------------------------
% This code is written by Amin Haghnejad, March 2021 (am.haghnejad@gmail.com)
%**************************************************************************
%%
clear
clc
%% === Read Data ===

data = xlsread('data.xlsx','1');          % An m-by-n matrix of observations with m = N*T and n = the number of all variables (dependent and regressors)
T = 49;                                   % Time dimension
N = 24;                                   % Cross-sectional dimension
alpha = 0.05;                             % Significance level
%
y = data(:,1);                            % y denotes dependent variable
x = data(:,2:end);                        % x denotes a vector of regressors
k = size(x,2);                            % k denotes the number of slope coefficients
%
X = cell(1,N);
Y = cell(1,N);
for i = 1:N
    X{i} = x(((1-T)+(T*i)):(T*i),:);
    Y{i} = y(((1-T)+(T*i)):(T*i),:);
end
I = eye(T);
tau = ones(T,1);
M = I-(tau*tau'/T);                                  
betahat_ols = cell(1,N);
for i = 1:N
    betahat_ols{i} = (X{i}'*M*X{i})\(X{i}'*M*Y{i});  % The ordinary least squares (OLS) estimate of slope coefficients for the individual units (individual regressions)
end
sigsq_hat = zeros(1,N);
for i = 1:N
    sigsq_hat(i) = ((Y{i}-(X{i}*betahat_ols{i}))'*M*(Y{i}-(X{i}*betahat_ols{i})))/(T-k-1);  % The error variances for the individual units based on the OLS estimator
end
L1 = cell(1,N);
L2 = cell(1,N);
for i = 1:N
    L1{i} = (X{i}'*M*X{i})/sigsq_hat(i);
    L2{i} = (X{i}'*M*Y{i})/sigsq_hat(i);
end
sum1 = zeros(k);
sum2 = zeros(k,1);
for i = 1:N
    sum1 = sum1+L1{i};
    sum2 = sum2+L2{i};
end
betahat_wfe = (sum1^-1)*sum2;   % The weighted fixed effects (WFE) estimate of slope coefficients using the error variances based on the OLS estimator
swamyhat = 0;
for i = 1:N
    swamyhat = swamyhat+(betahat_ols{i}-betahat_wfe)'*((X{i}'*M*X{i})/sigsq_hat(i))*(betahat_ols{i}-betahat_wfe);  % Swamy's (1970) test statistic
end
df = k*(N-1);
pvalue_swamyhat = 1-(chi2cdf(swamyhat,df));
CV_swamyhat = chi2inv((1-alpha),df);
%%
L1 = cell(N,1);
L2 = cell(N,1);
for i = 1:N
    L1{i} = X{i}'*M*X{i};
    L2{i} = X{i}'*M*Y{i};
end
sum1 = zeros(k);
sum2 = zeros(k,1);
for i = 1:N
    sum1 = sum1+L1{i};
    sum2 = sum2+L2{i};
end
betahat_fe = (sum1^-1)*sum2;    % The fixed effects (FE) estimate of slope coefficients
sigsq_tilde = zeros(1,N);
for i = 1:N
    sigsq_tilde(i) = ((Y{i}-(X{i}*betahat_fe))'*M*(Y{i}-(X{i}*betahat_fe)))/(T-1);   % The error variances for the individual units based on the standard fixed effects (FE) estimator
end
L1 = cell(1,N);
L2 = cell(1,N);
for i = 1:N
    L1{i} = (X{i}'*M*X{i})/sigsq_tilde(i);
    L2{i} = (X{i}'*M*Y{i})/sigsq_tilde(i);
end
sum1 = zeros(k);
sum2 = zeros(k,1);
for i = 1:N
    sum1 = sum1+L1{i};
    sum2 = sum2+L2{i};
end
betatilde_wfe = (sum1^-1)*sum2;   % The weighted fixed effects (WFE) estimate of slope coefficients using the error variances based on the FE estimator
swamytilde = 0;
for i = 1:N
    swamytilde = swamytilde+(betahat_ols{i}-betatilde_wfe)'*((X{i}'*M*X{i})/sigsq_tilde(i))*(betahat_ols{i}-betatilde_wfe);
end
delta = (N^(1/2))*(((N^-1)*swamytilde-k)/((2*k)^(1/2)));    % Delta statistic proposed by Pesaran and Yamagata (2008)
EZ = k;
VZ = (2*k)*(T-k-1)/(T+1);
deltaadj = (N^(1/2))*(((N^-1)*swamytilde-EZ)/(VZ^(1/2)));   % Adjusted Delta statistic proposed by Pesaran and Yamagata (2008)
if delta > 0
    pvalue_delta = 2*(1-(normcdf(delta,0,1)));
else
    pvalue_delta = 2*(normcdf(delta,0,1));
end
UCV_delta = norminv((1-(alpha/2)),0,1);
LCV_delta = norminv((alpha/2),0,1);

if deltaadj > 0
    pvalue_deltaadj = 2*(1-(normcdf(deltaadj,0,1)));
else
    pvalue_deltaadj = 2*(normcdf(deltaadj,0,1));
end
UCV_deltaadj = norminv((1-(alpha/2)),0,1);
LCV_deltaadj = norminv((alpha/2),0,1);
%%
disp('  ------------------------------------------------------')
disp('The results of tests for slope heterogeneity')
disp('  ------------------------------------------------------')
fprintf('\n The number of time periods: %2.0f \n',T);
fprintf('\n The number of cross-sections: %2.0f \n',N);
fprintf('\n Null hypothesis: Slope homogeneity \n');
disp(' ')
%
disp('  ------------------------------------------------------')
disp('  Swamy''s test ---- Swamy (1970)')
disp('  ------------------------------------------------------')
fprintf('\n  Test Statistic \t\t\t d.f.       \t   p-value \n');
fprintf('\n  %4.6f          \t\t %1.0f        \t       %4.6f \n',swamyhat,df,pvalue_swamyhat);
disp(' ')
%%
disp('  ------------------------------------------------------')
disp('  Delta tests ---- Pesaran and Yamagata (2008))')
disp('  ------------------------------------------------------')
fprintf('\n  Test                \t  Statistic    \t\t p-value \n\n');
fprintf('\n  Delta               \t  %4.6f         \t %4.6f \n\n',delta,pvalue_delta);
fprintf('\n  Adjusted-Delta      \t  %4.6f         \t %4.6f \n\n',deltaadj,pvalue_deltaadj);
%%





