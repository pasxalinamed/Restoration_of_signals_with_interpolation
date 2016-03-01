function theta = generate_theta(sigmafixed)

% The function GENERATE_THETA is used during the procedure of Gibbs sampling to
% draw samples from a multivariate Gaussian density. These samples considered representative of
% the conditional density of theta given the variable sigma (sigmafixed - st. deviation)

%INPUT
% st.deviation (sigmafixed), size: [1 x 1]

%OUTPUT
% AR parameters (theta), size: [order x 1]

% Reference
% M. Wu and W. Fitzgerald. Analytical approach to
% changepoint detection in laplacian noise. Vision, Image and Signal Processing,
% IEEE Proceedings, 142(3):174-180, Jun 1995.

global order w L

% inverse covariance matrix
Cthetainv=(L'*L)/(sigmafixed^2);

% compute Cthetainv square rout S via Cholesky decomposition
S=chol(Cthetainv);

% draw AR parameters (theta) samples from a multivariate Gaussian density
mu=0;
stan_dev=1;
theta = mu + stan_dev.*randn(order,1);
n=theta-mu;
e=S*n;
u=S\e+mu;

%compute n with Band LU Decomposition
[Lambda,U] = lu(S');
n=U\(Lambda\(u));

%compute theta with Band LU Decomposition
[Lambda,U] = lu(L'*L);
estim_theta=U\(Lambda\(L'*w));

theta=estim_theta+n;

end