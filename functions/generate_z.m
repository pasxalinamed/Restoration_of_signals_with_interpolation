function z = generate_z(thetafixed,sigmafixed)

% The function GENERATE_Z is used during the procedure of Gibbs sampling to
% draw samples from a multivariate Gaussian density. These samples considered representative of
% the conditional density of z given the variables theta (thetafixed - AR paratemers)
% and sigma (sigmafixed - st. deviation)

%INPUT
% st.deviation (sigmafixed), size: [1 x 1]
% AR parameters (thetafixed), size: [order x 1]

%OUTPUT
% missing data (z): size [zSize x 1]

% Reference
% M. Wu and W. Fitzgerald. Analytical approach to
% changepoint detection in laplacian noise. Vision, Image and Signal Processing,
% IEEE Proceedings, 142(3):174-180, Jun 1995.

global N zSize D B mstart y1 y2

col=zeros(1,N);
col(1,1)=1;
col(1,2:1+length(thetafixed))=-thetafixed;
r=zeros(1,N);
r(1,1)=1;
K=sptoeplitz(col,r);
KTK=K*K';
B1=KTK(1:length(y1),mstart:mstart+zSize-1);
D=KTK(length(y1)+1:length(y1)+zSize,mstart:mstart+zSize-1);
B2=KTK(length(y1)+zSize+1:end,mstart:mstart+zSize-1);
B=[B1;B2];

% inverse covariance matrix
C_inv=D/(sigmafixed^2);

% compute C_inv square rout S via Cholesky decomposition
S=chol(C_inv);

% draw missing data (z) samples from a multivariate Gaussian density
mu=0;
stan_dev=1;
z = mu + stan_dev.*randn(zSize,1);
n=z-mu;
e=S*n;
u=S\e+mu;

%compute n with Band LU Decomposition
[L,U] = lu(S');
n=U\(L\(u));

%compute z with Band LU Decomposition
[L,U] = lu(D);
z_estim=U\(L\(-B'*[y1 y2]'));

z=z_estim+n;

end