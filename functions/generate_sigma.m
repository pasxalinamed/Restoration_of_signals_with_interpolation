function sigma = generate_sigma(thetafixed)

% The function GENERATE_SIGMA is used during the procedure of Gibbs sampling to
% draw samples from the conditional density of sigma given the variable theta
% (thetafixed - AR paratemers). Samples are drawn at first step from a gamma
% prob. density and then they are scaled to fit the required sigma
% formulation.

%INPUT
% AR parameters (thetafixed), size: [1 x order]

%OUTPUT
% st.deviation (sigma), size: [1 x 1]

% Reference
% M. Wu and W. Fitzgerald. Analytical approach to
% changepoint detection in laplacian noise. Vision, Image and Signal Processing,
% IEEE Proceedings, 142(3):174-180, Jun 1995.

global  N w L

%generate sigma from a gamma prob. density
a = (N-1)/2;
b = 1;
u = gamrnd(a,b,1,1);

e=w-L*thetafixed';

%scale the result to fit the required sigma formulation
sigma=sqrt((e'*e)/2)*(1/sqrt(u));

end