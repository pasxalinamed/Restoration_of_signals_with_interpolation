function   p_sigma = conditional_sigma_PDF(variable,thetafixed)

% The function CONDITIONAL_SIGMA_PDF computes the conditional prob. density of
% std. deviation (sigma) given an instant of sigma (variable) and  AR parameters 
% (thetafixed)

%INPUT
% st.deviation (sigmafixed), size: [1 x 1]
% AR parameters (variable): size [order x 1]

%OUTPUT
% conditional prob. density (p_sigma), size: [1 x 1]

% Reference
% M. Wu and W. Fitzgerald. Analytical approach to
% changepoint detection in laplacian noise. Vision, Image and Signal Processing,
% IEEE Proceedings, 142(3):174-180, Jun 1995.

global  N w L

e=w-L*thetafixed';
p_sigma=(e'*e)^(N/2)*variable^(-N)*exp(-e'*e*1/(2*variable^2));

end