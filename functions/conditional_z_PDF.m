function p_z = conditional_z_PDF(variable,sigmafixed)

% The function CONDITIONAL_Z_PDF computes the conditional prob. density of
% missing samples (z) given an instant of z (variable) and std. deviation (sigmafixed)

%INPUT
% missing data (variable): size [zSize x 1]
% st.deviation (sigmafixed), size: [1 x 1]

%OUTPUT
% conditional prob. density (p_z), size: [1 x 1]

% Reference
% M. Wu and W. Fitzgerald. Analytical approach to
% changepoint detection in laplacian noise. Vision, Image and Signal Processing,
% IEEE Proceedings, 142(3):174-180, Jun 1995.

global y1 y2 L D B

% inverse covariance matrix
C_inv=D/(sigmafixed^2);

%compute z with Band LU Decomposition
[L,U] = lu(D);
z_estim=U\(L\(-B'*[y1 y2]'));

numerator=(variable'-z_estim)'*C_inv*(variable'-z_estim);
denominator=(2*sigmafixed^2);
p_z=exp(-numerator*(1/denominator));

end