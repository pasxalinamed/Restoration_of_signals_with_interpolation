function    p_theta = conditional_theta_PDF(variable,sigmafixed)

% The function CONDITIONAL_THETA_PDF computes the conditional prob. density of
% AR parameters (theta) given an instant of theta (variable) and std. deviation (sigmafixed)

%INPUT
% AR parameters (variable): size [order x 1]
% st.deviation (sigmafixed), size: [1 x 1]

%OUTPUT
% conditional prob. density (p_theta), size: [1 x 1]

% Reference
% M. Wu and W. Fitzgerald. Analytical approach to
% changepoint detection in laplacian noise. Vision, Image and Signal Processing,
% IEEE Proceedings, 142(3):174-180, Jun 1995.

global  order N y1 y2 w L zSize

x_damaged=[y1 sparse(1,zSize) y2];
col=zeros(1,N+order);
col(1,2:1+length(x_damaged))=x_damaged;
r=zeros(1,order);
L=sptoeplitz(col,r);
w=zeros(N+order,1);
w(1:N)=x_damaged;
estim_theta=(L'*L)\L'*w;
Cthetainv=(L'*L)/(sigmafixed^2);
Q_theta=w'*L*estim_theta-2*w'*L*variable'+...
    variable*Cthetainv*(sigmafixed^2)*variable';
p_theta=(sigmafixed^(-order))*exp(-Q_theta*(1/(2*sigmafixed^2)));

end