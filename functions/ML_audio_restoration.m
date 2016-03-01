function x_restored=ML_audio_restoration(x_damaged)

% The function ML_AUDIO_RESTORATION implements Maximum Likelihood algorithm
% in order to estimate the restored signal

%INPUT
% x_damaged: the damaged signal that needs to be restored, size: [1 x N]

%OUTPUT
% x_restored: the restored signal, size: [1 x N]

% Reference
% M. Wu and W. Fitzgerald. Analytical approach to
% changepoint detection in laplacian noise. Vision, Image and Signal Processing,
% IEEE Proceedings, 142(3):174–180, Jun 1995.

global N zSize y D B L order w y1 y2 mstart mplusl ylims

%update w, L (needed every time that data x changes)
col=zeros(1,N+order);
col(1,2:1+length(x_damaged))=x_damaged;
r=zeros(1,order);
L=sptoeplitz(col,r);
w=zeros(N+order,1);
w(1:N)=x_damaged;

%initialize sigma value
init_a = 0;
b = var(x_damaged);
sigma_0 = (b-init_a).*randn + init_a;

theta_0 = generate_theta(sigma_0);%returns theta_0 as column
sigma_0 = generate_sigma(theta_0');%scalar
z_0 = generate_z(theta_0,sigma_0);%returns z_0 as column

m(1,:)=[z_0' theta_0' sigma_0];

%update w, L (needed every time that data x changes)
x_damaged(mstart:mplusl-1)=z_0;
col=zeros(1,N+order);
col(1,2:1+length(x_damaged))=x_damaged;
r=zeros(1,order);
L=sptoeplitz(col,r);
w=zeros(N+order,1);
w(1:N)=x_damaged;

t=1;
difference(t)=1;
partial_likelihood=ones(2,1);

figure;
while difference>0.0001
    
    zfixed=m(t,1:zSize);
    
    x_damaged(mstart:mplusl-1)=zfixed;
    col=zeros(1,N+order);
    col(1,2:1+length(x_damaged))=x_damaged;
    r=zeros(1,order);
    L=sptoeplitz(col,r);
    w=zeros(N+order,1);
    w(1:N)=x_damaged;
    
    [Lambda,U] = lu(L'*L);
    estim_theta=U\(Lambda\(L'*w));
    m(t+1,zSize+1:zSize+order)=estim_theta;
    
    fixedtheta=estim_theta;
    
    % update B,D (needed every time that theta, which expresses the AR coefficients, changes)
    col=zeros(1,N);
    col(1,1)=1;
    col(1,2:1+length(fixedtheta))=-fixedtheta;
    r=zeros(1,N);
    r(1,1)=1;
    K=sptoeplitz(col,r);
    KTK=K*K';
    B1=KTK(1:length(y1),mstart:mstart+zSize-1);
    D=KTK(length(y1)+1:length(y1)+zSize,mstart:mstart+zSize-1);
    B2=KTK(length(y1)+zSize+1:end,mstart:mstart+zSize-1);
    B=[B1;B2];
    
    [Lam,U] = lu(D);
    z_estim=U\(Lam\(-B'*[y1 y2]'));
    m(t+1,1:zSize)=z_estim;
    
    %update w, L (needed every time that data x changes)
    x_damaged(mstart:mplusl-1)=z_estim;
    col=zeros(1,N+order);
    col(1,2:1+length(x_damaged))=x_damaged;
    r=zeros(1,order);
    L=sptoeplitz(col,r);
    w=zeros(N+order,1);
    w(1:N)=x_damaged;
    e=w-L*fixedtheta;
    
    %MLE for sigma
    sigma_estim=sqrt((e'*e)/N);%scalar
    m(t+1,end)=sigma_estim;
    
    %to minimize (it is equivalent to maximize prob.density)
    partial_likelihood(t)=e'*e;
    
    if t==1
        difference(t)=abs(partial_likelihood(t)-1);
    else
        difference(t)=abs(partial_likelihood(t)-partial_likelihood(t-1));
    end
    
    clf
    
    plot([y1 m(t+1,1:zSize) y2],'LineWidth',1);
    xlim([0 N])
    ylim([ylims])
    set(gca,'FontSize',13)
    ylims = get(gca,'YLim');
    hold on
    plot([mstart mstart],ylims, 'r:','LineWidth',2);
    plot([mplusl mplusl],ylims, 'r:','LineWidth',2);
    xlabel('i');
    ylabel('x_i');
    text(mstart+0.2,ylims(2)-ylims(2)/2,sprintf('Data restored with ML'),'FontSize',12);
    
    t=t+1;
end

str=sprintf('Converged after %i iterations', t)

x_restored=[y1 m(end,1:zSize) y2];

end