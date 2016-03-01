function x_restored=EM_audio_restoration(x_damaged,iters)

% The function EM_AUDIO_RESTORATION implements Expectation Maximization
% algorithm in order to estimate the restored signal

%INPUT
% x_damaged: the damaged signal that needs to be restored, size: [1 x N]
% iters: the number of iterations of the gibbs sampler, scalar

%OUTPUT
% x_restored: the restored signal, size: [1 x N]

% Reference
% M. Wu and W. Fitzgerald. Analytical approach to
% changepoint detection in laplacian noise. Vision, Image and Signal Processing,
% IEEE Proceedings, 142(3):174–180, Jun 1995.

global N zSize y D B order w L mstart mplusl y1 y2 ylims

%w L
col=zeros(1,N+order);
col(1,2:1+length(x_damaged))=x_damaged;
r=zeros(1,order);
L=sptoeplitz(col,r);
w=zeros(N+order,1);
w(1:N)=x_damaged;

%initialize sigma
init_a = min(x_damaged);
b = max(x_damaged);
sigma_0 = (b-init_a).*randn(1,1) + init_a;

theta_0 = generate_theta(sigma_0);%returns theta_0 as column
sigma_0 = generate_sigma(theta_0');%scalar
z_0 = generate_z(theta_0,sigma_0);%returns z_0 as column

m=zeros(iters,zSize+order+1);
m(1,:)=[z_0' theta_0' sigma_0];

%update w,L
x_damaged(mstart:mplusl-1)=z_0;
col=zeros(1,N+order);
col(1,2:1+length(x_damaged))=x_damaged;
r=zeros(1,order);
L=toeplitz(col,r);
w=zeros(N+order,1);
w(1:N)=x_damaged;

%update D,B
col=zeros(1,N);
col(1,1)=1;
col(1,2:1+length(theta_0))=-theta_0;
r=zeros(1,N);
r(1,1)=1;
K=sptoeplitz(col,r);
KTK=K*K';
B1=KTK(1:length(y1),mstart:mstart+zSize-1);
D=KTK(length(y1)+1:length(y1)+zSize,mstart:mstart+zSize-1);
B2=KTK(length(y1)+zSize+1:end,mstart:mstart+zSize-1);
B=[B1;B2];

%sigma
e=w-L*theta_0;

%MLE for sigma
sigma_estim=sqrt((e'*e)/N);%scalar

h = waitbar(0,'Please wait...');

figure;
for t=1:iters
    
    waitbar(t / iters)
    
    %for fixed z compute theta
    zfixed=m(t,1:zSize);
    
    [Lambda,U] = lu(L'*L);
    estim_theta=U\(Lambda\(L'*w));
    
    col=zeros(1,N);
    col(1,1)=1;
    col(1,2:1+length(estim_theta))=-estim_theta;
    r=zeros(1,N);
    r(1,1)=1;
    K=sptoeplitz(col,r);
    KTK=K*K';
    B1=KTK(1:length(y1),mstart:mstart+zSize-1);
    D=KTK(length(y1)+1:length(y1)+zSize,mstart:mstart+zSize-1);
    B2=KTK(length(y1)+zSize+1:end,mstart:mstart+zSize-1);
    B=[B1;B2];
    
    m(t+1,zSize+1:zSize+order)=estim_theta;
    
    %for fixed theta compute z
    fixedtheta=m(t+1,zSize+1:zSize+order);
    
    %update sigma
    e=w-L*fixedtheta';
    
    %MLE for sigma
    sigma_estim=sqrt((e'*e)/N);%scalar
    
    %compute T
    temp=inv(L'*L);
    temp_size=length(temp);
    
    for i=0:temp_size-1
        diagonal=diag(temp,i);
        column(i+1)=sum(diagonal);
    end
    tcol=[column zeros(1,zSize-length(column))];
    T=toeplitz(tcol,tcol);
    
    %compute q
    q=zeros(zSize,1);
    aug_data=[y1 zeros(1,zSize) y2];
    conv_result=conv(aug_data,T(1,:));
    
    q=conv_result(mstart:mplusl-1);
    
    [Lam,U] = lu(sigma_estim^2*T+D);
    z_estimiplus1=U\(Lam\(-B'*[y1 y2]'-sigma_estim^2*q'));
    m(t+1,1:zSize)=z_estimiplus1;
    
    m(t+1,end)=sigma_estim;
    
    %generate new couple of z,theta
    theta = generate_theta(sigma_0);%returns theta_0 as column
    z = generate_z(theta,sigma_0);%returns z_0 as column
    
    m(t+2,1:end-1)=[z' theta'];
    
    %update w,L
    x_damaged(mstart:mplusl-1)=z;
    col=zeros(1,N+order);
    col(1,2:1+length(x_damaged))=x_damaged;
    r=zeros(1,order);
    L=toeplitz(col,r);
    w=zeros(N+order,1);
    w(1:N)=x_damaged;
    
    e=w-L*theta;
    
    %MLE for sigma
    m(t+2,end)=sqrt((e'*e)/N);%scalar
    
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
    text(mstart+0.2,ylims(2)-ylims(2)/3,sprintf('Data restored with EM'),'FontSize',12);
    number_of_iteration=t
end

close(h)

x_restored=[y1 m(end,1:zSize) y2];

end