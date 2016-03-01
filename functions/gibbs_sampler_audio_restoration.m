function [x_restored,m]=gibbs_sampler_audio_restoration(x_damaged,iters,burnIn)

% The function GIBBS_SAMPLER_AUDIO_RESTORATION implements the gibbs sampler
% to approach the conditional density of the missing data (z)

%INPUT
% x_damaged: the damaged signal that needs to be restored, size: [1 x N]
% iters: the number of iterations of the gibbs sampler, scalar
% burnIn: the number of the samples of the chain that will be taken into
% account, because they are considered more reliable samples of the distribution

%OUTPUT
% x_restored: the restored signal, size: [1 x N]
% m: results of the last burnIn estimations [burnIn x (zSize+order+1)]

% Reference
% M. Wu and W. Fitzgerald. Analytical approach to
% changepoint detection in laplacian noise. Vision, Image and Signal Processing,
% IEEE Proceedings, 142(3):174-180, Jun 1995.

global order N zSize w L mstart mplusl y1 y2 D B ylims

% initialize chain
w=[x_damaged'; zeros(order,1)];
col=zeros(1,N+order);
col(1,2:1+length(x_damaged))=x_damaged';
r=zeros(1,order);
L=sptoeplitz(col,r);

sigma_0 = gamrnd((N-1)/2,1,1,1);
theta_0 = generate_theta(sigma_0);%returns theta_0 as column
sigma_0 = generate_sigma(theta_0');%scalar
z_0 = generate_z(theta_0,sigma_0);%returns z_0 as column

x_damaged(mstart:mplusl-1)=z_0;

h = waitbar(0,'Please wait...');
a_zet=zeros(iters,1);

currentm=zeros(1,zSize+order+1);
nextm=zeros(1,zSize+order+1);

currentm=[z_0' theta_0' sigma_0];
for t=1:iters-burnIn
    
    waitbar(t/iters)
    
    %Gibbs sampler for z
    thetafixed=currentm(1,zSize+1:zSize+order);%consider theta fixed
    sigmafixed=currentm(1,end);%sigma fixed
    conditionalVar=currentm(1,1:zSize);%old z to be checked
    
    % generate new candidate z
    candidateVar=generate_z(thetafixed,sigmafixed)';
    
    while conditionalVar==candidateVar
        candidateVar=generate_z(thetafixed,sigmafixed)';
    end
    
    U=rand();
    
    p_z1 = conditional_z_PDF(candidateVar,sigmafixed);
    p_z2 = conditional_z_PDF(conditionalVar,sigmafixed);
    
    a_zet(t,1)=min([1 p_z1/p_z2]);
    
    if U<=a_zet(t,1)
        nextm(1,1:zSize)=candidateVar;
        x_damaged(mstart:mplusl-1)=candidateVar;
        col=zeros(1,N+order);
        col(1,2:1+length(x_damaged))=x_damaged;
        r=zeros(1,order);
        L=sptoeplitz(col,r);
        w=zeros(N+order,1);
        w(1:N)=x_damaged;
    else
        nextm(1,1:zSize)=conditionalVar;
        x_damaged(mstart:mplusl-1)=conditionalVar;
        col=zeros(1,N+order);
        col(1,2:1+length(x_damaged))=x_damaged;
        r=zeros(1,order);
        L=sptoeplitz(col,r);
        w=zeros(N+order,1);
        w(1:N)=x_damaged;
    end
    
    %Gibbs sampler for theta
    sigmafixed=currentm(1,end);
    conditionalVar=currentm(1,zSize+1:zSize+order);
    
    candidateVar=generate_theta(sigmafixed)';
    
    while conditionalVar==candidateVar
        candidateVar=generate_theta(sigmafixed)';
    end
    
    U=rand();
    
    p_theta1 = conditional_theta_PDF(candidateVar,sigmafixed);
    p_theta2 = conditional_theta_PDF(conditionalVar,sigmafixed);
    
    a=min([1 p_theta1/p_theta2]);
    
    if U<=a
        nextm(1,1+zSize:zSize+order)=candidateVar;
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
    else
        nextm(1,1+zSize:zSize+order)=conditionalVar;
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
    end
    
    %Gibbs sampler for sigma
    thetafixed=currentm(1,zSize+1:zSize+order);
    conditionalVar=currentm(1,end);
    candidateVar=generate_sigma(thetafixed);
    
    while conditionalVar==candidateVar
        candidateVar=generate_sigma(thetafixed);
    end
    
    U=rand();
    
    p_sigma1 = conditional_sigma_PDF(candidateVar,thetafixed);
    p_sigma2 = conditional_sigma_PDF(conditionalVar,thetafixed);
    
    a=min([1 p_sigma1/p_sigma2]);
    
    if U<=a
        nextm(1,end)=candidateVar;
    else
        nextm(1,end)=conditionalVar;
    end
    
    currentm=nextm;
end

m(1,:)=currentm;

for t=1:burnIn
    
    waitbar(t/iters)
    
    %Gibbs sampler for z
    thetafixed=m(t,zSize+1:zSize+order);%consider theta fixed
    sigmafixed=m(t,end);%sigma fixed
    conditionalVar=m(t,1:zSize);%old z to be checked
    
    % generate new candidate z
    candidateVar=generate_z(thetafixed,sigmafixed)';
    
    while conditionalVar==candidateVar
        candidateVar=generate_z(thetafixed,sigmafixed)';
    end
    
    U=rand();
    
    p_z1 = conditional_z_PDF(candidateVar,sigmafixed);
    p_z2 = conditional_z_PDF(conditionalVar,sigmafixed);
    
    a_zet(t,1)=min([1 p_z1/p_z2]);
    
    if U<=a_zet(t,1)
        m(t+1,1:zSize)=candidateVar;
        x_damaged(mstart:mplusl-1)=candidateVar;
        col=zeros(1,N+order);
        col(1,2:1+length(x_damaged))=x_damaged;
        r=zeros(1,order);
        L=sptoeplitz(col,r);
        w=zeros(N+order,1);
        w(1:N)=x_damaged;
    else
        m(t+1,1:zSize)=conditionalVar;
        x_damaged(mstart:mplusl-1)=conditionalVar;
        col=zeros(1,N+order);
        col(1,2:1+length(x_damaged))=x_damaged;
        r=zeros(1,order);
        L=sptoeplitz(col,r);
        w=zeros(N+order,1);
        w(1:N)=x_damaged;
    end
    
    %Gibbs sampler for theta
    sigmafixed=m(t,end);
    conditionalVar=m(t,zSize+1:zSize+order);
    
    candidateVar=generate_theta(sigmafixed)';
    
    while conditionalVar==candidateVar
        candidateVar=generate_theta(sigmafixed)';
    end
    
    U=rand();
    
    p_theta1 = conditional_theta_PDF(candidateVar,sigmafixed);
    p_theta2 = conditional_theta_PDF(conditionalVar,sigmafixed);
    
    a=min([1 p_theta1/p_theta2]);
    
    if U<=a
        m(t+1,1+zSize:zSize+order)=candidateVar;
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
    else
        m(t+1,1+zSize:zSize+order)=conditionalVar;
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
    end
    
    %Gibbs sampler for sigma
    thetafixed=m(t+1,zSize+1:zSize+order);
    conditionalVar=m(t,end);
    candidateVar=generate_sigma(thetafixed);
    
    while conditionalVar==candidateVar
        candidateVar=generate_sigma(thetafixed);
    end
    
    U=rand();
    
    p_sigma1 = conditional_sigma_PDF(candidateVar,thetafixed);
    p_sigma2 = conditional_sigma_PDF(conditionalVar,thetafixed);
    
    a=min([1 p_sigma1/p_sigma2]);
    
    if U<=a
        m(t+1,end)=candidateVar;
    else
        m(t+1,end)=conditionalVar;
    end
end

close(h)

x_restored=[y1 mean(m(end-burnIn:end,1:zSize)) y2];

figure;
plot(x_restored,'k','LineWidth',1);
xlim([0 N]);
ylim(ylims)
set(gca,'FontSize',13)
hold on
plot([mstart mstart],ylims, 'r:','LineWidth',2);
plot([mplusl mplusl],ylims, 'r:','LineWidth',2);
xlabel('i');
ylabel('x_i');
text(mstart+0.5,ylims(2)-ylims(2)/2+0.1,sprintf('Data restored with Gibbs sampling'),'FontSize',12);

end