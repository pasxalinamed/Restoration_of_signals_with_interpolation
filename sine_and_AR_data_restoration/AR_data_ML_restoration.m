%Tests the ML algorithm on AR modeled data of 4th order

clear all;close all;
global N x_damaged zSize y order y1 y2 mstart mplusl ylims

order=4; %order of AR model
N=1000; %size of data
mstart=300; %data removal starts from here
mplusl=499; %data removed until here
zSize=mplusl-mstart;

%Generate original AR modeled data
c = fir1(256, 0.5);
[alpha,p0] = lpc(c,order);
e = sqrt(p0)*randn(1,N); 
x = filter(1,alpha,e);
original_z=x(mstart:mplusl-1);

% plot original data
figure;
plot(x,'LineWidth',1);
xlim([0 N])
ylim([-0.2 0.2])
set(gca,'FontSize',13)
ylims = get(gca,'YLim');
hold on
plot([mstart mstart],ylims, 'r:','LineWidth',2);
plot([mplusl mplusl],ylims, 'r:','LineWidth',2);
xlabel('i');
ylabel('x_i');
text(mstart+0.2,ylims(1)+0.03,sprintf('Removed\n data'),'FontSize',12);
set(gca,'XTick',0:N/10:N)

%remove part of the original data
y1=x(1:mstart-1);
y2=x(mplusl:end);
y=[y1 y2];
x_damaged=x;
x_damaged(mstart:mplusl-1)=0;

%restore data
x_restored=ML_audio_restoration(x_damaged);

print -deps mlrestoredARData