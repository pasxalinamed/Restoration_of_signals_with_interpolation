%Tests the Gibbs sampling algorithm on AR modeled data of 4th order

clear;close all;
global order N x_damaged zSize y mstart mplusl y1 y2 ylims

order=4; %order of AR model
N=1000; %size of data
mstart=300; %data removal starts from here
mplusl=499; %data removed until here
zSize=mplusl-mstart;

%Generate AR modeled data
c = fir1(256, 0.5);
[alpha,p0] = lpc(c,order);
e = randn(1,N);
x = filter(1,alpha,e);
original_z=x(mstart:mplusl-1);

figure;
plot(x,'LineWidth',1);
xlim([0 N])
ylim([-15 15])
set(gca,'FontSize',13)
ylims = get(gca,'YLim');
hold on
plot([mstart mstart],ylims, 'r:','LineWidth',2);
hold on
plot([mplusl mplusl],ylims, 'r:','LineWidth',2);
xlabel('i');
ylabel('x_i');
text(mstart+0.5,ylims(1)+3,sprintf('Removed\n data'),'FontSize',12);
set(gca,'XTick',0:N/10:N)

%remove part of the data
y1=x(1:mstart-1);
y2=x(mplusl:end);
y=[y1 y2];
x_damaged=x;
x_damaged(mstart:mplusl-1)=0;

%restore data
iters=2*zSize;
burnIn=50;
x_restored=gibbs_sampler_audio_restoration(x_damaged,iters,burnIn);

print -deps gibbrestoredARData