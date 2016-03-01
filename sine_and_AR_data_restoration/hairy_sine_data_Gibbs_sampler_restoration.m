% Tests Gibbs sampler algorithm on hairy sine data 

clear;close all;
global order N x_damaged zSize y mstart mplusl y1 y2 ylims

order=80; %order of AR models
N=1000; %size of data
mstart=300; %data removal starts from here
mplusl=499; %data removed until here
zSize=mplusl-mstart;

% Hairy sine wave data
vect=0:1:N-1;% Time Samples
f=5;% Input Signal Frequency
fs=1000;% Sampling Frequency
x=sin(2*pi*f/fs*vect);% Generate Sine Wave  
stand_dev=0.2;
e= random('norm',0, stand_dev,1,N);
x=x+e;

original_z=x(mstart:mplusl-1);

figure;
plot(x,'LineWidth',1);
xlim([0 N])
ylim([-2 2])
set(gca,'FontSize',13)
ylims = get(gca,'YLim');
hold on
plot([mstart mstart],ylims, 'r:','LineWidth',2);
hold on
plot([mplusl mplusl],ylims, 'r:','LineWidth',2);
xlabel('i');
ylabel('x_i');
text(mstart+0.5,ylims(2)-ylims(2)/2,sprintf('Removed\n data'),'FontSize',12);
set(gca,'XTick',0:N/10:N)

%remove part of the data
y1=x(1:mstart-1);
y2=x(mplusl:end);
y=[y1 y2];
x_damaged=x;
x_damaged(mstart:mplusl-1)=0;

%restore data
iters=100;
burnIn=10;
x_restored=gibbs_sampler_audio_restoration(x_damaged,iters,burnIn);