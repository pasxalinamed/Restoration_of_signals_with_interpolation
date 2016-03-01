%Tests Gibbs sampling, EM, ML on an orchestral song of 1000 samples

clear;close all;
global order N zSize mstart mplusl y1 y2 ylims y B

order=40; %order of AR model
mstart=500; %data removal starts from here
mplusl=mstart+200; %data removed until here
zSize=mplusl-mstart;

%audio signal
[x,fs] = audioread('Kalonixtia_6sec.wav');
x = resample(x,1,5);
fs=fs/5;
x=x(:,1);
N=length(x);

figure;
plot(x,'LineWidth',1);
xlim([0 N])
set(gca,'FontSize',13)
xlabel('i');
ylabel('x_i');
print -deps wholeoriginal_orch_song

x=x(30000:31000);
x=x';
sound(x,fs)

N=length(x); %size of data
original_z=x(mstart:mplusl-1);

figure;
plot(x,'LineWidth',1);
xlim([0 N])
set(gca,'FontSize',13)
ylims = get(gca,'YLim');
hold on
plot([mstart mstart],ylims, 'r:','LineWidth',2);
hold on
plot([mplusl mplusl],ylims, 'r:','LineWidth',2);
xlabel('i');
ylabel('x_i');
text(mstart+0.5,ylims(1)+0.3,sprintf('Removed data'),'FontSize',12);
print -deps original_orch_song
%%
%remove part of the data
y1=x(1:mstart-1);
y2=x(mplusl:end);
% [y1 sparse(1,zSize) y2] is x_damaged

sound([y1 zeros(1,zSize) y2],fs)

%%
%restore data with Gibbs sampling
iters=3*zSize;
burnIn=50;
[x_restored,m]=gibbs_sampler_audio_restoration([y1 sparse(1,zSize) y2],iters,burnIn);

sound(x_restored,fs)

print -deps restored_song2_Gibbs

%%
%restore data with EM
iters=100;

x_restored_EM=EM_audio_restoration([y1 sparse(1,zSize) y2],iters);

sound(x_restored_EM,fs)

print -deps restored_song2_EM
%%
%restore data with ML
x_restored_ML=ML_audio_restoration([y1 sparse(1,zSize) y2]);

sound(x_restored_ML,fs)

print -deps restored_song2_ML