%Tests Gibbs sampling, EM, ML on an acapella song of 1000 samples

clear;close all;
global order N zSize mstart mplusl y1 y2 ylims y B

order=40; %order of AR model
mstart=400; %data removal starts from here
mplusl=mstart+200; %data removed until here
zSize=mplusl-mstart;

%audio signal
[x,fs] = audioread('kalanta_4sec.wav');
x=x';
x = resample(x,1,5);
fs=fs/5;
x=x(4000:5000);
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
text(mstart+0.5,ylims(1)+0.3-0.2,sprintf('Removed data'),'FontSize',12);
print -deps original_acap_song
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

print -deps restored_song_Gibbs

%%
%restore data with EM
iters=100;
x_restored_EM=EM_audio_restoration([y1 sparse(1,zSize) y2],iters);
sound(x_restored_EM,fs)

print -deps restored_song_EM
%%
%restore data with ML
x_restored_ML=ML_audio_restoration([y1 sparse(1,zSize) y2])
sound(x_restored_ML,fs)

print -deps restored_song_ML