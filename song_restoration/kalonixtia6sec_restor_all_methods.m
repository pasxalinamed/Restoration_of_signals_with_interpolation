%Tests Gibbs sampling, EM, ML on an orchestral song of 6 sec.

clear;close all;
global order N zSize mstart mplusl y1 y2 ylims y B

order=40; %order of AR model
mstart=300; %data removal starts from here
mplusl=mstart+800; %data removed until here
zSize=mplusl-mstart;

%audio signal
[x_initial,fs] = audioread('Kalonixtia_6sec.wav');
x_initial=x_initial(:,1);
x_initial=x_initial';
x_resampled = resample(x_initial,1,5);
start_part=30000;
ending_part=31500;
current_x=x_resampled(1,start_part:ending_part);
fs=fs/5;
sound(x_resampled,fs)

N=length(current_x); %size of data
original_z=current_x(mstart:mplusl-1);

figure;
plot(x_resampled,'LineWidth',1);
set(gca,'FontSize',13)
ylims = get(gca,'YLim');
hold on
plot([start_part+mstart start_part+mstart],ylims, 'r:','LineWidth',2);
hold on
plot([start_part+mplusl start_part+mplusl],ylims, 'r:','LineWidth',2);
xlabel('i');
ylabel('x_i');
print -deps kalonixtiaremoveddata
%%
%remove part of the data
y1=current_x(1:mstart-1);
y2=current_x(mplusl:end);
x_damaged=[x_resampled(1:start_part) y1 sparse(1,zSize) y2 x_resampled(ending_part:end)];

sound(full(x_damaged),fs)
audiowrite('kalonixtiaDamaged.wav',full(x_damaged),fs);
%%
%restore data with Gibbs sampling
iters=3*zSize;
burnIn=500;
[x_restored,m]=gibbs_sampler_audio_restoration([y1 sparse(1,zSize) y2],iters,burnIn);

x_resampled(1,start_part:ending_part)=x_restored;
x_restoredgibb=x_resampled;

sound(x_restoredgibb,fs)
figure;
plot(x_restoredgibb)
set(gca,'FontSize',13)
ylims = get(gca,'YLim');
hold on
plot([start_part+mstart start_part+mstart],ylims, 'r:','LineWidth',2);
hold on
plot([start_part+mplusl start_part+mplusl],ylims, 'r:','LineWidth',2);
xlabel('i');
ylabel('x_i');
audiowrite('kalonixtiaGibbrest.wav',x_restoredgibb,fs);

print -deps kalonixtiarestored_song_Gibbs
%%
%restore data with EM
iters=100;
x_restored=EM_audio_restoration([y1 sparse(1,zSize) y2],iters);

x_resampled(1,start_part:ending_part)=x_restored;
x_restored_EM=x_resampled;

sound(x_restored_EM,fs)
figure;
plot(x_restored_EM)
set(gca,'FontSize',13)
ylims = get(gca,'YLim');
hold on
plot([start_part+mstart start_part+mstart],ylims, 'r:','LineWidth',2);
hold on
plot([start_part+mplusl start_part+mplusl],ylims, 'r:','LineWidth',2);
xlabel('i');
ylabel('x_i');
audiowrite('kalonixtiaEMrest.wav',x_restored_EM,fs);
print -deps kalonixtiarestored_song_EM
%%
%restore data with ML
x_restored=ML_audio_restoration([y1 sparse(1,zSize) y2]);

x_resampled(1,start_part:ending_part)=x_restored;
x_restored_ML=x_resampled;

sound(x_restored_ML,fs)
figure;
plot(x_restored_ML)
set(gca,'FontSize',13)
ylims = get(gca,'YLim');
hold on
plot([start_part+mstart start_part+mstart],ylims, 'r:','LineWidth',2);
hold on
plot([start_part+mplusl start_part+mplusl],ylims, 'r:','LineWidth',2);
xlabel('i');
ylabel('x_i');
audiowrite('kalonixtiaMLest.wav',x_restored_ML,fs);
print -deps kalonixtiarestored_song_ML
save kalonixtia