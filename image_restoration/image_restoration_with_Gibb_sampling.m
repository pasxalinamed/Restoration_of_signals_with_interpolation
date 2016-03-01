%Restores an image area using Gibbs sampling

clear;close all;
global order N zSize mstart mplusl y1 y2

order=10; %order of AR model

x = imread('man.jpg');
imshow(x)
x = rgb2gray(x);
x = im2double(x);
figure;
imshow(x);
display('Please select area to restore')
rect= getrect;

x_coord=uint64(abs(rect(1)));
y_coord=uint64(abs(rect(2)));
width=uint64(abs(rect(3)));
height =uint64(abs(rect(4)));

x_damaged=x;
x_damaged(y_coord:y_coord+height, x_coord:x_coord+width)=1;
figure;
imshow(x_damaged);
x_new=x_damaged;

figure
for i=y_coord:y_coord+height
    
    to_restore=x_damaged(i,:);
    mstart=x_coord; %data removal starts from here
    mplusl=x_coord+width+1; %data removed until here
    zSize=abs(mplusl-mstart);
    zSize=double(zSize);
    N=length(to_restore); %size of data
    
    %remove part of the data
    y1=to_restore(1:mstart-1);
    y2=to_restore(mplusl:end);
    
    %restore data with Gibbs sampling
    iters=zSize*3;
    burnIn=30;
    
    x_restored=gibbs_sampler_image_restoration(to_restore,iters,burnIn);
    
    x_new(i,:)=x_restored;
    imshow(x_new)
end

imwrite(x_new,'gibbrestored_image.jpg')
print -deps gibbrestored_man_image