clc;
clear all;
close all;
I=imread('mean_var.bmp');
I=double(I(:,:,1));

[m,n]=size(I);
for i=1:m
    for j=20:n
        if (I(i,j)==114)
            I(i,j)=I(i,j)+(j-20)/1.5;
        end 
    end
end
figure;
imagesc(I,[0,255]);colormap(gray);
img=uint8(I);
imwrite(img,'inhomenity.bmp');