%*********************图像匹配拉氏模型程序****************
%by  重庆大学  山金孝 data：2011-8-25        update：2011—11—06

clear all;
close all;
clc;

Img = imread('twoObj.bmp');
% Img=rgb2gray(Img);
Img = double(Img(:,:,1));


%*****************专门针对‘mean_var.bmp’图像的平滑*******************
%  K_temp= fspecial('gaussian',5,.8);
% smooth=conv2(Img,K_temp,'same');
% Img=smooth;
%**************************************************************

sigma =2;% 此参数的大小决定了分割何种图像.

K = fspecial('gaussian',2*round(2*sigma)+1,sigma);

%***********************初始水平集选取****************************
figure(1);
imagesc(Img,[0,255]);colormap('gray');
hold on;title('choose Initial contour');
axis off;
MC=roipoly;
c0=2;
ILS=c0*(0.5-MC);
phi = ILS;
[c, h] = contour(phi, [0 0], 'r');
%*************************************************************************

%*****************参数设置************************************************
timestep =1;
epsilon =1;
area_alfa=0;
dis_mu =.04;                  %timestep*dis_mu<h^2/4=1/4;满足CFL条件
% dis_mu =0;
length_u=1;
lamda1_add_lamda2=1;%正负大小根据情况适当选取
%********************************************************************
time = cputime;
%***********************开始循环迭代*********************************************
for n = 1:200     
      [phi,f1,f2,Hphi]= Evolution(Img,phi,timestep,epsilon,K,lamda1_add_lamda2,length_u,dis_mu,area_alfa);
%       phi = conv2(phi,K_phi,'same');      
      if mod(n,20) == 0
      pause(0.0001);
      imagesc(Img,[0 255]);colormap(gray);axis off;
      hold on;contour(phi,[0 0],'r','LineWidth',2);
      iterNum=[num2str(n), ' iterations'];        
%       title(iterNum);         
      hold off;
      end
end
%************************************************************************
totaltime = cputime - time
%***************最终的拟合图像域*****************************************
fitting_img=f1.*phi+f2.*(1-phi);
%********************************************************************

%***************水平集三维显示***************************************
figure(2);
mesh(-phi);
%*****************************************************************
 

%********************下面为分割部分程序******************************
phi=phi./255;
mask=phi.*(Img+1);%削去原图中为零的部分
[r,c]=size(mask);
for i=1:r
    for j=1:c
        if mask(i,j)<=0
            mask(i,j)=255;
        else
            mask(i,j)=0;
        end
    end
end
%mask=(mask~=0);
sigma=.5;      %  标准差愈大俞平滑（图小而尖），即俞模糊
M=fspecial('gaussian',round(2*sigma)*2+1,sigma); % Gaussian kernel
mask=conv2(mask,M,'same');  %图像卷积
figure(3);imagesc(mask,[0,255]);colormap('gray');hold on;axis off;

figure(4);
imagesc(fitting_img,[0,255]);colormap(gray);axis off;
%**************************************************************************



