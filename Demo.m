%*********************ͼ��ƥ������ģ�ͳ���****************
clear all;
close all;
clc;

Img = imread('twoObj.bmp');
% Img=rgb2gray(Img);
Img = double(Img(:,:,1));


%*****************ר����ԡ�mean_var.bmp��ͼ���ƽ��*******************
%  K_temp= fspecial('gaussian',5,.8);
% smooth=conv2(Img,K_temp,'same');
% Img=smooth;
%**************************************************************

sigma =2;% �˲����Ĵ�С�����˷ָ����ͼ��.

K = fspecial('gaussian',2*round(2*sigma)+1,sigma);

%***********************��ʼˮƽ��ѡȡ****************************
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

%*****************��������************************************************
timestep =1;
epsilon =1;
area_alfa=0;
dis_mu =.04;                  %timestep*dis_mu<h^2/4=1/4;����CFL����
% dis_mu =0;
length_u=1;
lamda1_add_lamda2=1;%������С��������ʵ�ѡȡ
%********************************************************************
time = cputime;
%***********************��ʼѭ������*********************************************
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
%***************���յ����ͼ����*****************************************
fitting_img=f1.*phi+f2.*(1-phi);
%********************************************************************

%***************ˮƽ����ά��ʾ***************************************
figure(2);
mesh(-phi);
%*****************************************************************
 

%********************����Ϊ�ָ�ֳ���******************************
phi=phi./255;
mask=phi.*(Img+1);%��ȥԭͼ��Ϊ��Ĳ���
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
sigma=.5;      %  ��׼��������ƽ����ͼС���⣩������ģ��
M=fspecial('gaussian',round(2*sigma)*2+1,sigma); % Gaussian kernel
mask=conv2(mask,M,'same');  %ͼ�����
figure(3);imagesc(mask,[0,255]);colormap('gray');hold on;axis off;

figure(4);
imagesc(fitting_img,[0,255]);colormap(gray);axis off;
%**************************************************************************



