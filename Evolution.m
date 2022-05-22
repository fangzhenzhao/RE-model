function [phi,f1,f2,Hphi] = sjx_evolution(I,phi,timestep,epsilon,K,lamda1_add_lamda2,length_u,dis_mu,area_alfa)

%****************************各个函数的计算********************
phi = NeumannBoundCond(phi);  
Hphi = Heaviside(phi,epsilon);
DiracPhi = Delta(phi,epsilon);
[f1,f2] = Local_Avr(I,Hphi,K);
%*****************************************************************

%****************计算边缘停止函数****************
fitting=f1.*Hphi+f2.*(1-Hphi);
[fitting_nx,fitting_ny]=gradient(fitting);
mo=sqrt(fitting_nx.^2+fitting_ny.^2);
g=1./(1+mo.^3);
%*************************************************

%****************长度约束项（水平集梯度）计算*************
[phi_x,phi_y]=gradient(phi);
s1=sqrt(phi_x.^2+phi_y.^2);
s=s1+eps.*(s1==0);
Nx=phi_x./s;
Ny=phi_y./s;
cuaruve=div(Nx,Ny);
%**************************************************************

%*********************数据力*******************************
delta_I=4*del2(f1).*Hphi+4*del2(f2).*(1-Hphi);
%********************************************************

%******************演化方程中各项计算*********************************
enegry_term=lamda1_add_lamda2*delta_I.*DiracPhi;
length_term=length_u*DiracPhi.*cuaruve;
DistanceRegul_Term=dis_mu*DRT(phi);
area_term=area_alfa*g.*DiracPhi;
%***********************************************************************


%***********************最终迭代方程**************************************
phi = phi + timestep*(enegry_term+length_term+DistanceRegul_Term+area_term);
% phi = phi + timestep*DiracPhi.*((I - f1.*Hphi - f2.*(1 - Hphi)).*(f1 - f2));
%**************************************************************************



%***********************Heavside函数和Driac函数计算***********************
function H = Heaviside(phi,epsilon)
H = 0.5*(1+(2/pi)*atan(phi./epsilon));

function Delta_h = Delta(phi,epsilon)
Delta_h = (epsilon/pi)./(epsilon^2+phi.^2);
%*********************************************************************


%*********************局部拟合函数f1，f2的计算***********************
function [f1,f2] = Local_Avr(I,H,K)
f1 = conv2(I.*H,K,'same');
c1 = conv2(H,K,'same');
f1 = f1./c1;
f2 = conv2(I.*(1-H),K,'same');
c2 = conv2(1-H,K,'same');
f2 = f2./c2;
%********************************************************************

% **************计算 Neumann边界********************************
function g = NeumannBoundCond(f)
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);
%**************************************************************

%****************************计算水平集函数距离规则项****************
function f=DRT(phi)
[nx,ny]=gradient(phi);
s2=sqrt(nx.^2+ny.^2);
a=(s2>=0) & (s2<=1);
b=(s2>1);
ps=(s2-1./((s2==0)+(s2~=0).*s2)).*b+a.*(2*s2.^3-3*s2.^2+s2);
dps=((ps~=0).*ps+(ps==0))./((s2~=0).*s2+(s2==0)); 
f = div(dps.*nx - nx, dps.*ny - ny) + 4*del2(phi);  
%****************************************************************

%**********************散度计算***********************************
function k=div(x,y)
[nxx,nxy]=gradient(x);
[nyx,nyy]=gradient(y);
k=nxx+nyy;
%******************************************************************












