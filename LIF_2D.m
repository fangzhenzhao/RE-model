function [phi,f1,f2,Hphi] = LIF_2D(I,phi,timestep,epsilon,K,lamda1_add_lamda2,length_u,dis_mu,area_alfa,g)

phi = NeumannBoundCond(phi);
  
Hphi = Heaviside(phi,epsilon);
DiracPhi = Delta(phi,epsilon);

[f1,f2] = Local_Avr(I,Hphi,K);

[phi_x,phi_y]=gradient(phi);
s1=sqrt(phi_x.^2+phi_y.^2);
s=s1+eps.*(s1==0);
Nx=phi_x./s;
Ny=phi_y./s;
cuaruve=div(Nx,Ny);

delta_I=4*del2(f1).*Hphi+4*del2(f2).*(1-Hphi);

enegry_term=lamda1_add_lamda2*delta_I.*DiracPhi;
length_term=length_u*DiracPhi.*cuaruve;
DistanceRegul_Term=dis_mu*DRT(phi);
area_term=area_alfa*g.*DiracPhi;

phi = phi + timestep*(enegry_term+length_term+DistanceRegul_Term+area_term);
% phi = phi + timestep*DiracPhi.*((I - f1.*Hphi - f2.*(1 - Hphi)).*(f1 - f2));

function H = Heaviside(phi,epsilon)
H = 0.5*(1+(2/pi)*atan(phi./epsilon));

function Delta_h = Delta(phi,epsilon)
Delta_h = (epsilon/pi)./(epsilon^2+phi.^2);

function [f1,f2] = Local_Avr(I,H,K)

f1 = conv2(I.*H,K,'same');
c1 = conv2(H,K,'same');
f1 = f1./c1;
f2 = conv2(I.*(1-H),K,'same');
c2 = conv2(1-H,K,'same');
f2 = f2./c2;

function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);

function f=DRT(phi)
[nx,ny]=gradient(phi);
s2=sqrt(nx.^2+ny.^2);
a=(s2>=0) & (s2<=1);
b=(s2>1);
ps=(s2-1./((s2==0)+(s2~=0).*s2)).*b+a.*(2*s2.^3-3*s2.^2+s2);
dps=((ps~=0).*ps+(ps==0))./((s2~=0).*s2+(s2==0)); 
f = div(dps.*nx - nx, dps.*ny - ny) + 4*del2(phi);  

function k=div(x,y)
[nxx,nxy]=gradient(x);
[nyx,nyy]=gradient(y);
k=nxx+nyy;