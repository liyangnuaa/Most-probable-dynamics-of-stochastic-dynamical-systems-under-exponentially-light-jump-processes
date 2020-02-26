clear;
clc;

xmin=-1.8;
xmax=1.8;
ymin=-1.5;
ymax=1.5;

global r1 r2 alpha1 lambda1 alpha2 lambda2 gamma miu
r1=0;
r2=1;
alpha1=1.5;
lambda1=0.1;
alpha2=1.5;
lambda2=0.1;
miu=1;
gamma=1;

fxy=@(x)[x(1)-x(1)^3-gamma*x(1)*x(2)^2;-miu*(1+x(1)^2)*x(2)];
% xnode1=fsolve(fxy,225);
% xnode2=fsolve(fxy,300);
% xsad=fsolve(fxy,260);
xnode1=[-1;0];
xnode2=[1;0];
xsad=[0;0];
xnode=xnode1;

A=[1-3*xnode(1)^2-gamma*xnode(2)^2, -2*gamma*xnode(1)*xnode(2);
    -2*miu*xnode(1)*xnode(2), -miu*(1+xnode(1)^2)];
f10=@(x)exp(-lambda1*abs(x).^alpha1);
f11=@(x)x.*exp(-lambda1*abs(x).^alpha1);
f12=@(x)x.^2.*exp(-lambda1*abs(x).^alpha1);
f20=@(x)exp(-lambda2*abs(x).^alpha2);
f21=@(x)x.*exp(-lambda2*abs(x).^alpha2);
f22=@(x)x.^2.*exp(-lambda2*abs(x).^alpha2);
C=zeros(2,2);
C(1,1)=r1^2+r2^2*quadgk(f12,-inf,inf)*quadgk(f20,-inf,inf);
C(1,2)=r2^2*quadgk(f11,-inf,inf)*quadgk(f21,-inf,inf);
C(2,1)=C(1,2);
C(2,2)=r1^2+r2^2*quadgk(f10,-inf,inf)*quadgk(f22,-inf,inf);

B=zeros(4,4);
B(1:2,1:2)=A;
B(1:2,3:4)=C;
B(3:4,3:4)=-A';
[Bv,Be]=eig(B);
Bv1=Bv(1:2,3:4);
Bv2=Bv(3:4,3:4);
M=real(Bv2/Bv1);

A0=[2*A(1,1) 2*A(1,2) 0;A(2,1) A(1,1)+A(2,2) A(1,2);0 2*A(2,1) 2*A(2,2)];
B0=-[C(1,1);C(1,2);C(2,2)];
s=A0\B0;
Z=inv([s(1) s(2);s(2) s(3)]);

Nphi=1000;            % 环上划分精度
Sphi=zeros(1,Nphi);    % 所有的phi中S的最小值
phi=linspace(0,2*pi,Nphi);
% phi1=linspace(0,0.01,Nphi/2);
% phi2=linspace(2*pi-0.01,2*pi,Nphi/2);
% phi=[phi1 phi2];

R=0.01;
tf=1000;
h=0.005;
nT=floor(tf/h);

Np=zeros(1,Nphi);
xlamS=zeros(5,Nphi);
xlamS(1:2,:)=[xnode(1)+R*cos(phi);xnode(2)+R*sin(phi)];
% xlamS(3:4,:)=0;
xlamS(3:4,:)=M*[R*cos(phi);R*sin(phi)];
xlamS(5,:)=1/2*(Z(1,1)*(xlamS(1,:)-xnode(1)).^2+Z(2,2)*(xlamS(2,:)-xnode(2)).^2+2*Z(1,2)*(xlamS(1,:)-xnode(1)).*(xlamS(2,:)-xnode(2)));
% x1=zeros(Nphi,nT);
% x2=zeros(Nphi,nT);
% x3=zeros(Nphi,nT);
% x4=zeros(Nphi,nT);
% x5=zeros(Nphi,nT);
pos=1:1:Nphi;
delta=0;

for i=1:Nphi
    t0=0;
    X0=xlamS(:,i);
%     x1(i,1)=xnode(1);
%     x2(i,1)=xnode(2);
%     x3(i,1)=0;
%     x4(i,1)=0;
%     x5(i,1)=0;
%     x1(i,2)=X0(1);
%     x2(i,2)=X0(2);
%     x3(i,2)=X0(3);
%     x4(i,2)=X0(4);
%     x5(i,2)=X0(5);
    
    for j=1:nT-2
        X1=rk4(t0,h,X0);
        if (X1(1)>0)||(X1(1)<xmin)||(X1(2)>ymax)||(X1(2)<ymin)
            Np(i)=j+1;
            break;
        end
        
        X0=X1;
%         x1(i,j+2)=X0(1);
%         x2(i,j+2)=X0(2);
%         x3(i,j+2)=X0(3);
%         x4(i,j+2)=X0(4);
%         x5(i,j+2)=X0(5);
    end
    Sphi(i)=X1(5);
end

I=Sphi>100;
if ~isempty(I)
    Sphi(I)=100;
end

t=(0:nT)*h;
Sm=min(Sphi);

figure;
plot(phi,Sphi);
% axis([0-0.1 2*pi+0.1 0 0.07])
axis([0-0.1 1.2 0 0.004])

figure;
plot([-0.1 1.2],[0.004 0.004]);
hold on
plot([1.2 1.2],[0 0.004]);
hold off


function xout=rk4(t0,h,x0)
k1=h*fun(t0,x0);
k2=h*fun(t0+h/2,x0+0.5*k1);
k3=h*fun(t0+h/2,x0+0.5*k2);
k4=h*fun(t0+h,x0+k3);
xout=x0+(k1+2*k2+2*k3+k4)/6;


function y=fun(~,x)

global r1 r2 alpha1 lambda1 alpha2 lambda2 gamma miu

f10=@(t)exp(r2*x(3)*t-lambda1*abs(t).^alpha1);
f11=@(t)t.*exp(r2*x(3)*t-lambda1*abs(t).^alpha1);
f20=@(t)exp(r2*x(4)*t-lambda2*abs(t).^alpha2);
f21=@(t)t.*exp(r2*x(4)*t-lambda2*abs(t).^alpha2);

y=zeros(size(x));
y(1)=x(1)-x(1)^3-gamma*x(1)*x(2)^2+r1^2*x(3)+r2*quadgk(f11,-inf,inf)*quadgk(f20,-inf,inf);
y(2)=-miu*(1+x(1)^2)*x(2)+r1^2*x(4)+r2*quadgk(f10,-inf,inf)*quadgk(f21,-inf,inf);
y(3)=-(1-3*x(1)^2-gamma*x(2)^2)*x(3)+2*miu*x(1)*x(2)*x(4);
y(4)=2*gamma*x(1)*x(2)*x(3)+miu*(1+x(1)^2)*x(4);
y(5)=x(3)*y(1)+x(4)*y(2);
