clear;
clc;

global r1 r2 alpha lambda Ch gamma theta S0
r1=0;
r2=1;
alpha=1.5;
lambda=0.1;
Ch=46.8;
gamma=0.61;
theta=5.67e-8;
S0=1368;

fxy=@(x)(S0*(0.5+0.2*tanh((x-265)/10))/4-gamma*theta*x^4)/Ch;
xnode1=fsolve(fxy,225);
xnode2=fsolve(fxy,300);
xsad=fsolve(fxy,260);
xnode=xnode1;

% syms x r1 r2 alpha lambda Ch gamma theta S0
% y=(S0*(0.5+0.2*tanh((x-265)/10))/4-gamma*theta*x^4)/Ch;
% z=diff(y,x);

R=0.05;
x=xnode;
fder=-((S0*(tanh(x/10 - 53/2)^2/50 - 1/50))/4 + 4*gamma*theta*x^3)/Ch;
fB=@(x)r2^2*x.^2.*exp(-lambda*abs(x).^alpha);

B=zeros(2,2);
B(1,1)=fder;
B(1,2)=r1^2+quadgk(fB,-inf,inf);
B(2,2)=-fder;
[Bv,Be]=eig(B);
Bv1=Bv(1,2);
Bv2=Bv(2,2);
M=real(Bv2/Bv1);

C=fder;
D=B(1,2);
IR2=D-r1^2;
Z1=-2*C*r1^2/D^2;
Z2=-2*C*IR2/D^2;

tf=300;
h=0.005;
nT=floor(tf/h);
X=zeros(5,nT);
X(:,1)=[xnode+R;M*R;1/2*(Z1+Z2)*R^2;1/2*Z1*R^2;1/2*Z2*R^2];
% X(:,1)=[xnode-R;-M*R;1/2*(Z1+Z2)*R^2;1/2*Z1*R^2;1/2*Z2*R^2];

for j=1:nT-1
    t0=0;
    x1=X(:,j);
    x2=rk4(t0,h,x1);
    X(:,j+1)=x2;
end
Y=X;
% Xd=X(1,2:end)-X(1,1:end-1);
% % I=find(Xd<0,1,'first');
% I=find(Xd>0,1,'first');
% if ~isempty(I)
%     X=X(:,1:I);
% end
X0=[xnode;0;0;0;0];
X=[X0 X];

figure;
plot(X(1,:),X(2,:));
% n0=18500;
% figure;
% plot(X(1,1:n0),X(2,1:n0));
% figure;
% plot(X(1,1:n0),X(5,1:n0)+0.051468888507113);
% figure;
% plot([xnode1 xsad],[0 0]);
% figure;
% plot([xnode2 xsad],[0 0]);
% figure;
% plot(X(1,:),X(3,:));
% 
% figure;
% plot(X(1,:),X(4,:));

figure;
plot(X(1,:),X(3,:));

%%% È·¶¨ÐÔÊÆÚå
xdet=linspace(200,310,10000);
UT=(-1/4*S0*(0.5*xdet+2*log(cosh((xdet-265)/10)))+1/5*gamma*theta*xdet.^5)/Ch;
figure;
plot(xdet,UT);
figure;
plot(xdet,(UT-min(UT))/67,'r-');

% x=linspace(200,300,10000);
% y0=zeros(size(x));
% UT=(S0*(0.5+0.2*tanh((x-265)/10))/4-gamma*theta*x.^4)/Ch;
% figure;
% plot(x,UT);
% hold on
% plot(x,y0);
% grid on
% hold off


function xout=rk4(t0,h,x0)
k1=h*fun(t0,x0);
k2=h*fun(t0+h/2,x0+0.5*k1);
k3=h*fun(t0+h/2,x0+0.5*k2);
k4=h*fun(t0+h,x0+k3);
xout=x0+(k1+2*k2+2*k3+k4)/6;


function y=fun(~,x)

global r1 r2 alpha lambda Ch gamma theta S0

fx=@(t)r2*t.*exp(r2*x(2)*t-lambda*abs(t).^alpha);
f2=@(t)(r2*t*x(2).*exp(r2*t*x(2))-exp(r2*t*x(2))+1).*exp(-lambda*abs(t).^alpha);
y=zeros(size(x));
y(1)=(S0*(0.5+0.2*tanh((x(1)-265)/10))/4-gamma*theta*x(1)^4)/Ch+r1^2*x(2)+quadgk(fx,-inf,inf);
y(2)=-(-((S0*(tanh(x(1)/10 - 53/2)^2/50 - 1/50))/4 + 4*gamma*theta*x(1)^3)/Ch)*x(2);
y(3)=x(2)*y(1);
y(4)=1/2*r1^2*x(2)^2;
y(5)=quadgk(f2,-inf,inf);
