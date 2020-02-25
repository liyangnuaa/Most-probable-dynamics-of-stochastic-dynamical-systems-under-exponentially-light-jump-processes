clear;
clc;

%%% alpha=1.5
lambda15=[0.5 0.25 0.1 0.075];
W15=[4.201327 1.066264 0.1710587 0.09623374];
% lambda15=[0.75 0.5 0.25 0.1 0.075 0.05 0.025 0.01];
% W15=[9.160059 4.201327 1.066264 0.1710587 0.09623374 0.04277391 0.01069384 0.0017110244];
loglambda15=log(lambda15);
logW15=log(W15);
N=length(lambda15);
xmean=sum(loglambda15)/N;
ymean=sum(logW15)/N;
phi15=sum((loglambda15-xmean).*(logW15-ymean))/sum((loglambda15-xmean).^2);
logC15=ymean-phi15*xmean;

%%% alpha=2
lambda2=[0.1 0.075 0.05 0.025];
W2=[0.8128709 0.528309 0.2877052 0.1017468];
loglambda2=log(lambda2);
logW2=log(W2);
N=length(lambda2);
xmean=sum(loglambda2)/N;
ymean=sum(logW2)/N;
phi2=sum((loglambda2-xmean).*(logW2-ymean))/sum((loglambda2-xmean).^2);
logC2=ymean-phi2*xmean;

%%% alpha=2.5
lambda25=[0.075 0.05 0.025 0.01];
W25=[1.384221 0.85191433 0.3711322 0.12363477];
% lambda25=[0.1 0.075 0.05 0.025 0.01 0.0075 0.005 0.0025];
% W25=[1.9521955 1.384221 0.85191433 0.3711322 0.12363477 0.087544659 0.053818651 0.02342643];
loglambda25=log(lambda25);
logW25=log(W25);
N=length(lambda25);
xmean=sum(loglambda25)/N;
ymean=sum(logW25)/N;
phi25=sum((loglambda25-xmean).*(logW25-ymean))/sum((loglambda25-xmean).^2);
logC25=ymean-phi25*xmean;

%%% alpha=3
lambda3=[0.05 0.025 0.01 0.0075];
W3=[1.706197 0.8545498 0.3420913 0.2565956];
loglambda3=log(lambda3);
logW3=log(W3);
N=length(lambda3);
xmean=sum(loglambda3)/N;
ymean=sum(logW3)/N;
phi3=sum((loglambda3-xmean).*(logW3-ymean))/sum((loglambda3-xmean).^2);
logC3=ymean-phi3*xmean;

%%% alpha=3.5
lambda35=[0.01 0.0075 0.005 0.0025];
W35=[0.6965322 0.54443585 0.38468282 0.21240368];
loglambda35=log(lambda35);
logW35=log(W35);
N=length(lambda35);
xmean=sum(loglambda35)/N;
ymean=sum(logW35)/N;
phi35=sum((loglambda35-xmean).*(logW35-ymean))/sum((loglambda35-xmean).^2);
logC35=ymean-phi35*xmean;

figure;
plot(loglambda15,logW15,'ko-');
hold on
plot(loglambda2,logW2,'ms-');
hold on
plot(loglambda25,logW25,'r*-');
hold on
plot(loglambda3,logW3,'b+-');
hold on
plot(loglambda35,logW35,'gdiamond-');
hold off

%%% lambda=0.01
alpha001=[1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 4 5];
W001=[1.8192583e-4 0.0017110244 0.008164133 0.0257421 0.06191863 0.12363477 0.21605158 0.34209129 0.50255539 1.1756906 2.4078445];
Calpha=W001./0.01.^(3./alpha001);
logalpha001=log(alpha001);
logW001=log(W001);
logCalpha=log(Calpha);

alpha001=[1.001 1.01 1.1 alpha001 7];
Calpha=[5.4361 5.8974 7.947 Calpha 38.5208];
figure;
plot(alpha001,Calpha,'ko');

% alpha1=alpha001(2);
% alpha2=alpha001(7);
% alpha3=alpha001(11);
% alpha4=alpha001(15);
% C1=Calpha(2);
% C2=Calpha(7);
% C3=Calpha(11);
% C4=Calpha(15);
% % fr=@(r)[(C2-C1)/(C3-C1)*((alpha2+r(1))^(r(2))*(alpha1+r(1))^(r(2))-(alpha2+r(1))^(r(2))*(alpha3+r(1))^(r(2)))-...
% %     ((alpha3+r(1))^(r(2))*(alpha1+r(1))^(r(2))-(alpha2+r(1))^(r(2))*(alpha3+r(1))^(r(2)));
% %     (C4-C1)/(C3-C1)*((alpha4+r(1))^(r(2))*(alpha1+r(1))^(r(2))-(alpha4+r(1))^(r(2))*(alpha3+r(1))^(r(2)))-...
% %     ((alpha3+r(1))^(r(2))*(alpha1+r(1))^(r(2))-(alpha4+r(1))^(r(2))*(alpha3+r(1))^(r(2)))];
% fr=@(r)[(C2-C1)/(C3-C1)*((alpha3+r(1))^(-r(2))-(alpha1+r(1))^(-r(2)))-((alpha2+r(1))^(-r(2))-(alpha1+r(1))^(-r(2)));
%     (C4-C1)/(C3-C1)*((alpha3+r(1))^(-r(2))-(alpha1+r(1))^(-r(2)))-((alpha4+r(1))^(-r(2))-(alpha1+r(1))^(-r(2)))];
% options = optimset('TolFun',1e-10);
% r=fsolve(fr,[2;3],options);
% k1=(C3-C1)/((alpha3+r(1))^(-r(2))-(alpha1+r(1))^(-r(2)));
% k2=C2-k1*(alpha2+r(1))^(-r(2));
% 
% x=linspace(1,7,1000);
% y=k1./(x+r(1)).^r(2)+k2;
% figure;
% plot(x,y);

x0=alpha001';
y0=Calpha';
syms t
% f=fittype('k1/(t+r1)^r2+k2','independent','t','coefficients',{'r1','r2','k1','k2'});
% cfun=fit(x0,y0,f,'StartPoint',[3,4,-1e5,38]); %显示拟合函数，数据必须为列向量形式
f=fittype('k1*exp(-r1*t^r2)+k2','independent','t','coefficients',{'r1','r2','k1','k2'});
cfun=fit(x0,y0,f,'StartPoint',[1,1,-1e2,38]); %显示拟合函数，数据必须为列向量形式
xi=1:0.01:7;
yi=cfun(xi);
figure;
plot(alpha001',Calpha','r*',xi,yi,'b-');
title('拟合函数图形');
