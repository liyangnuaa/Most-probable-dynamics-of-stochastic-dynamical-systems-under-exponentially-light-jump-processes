clear;
clc;

% %%% alpha2=1.5  lambda2=0.5  alpha1=1.5
% lambda15=[0.5 0.25 0.1 0.075];
% W15=[0.031838 0.00817807 0.00130868 7.3614e-4];
% loglambda15=log(lambda15);
% logW15=log(W15);
% N=length(lambda15);
% xmean=sum(loglambda15)/N;
% ymean=sum(logW15)/N;
% phi15=sum((loglambda15-xmean).*(logW15-ymean))/sum((loglambda15-xmean).^2);
% logC15=ymean-phi15*xmean;
% 
% figure;
% plot(loglambda15,logW15,'ko');
% 
% %%% alpha2=2  lambda2=1  alpha1=2
% lambda221=[1 0.75 0.5 0.25];
% W221=[0.30398491 0.20450573 0.11200057 0.03974673];
% loglambda221=log(lambda221);
% logW221=log(W221);
% N=length(lambda221);
% xmean=sum(loglambda221)/N;
% ymean=sum(logW221)/N;
% phi221=sum((loglambda221-xmean).*(logW221-ymean))/sum((loglambda221-xmean).^2);
% logC221=ymean-phi15*xmean;
% 
% 
% %%% alpha2=1.5  lambda1=1  alpha1=1.5
% lambda15=[1 0.75 0.5 0.25];
% W15=[0.19907001 0.14266048 0.11566199 0.07324356];
% loglambda15151=log(lambda15);
% logW15151=log(W15);
% N=length(lambda15);
% xmean=sum(loglambda15151)/N;
% ymean=sum(logW15151)/N;
% phi15151=sum((loglambda15151-xmean).*(logW15151-ymean))/sum((loglambda15151-xmean).^2);
% logC15151=ymean-phi15151*xmean;
% 
% %%% alpha2=2  lambda1=0.5  alpha1=2
% lambda15=[0.5 0.25 0.1];
% W15=[0.0773919 0.048315 0.0355768];
% loglambda225=log(lambda15);
% logW225=log(W15);
% N=length(lambda15);
% xmean=sum(loglambda225)/N;
% ymean=sum(logW225)/N;
% phi225=sum((loglambda225-xmean).*(logW225-ymean))/sum((loglambda225-xmean).^2);
% logC225=ymean-phi225*xmean;



%%% alpha2=1.5  alpha1=1.5  lambda1=lambda2
lambda1215=[1 0.75 0.5 0.25];
W1215=[0.19907001 0.09375062 0.031838 0.00501755];
loglambda1215=log(lambda1215);
logW1215=log(W1215);
N=length(lambda1215);
xmean=sum(loglambda1215)/N;
ymean=sum(logW1215)/N;
phi1215=sum((loglambda1215-xmean).*(logW1215-ymean))/sum((loglambda1215-xmean).^2);
logC1215=ymean-phi1215*xmean;

%%% alpha1=2  alpha2=2  lambda1=lambda2
lambda2=[0.75 0.5 0.25 0.1];
W2=[0.17194459 0.0773919 0.019376535 0.003100238];
loglambda2=log(lambda2);
logW2=log(W2);
N=length(lambda2);
xmean=sum(loglambda2)/N;
ymean=sum(logW2)/N;
phi2=sum((loglambda2-xmean).*(logW2-ymean))/sum((loglambda2-xmean).^2);
logC2=ymean-phi2*xmean;

%%% alpha1=3  alpha2=3  lambda1=lambda2
lambda3=[0.5 0.25 0.1 0.075];
W3=[0.16176319 0.06435060 0.01898396 0.01293624];
loglambda3=log(lambda3);
logW3=log(W3);
N=length(lambda3);
xmean=sum(loglambda3)/N;
ymean=sum(logW3)/N;
phi3=sum((loglambda3-xmean).*(logW3-ymean))/sum((loglambda3-xmean).^2);
logC3=ymean-phi3*xmean;

%%% alpha1=4  alpha2=4  lambda1=lambda2
lambda4=[0.25 0.1 0.075 0.05];
W4=[0.10945628 0.04383735 0.03287972 0.02192079];
loglambda4=log(lambda4);
logW4=log(W4);
N=length(lambda4);
xmean=sum(loglambda4)/N;
ymean=sum(logW4)/N;
phi4=sum((loglambda4-xmean).*(logW4-ymean))/sum((loglambda4-xmean).^2);
logC4=ymean-phi4*xmean;

figure;
plot(loglambda1215,logW1215,'ko-');
hold on
plot(loglambda2,logW2,'ms-');
hold on
plot(loglambda3,logW3,'b+-');
hold on
plot(loglambda4,logW4,'gdiamond-');
hold off


%%% lambda=0.5
alpha001=[1.25 1.5 1.75 2 2.25 2.5 2.75 3 4 5];
W001=[0.01431539 0.031838 0.053791617 0.0773919 0.10081822 0.12299605 0.14341376 0.16176319 0.2182865 0.2516869];
Calpha=W001./0.5.^(4./alpha001);
logalpha001=log(alpha001);
logW001=log(W001);
logCalpha=log(Calpha);

% alpha001=[1.001 1.01 1.1 alpha001 7];
% Calpha=[5.4361 5.8974 7.947 Calpha 38.5208];
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
cfun=fit(x0,y0,f,'StartPoint',[1,1,-1,0.1]); %显示拟合函数，数据必须为列向量形式
xi=1:0.01:7;
yi=cfun(xi);
figure;
plot(alpha001',Calpha','r*',xi,yi,'b-');
title('拟合函数图形');
