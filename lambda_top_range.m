%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       MCC_6G: Plot sigma_clip_square VS Lambda_top
%                 By: Ahmed Hussein
%                ahussein@albany.edu
% This code can be reused under the CC BY license
% "https://creativecommons.org/licenses/by/2.0/"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
x=0:0.001:10;
Z=x.*(1/sqrt(2*pi)*exp(-x.^2/2));
plot(x,Z,'linewidth',2)
hold all
Y=-2*(qfunc(x)).^2+qfunc(x).*(x.^2)+qfunc(x);
plot(x,Y,'linewidth',2)
xlabel('\lambda_{top}')

figure,
plot(x,Y-Z,'linewidth',2)
M1=find(Y-Z>0.0011);
lambda1=x(M1(end));
M2=find(Y-Z>1.9343e-08);
lambda2=x(M2(end));
M3=find(Y-Z>1.4529*10^-25);
lambda3=x(M3(end));

xlabel('\lambda_{top}')
ylabel('\sigma_{clip}^2')

set(gca,'FontSize',14,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')
%Y(M2(end))-Z(M2(end)) % check the difference value by using M=find(Y>Z) and putting target lambda as the end of
% the array example to calculate limit for lambda=5, so x=0:0.001:5