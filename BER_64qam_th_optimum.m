%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MCC_6G: Analytical BER over a wide range of Lambda_top
%                 By: Ahmed Hussein
%                ahussein@albany.edu
% This code can be reused under the CC BY license
% "https://creativecommons.org/licenses/by/2.0/"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Self_Notes: % For plot of figures 1 and 2, modify the step in Clip_top to be 20 and
% uncomment figs 1 and 2 but for figure 3 and calculations, return step to
% 1 
clc
clear all
close all
Mod=[64]; % QAM modulation, Can be changed according to which plot is requested 
Ps_dbm=10:0.1:45;
Ps=10.^(Ps_dbm./10); % Transmitted Electric Power in mW
No_dbm=-3;
No=10^(No_dbm/10);
Clip_top=[10:1:200]; % Upper clipping threshold for Optical power in mW
% Self_test: % Clip_top=89.21; % Upper clipping threshold for Optical power in mW


for loop=1:length(Clip_top)
BER=[];
gama_elec=[];

M=Mod(1);

lambda_top=Clip_top(loop)./sqrt(2.*Ps); % lambda arranged from high to low as power arranged from low to high

lambda_1=lambda_top(lambda_top<10); %Effect of clipping appears, as power increases and threshold is constant, the lambda decreases and ber increases
lambda_2=lambda_top(lambda_top>=10); % First range of sigma sq is corresponding and equal to zero

sigma_clip1_sq=zeros(1,length(lambda_2));
sigma_clip2_sq=Ps(1,length(lambda_2)+1:end).*(qfunc(lambda_1)-2.*(qfunc(lambda_1).^2)+qfunc(lambda_1).*lambda_1.^2-lambda_1.*1/sqrt(2*pi).*exp(-(lambda_1).^2./2));
sigma_clip_sq=[sigma_clip1_sq sigma_clip2_sq];

for ii=1:length(lambda_top)
% if lambda_top(ii) <= 10
     gama_elec(ii)=(1/2-qfunc(lambda_top(ii))).^2./((0.5*sigma_clip_sq(ii)./(Ps(ii)./log2(M)))+(0.5*No/(Ps(ii)./log2(M))));
%     
% elseif lambda_top(ii) > 10
%     
%     gama_elec(ii) = (Ps(ii)./log2(M))./No;
%     
% end
end
% figure(1)
% 
% plot(lambda_top,gama_elec,'Linewidth',2)
% hold all
lambdatop_peak1(loop)=find(gama_elec==max(gama_elec));
lambdatop_peak2(loop)=lambda_top(lambdatop_peak1(loop));
BER=(4*(sqrt(M)-1))/(sqrt(M)*log2(M))*(qfunc(sqrt((3*log2(M))/(M-1)*gama_elec))) + (4*(sqrt(M)-2))/(sqrt(M)*log2(M))*(qfunc(3*sqrt((3*log2(M))/(M-1)*gama_elec)));

P_opt_loc(loop)=find(BER==min(BER)); % Optimum Power for minimum BER (Needs more optimization with constraints on BER)
P_opt(loop)=Ps(P_opt_loc(loop));

minBER(loop)=min(BER);

L=loop-1*10;
if mod(L,10)==0
figure(2)
semilogy(Ps_dbm-No_dbm,BER,'Linewidth',2)    
hold all
end

end
P_opt_dBm=10*log10(P_opt);

figure(2)
semilogy(P_opt_dBm-No_dbm,minBER,'bo','Linewidth',2);
xlabel('E_{b}/N_{o} (dB)')
ylabel('BER')
set(gca,'FontSize',14,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')

Aclip_sigm=sqrt(lambdatop_peak2.*sqrt(2*pi)/2)

figure(3)
[Hy,Hx]=hist(Aclip_sigm,10);
bar(Hx,Hy);
xlabel('Clipping Threshold (V)')

set(gca,'FontSize',14,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')

Optimumclip_ratio=Hx(end) % Optimum Clip electrical level ratio to standard deviation
Optimumlambda_top=Optimumclip_ratio^2*2/sqrt(2*pi)

