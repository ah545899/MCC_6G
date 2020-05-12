%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MCC_6G: Analytical BER against optimum Lambda_top and SNR
%                 By: Ahmed Hussein
%                ahussein@albany.edu
% This code can be reused under the CC BY license
% "https://creativecommons.org/licenses/by/2.0/"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These files are saved and obtained from running BER_64qam_th_optimum.m
% for multiple modulation orders

A1=load('qam64_minBER.mat','lambdatop_peak2','minBER','P_opt_dBm','No_dbm')
A2=load('qam32_minBER.mat','lambdatop_peak2','minBER','P_opt_dBm','No_dbm')
A3=load('qam16_minBER.mat','lambdatop_peak2','minBER','P_opt_dBm','No_dbm')

figure(1)
semilogy(A1.lambdatop_peak2,A1.minBER,'o','linewidth',2)
hold all
semilogy(A2.lambdatop_peak2,A2.minBER,'s','linewidth',2)
hold all
semilogy(A3.lambdatop_peak2,A3.minBER,'*','linewidth',2)
grid
axis([1.6 3.5 10^-10 1])
legend('64-QAM','32-QAM','16-QAM')
xlabel('\lambda_{top}')
ylabel('BER')
set(gca,'FontSize',14,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')

figure(2)
semilogy(A1.P_opt_dBm-A1.No_dbm,A1.minBER,'o','linewidth',2)
hold all
semilogy(A2.P_opt_dBm-A2.No_dbm,A2.minBER,'s','linewidth',2)
hold all
semilogy(A3.P_opt_dBm-A3.No_dbm,A3.minBER,'*','linewidth',2)
grid
axis([15 35 10^-10 1])
legend('64-QAM','32-QAM','16-QAM')
xlabel('E_{b}/N_{o} (dB)')
ylabel('BER')
set(gca,'FontSize',14,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')