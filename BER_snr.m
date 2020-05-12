%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                MCC_6G: Experimental BER
%                 By: Ahmed Hussein
%                ahussein@albany.edu
% This code can be reused under the CC BY license
% "https://creativecommons.org/licenses/by/2.0/"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('SNR.mat','snr_x','snr_xdb'); % Measured SNR Values
MM=[64,32,16]; % Calculating BER based on EVM relation to SNR
for loop=1:3
M=MM(loop);
L=sqrt(M);
EVM=sqrt(1./snr_x);
BER=2*(1-1/L)/log2(L)*qfunc(sqrt((3*log2(L)./(L^2-1)).*(2./(EVM.^2*log2(M)))));
semilogy(snr_xdb,BER,'linewidth',2)
hold all
end
xlabel('SNR (dB)')
ylabel('BER')
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'marker','text'),'FontSize',14)
