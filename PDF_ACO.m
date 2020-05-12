%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                MCC_6G: Plot ACO-OFDM
%                 By: Ahmed Hussein
%                ahussein@albany.edu
% This code can be reused under the CC BY license
% "https://creativecommons.org/licenses/by/2.0/"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigmaA=1;
w=0:1/(100*2*pi):2*pi;
f=1/(sqrt(2*pi)*sigmaA)*exp(-w.^2./(2*sigmaA^2))+1/2*dirac(w);
plot(w,f,'linewidth',2)
xlabel('w (rad/sec)')
ylabel('PDF')
set(gca,'FontSize',14,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')
