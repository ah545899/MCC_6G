%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MCC_6G: Simulated BER for MCC under different interference scenarios
%                 By: Ahmed Hussein
%                ahussein@albany.edu
% This code can be reused under the CC BY license
% "https://creativecommons.org/licenses/by/2.0/"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interference and noise effect on BER
% Initial investigation of interference effect on BER curves VS SNR
% uses U-OFDM with positive samples in positive cycle and -ve samples in
% -ve cycle
% Distortion from the fundamental harmonic with proper and improper design


clc
close all 
clear all

ma =6;
Ma = 2^ma; % modulation order
Nffta = 32;
T=1; % Can't be equal T_pwm in this code (maximum value half T_pwm)
T_pwm=2;
deltaf=1/T;
fs=Nffta*deltaf;
L=100000; % number of PWM cycles in a frame (gets 100000 bits by 16-QAM)
t=0:1/fs:T_pwm-1/fs; % because unipolar DCO is used
A=2; % Amplitude of square wave where signal is overlapped (integrated/superimposed)
% fund=0.6366; % based on spectrum analysis of a square wave in another program (obtained from graph based on A=2)
k=-Nffta/2:Nffta/2-1; % OFDM subcarriers
TX=[];
TX2pos=[];
TX2negs=[];
D_C=1/2; % pick up values (0.1 to 0.9) if T_pwm is 10 times T
% Change duty cycle based on plot
C_1=D_C*T_pwm/T;
C_2=T_pwm/T;
Bits=[];
SYMB1A=[];
Sign=[];
snr=15:1:45;
power_factor=1/4; % P_i/Pn
sinr=snr+10*log10(1+power_factor);
% har=0:1:12;
% tot_hmnc=fund./(10.^(har/20));
% num_hmncs=3;
% sngl_hmnc=tot_hmnc/num_hmncs;
err=zeros(1,length(snr));
err_prp=zeros(1,length(snr));
err_imprp_distort=zeros(1,length(snr));
err_imprp_25dc=zeros(1,length(snr));
err_imprp_10dc=zeros(1,length(snr));
err_distort=zeros(1,length(snr));
% err_imprp_one=zeros(1,length(snr));
Rpos=A*ones(1,Nffta);

for l=1:L      
        bit_streama=randi(Ma-1,Nffta/2-1,1); % Conventional OFDM
        %Bits=[Bits bit_streama];
        symb1a=qammod(bit_streama,Ma);
        %SYMB1A=[SYMB1A symb1a];        
        tx_syma = [0 transpose(symb1a) 0 (fliplr(conj(transpose(symb1a))))];
        ya = ifft(tx_syma);
        yapos=zeros(1,length(ya));
        ynegs=zeros(1,length(ya));

        yapos(ya>0)=ya(ya>0);
        ynegs(ya<0)=ya(ya<0);

        Txpos=-yapos+Rpos;
        Txnegs=-ynegs;
        TX=[Txpos Txnegs];
        TX_imprp=[Txpos(1:end/2) Txnegs(1:end/2) Txpos(end/2+1:end) Txnegs(end/2+1:end)];
        %Sign=[Sign TX];
        for lp=1:length(snr)
            
        TX_n=awgn(TX,snr(lp),'measured','db');
        Noise=TX_n-TX;
        Pn=sum(Noise.^2)/length(Noise);
        P_i=power_factor*Pn;
        A_i=sqrt(2*P_i);
        
        TX_n_imprp=awgn(TX_imprp,snr(lp),'measured','db');
        Noise_imprp=TX_n_imprp-TX_imprp;
        Pn_imprp=sum(Noise_imprp.^2)/length(Noise_imprp);
        P_i_imprp=power_factor*Pn_imprp;
        A_i_imprp=sqrt(2*P_i_imprp);
    
%         A_i=0;

        TX_prp=TX_n+A_i*[ones(1,Nffta) zeros(1,Nffta)];
        TX_imprp_distort=TX_n_imprp+A_i*[ones(1,Nffta) zeros(1,Nffta)];
        TX_imprp_distort_DC25=TX_n+A_i*[ones(1,Nffta/2) zeros(1,3*Nffta/2)];
        TX_imprp_distort_DC10=TX_n+A_i*[ones(1,Nffta/5-0.4) zeros(1,9*Nffta/5+0.4)];
        TX_distort=TX_n+A_i*(cos(2*pi*2/T_pwm*t));
%         TX_imprp_one=TX_n+tot_hmnc(lp)*cos(2*pi*2/T_pwm*t);
%         TX_imprp_one=TX_n+A/10*[ones(1,Nffta/2) zeros(1,3*Nffta/2)];

        RX=-(TX_n(1:end/2)-A)-TX_n(end/2+1:end);
        RX_fft=fft(RX);
        RX_demod=qamdemod(RX_fft,Ma);
        bits_RX=RX_demod(2:16);
        err(lp)=err(lp)+sum(biterr(bits_RX,bit_streama'));
        
        RX_prp=-(TX_prp(1:end/2)-A)-TX_prp(end/2+1:end);
        RX_fft_prp=fft(RX_prp);
        RX_demod_prp=qamdemod(RX_fft_prp,Ma);
        bits_RX_prp=RX_demod_prp(2:16);
        err_prp(lp)=err_prp(lp)+sum(biterr(bits_RX_prp,bit_streama'));
        
        RX_imprp_distort=[-(TX_imprp_distort(1:end/4)-A) - TX_imprp_distort(end/4+1:end/2),-(TX_imprp_distort(end/2+1:3*end/4)-A) - TX_imprp_distort(3*end/4+1:end)];
        RX_fft_imprp_distort=fft(RX_imprp_distort);
        RX_demod_imprp_distort=qamdemod(RX_fft_imprp_distort,Ma);
        bits_RX_imprp_distort=RX_demod_imprp_distort(2:16);
        err_imprp_distort(lp)=err_imprp_distort(lp)+sum(biterr(bits_RX_imprp_distort,bit_streama'));
        
        RX_imprp_25dc=-(TX_imprp_distort_DC25(1:end/2)-A)-TX_imprp_distort_DC25(end/2+1:end);
        RX_fft_imprp_25dc=fft(RX_imprp_25dc);
        RX_demod_imprp_25dc=qamdemod(RX_fft_imprp_25dc,Ma);
        bits_RX_imprp_25dc=RX_demod_imprp_25dc(2:16);
        err_imprp_25dc(lp)=err_imprp_25dc(lp)+sum(biterr(bits_RX_imprp_25dc,bit_streama'));
        
        RX_imprp_10dc=-(TX_imprp_distort_DC10(1:end/2)-A)-TX_imprp_distort_DC10(end/2+1:end);
        RX_fft_imprp_10dc=fft(RX_imprp_10dc);
        RX_demod_imprp_10dc=qamdemod(RX_fft_imprp_10dc,Ma);
        bits_RX_imprp_10dc=RX_demod_imprp_10dc(2:16);
        err_imprp_10dc(lp)=err_imprp_10dc(lp)+sum(biterr(bits_RX_imprp_10dc,bit_streama'));
        
        RX_distort=-(TX_distort(1:end/2)-A)-TX_distort(end/2+1:end);
        RX_fft_distort=fft(RX_distort);
        RX_demod_distort=qamdemod(RX_fft_distort,Ma);
        bits_RX_distort=RX_demod_distort(2:16);
        err_distort(lp)=err_distort(lp)+sum(biterr(bits_RX_distort,bit_streama'));
        
        
%         
%         RX_imprp_one=-(TX_imprp_one(1:end/2)-A)-TX_imprp_one(end/2+1:end);
%         RX_fft_imprp_one=fft(RX_imprp_one);
%         RX_demod_imprp_one=qamdemod(RX_fft_imprp_one,Ma);
%         bits_RX_imprp_one=RX_demod_imprp_one(2:16);
%         err_imprp_one(lp)=err_imprp_one(lp)+sum(biterr(bits_RX_imprp_one,bit_streama'));
        end
        l
end

BER=err/((Nffta/2-1)*L*ma);
BER_prp=err_prp/((Nffta/2-1)*L*ma);
BER_imprp_distort=err_imprp_distort/((Nffta/2-1)*L*ma);
BER_imprp_25dc=err_imprp_25dc/((Nffta/2-1)*L*ma);
BER_imprp_10dc=err_imprp_10dc/((Nffta/2-1)*L*ma);
BER_distort=err_distort/((Nffta/2-1)*L*ma);
% BER_imprp_one=err_imprp_one/((Nffta/2-1)*L*ma);

save('simber2.mat')

BER(BER==0)=10^-12;
BER_prp(BER_prp==0)=10^-12;
BER_imprp_distort(BER_imprp_distort==0)=10^-12;
BER_imprp_25dc(BER_imprp_25dc==0)=10^-12;
BER_imprp_10dc(BER_imprp_10dc==0)=10^-12;
BER_distort(BER_distort==0)=10^-12;

%%%%%%%%%%%%%%%%%% Modify x-axis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l=1:1    
        bit_streama=randi(Ma-1,Nffta/2-1,1); % Conventional OFDM
        %Bits=[Bits bit_streama];
        symb1a=qammod(bit_streama,Ma);
        %SYMB1A=[SYMB1A symb1a];
        
        tx_syma = [0 transpose(symb1a) 0 (fliplr(conj(transpose(symb1a))))];       
        
        ya = ifft(tx_syma);
        yapos=zeros(1,length(ya));
        ynegs=zeros(1,length(ya));

        yapos(ya>0)=ya(ya>0);
        ynegs(ya<0)=ya(ya<0);

        Txpos=-yapos+Rpos;
        Txnegs=-ynegs;
        TX=[Txpos Txnegs];
        TX_imprp=[Txpos(1:end/2) Txnegs(1:end/2) Txpos(end/2+1:end) Txnegs(end/2+1:end)];
        %Sign=[Sign TX];
        for lp=1:length(snr)
            
        TX_n=awgn(TX,snr(lp),'measured','db');
        RX=-(TX_n(1:end/2)-A)-TX_n(end/2+1:end);
        Ps=sum(RX.^2)/length(RX);
        Noise=TX_n-TX;
        Pn=sum(Noise.^2)/length(Noise);
        P_i=power_factor*Pn;
        
        SINR(lp)=10*log10(Ps/(Pn+P_i));
      
       end
        l
end

% Difference between axis with square wave and without square wave
diff=SINR-sinr; % nearly 5dB
sinr=sinr-5;
%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogy(sinr,BER,'->','linewidth',3)
hold all
semilogy(sinr,BER_prp,'--^','linewidth',2.5)
hold all
% semilogy(snr,BER_imprp_one,'-o','linewidth',2)
% hold all
semilogy(sinr,BER_imprp_distort,'->','linewidth',2)
hold all
semilogy(sinr,BER_imprp_25dc,'-s','linewidth',2)

hold all
semilogy(sinr,BER_imprp_10dc,'-*','linewidth',2)

hold all
semilogy(sinr,BER_distort,'-o','linewidth',2)

legend('U-OFDM without interference','Orthogonal U-OFDM with interference','Overlapped U-OFDM with different PWM frequency interference','Overlapped U-OFDM with 25% duty cycle PWM interference','Overlapped U-OFDM with 10% duty cycle PWM interference','Overlapped U-OFDM with single harmonic interference','location','SouthWest')
xlabel('SINR (dB)')
ylabel('BER')
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
set(findall(gcf,'marker','text'),'FontSize',12)

