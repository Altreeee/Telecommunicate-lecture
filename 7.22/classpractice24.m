function classpractice24()
clc
clear

disp('----------------------------------------------------------------------')
disp('-------------------------class practice 2.4----------------------------')
disp('-------------------by Amber-GaoQi on 23/7/26---------------------')
disp('-----------------------------------------------------------------------')

mTYPE=input('choose the mod 0,1,2 or 3: ');
snrdB_bit=input('the range of SNR in dB for simulations: ');
snr_bit=10.^(snrdB_bit/10);
Nsim=input('the number of random samples: ');


switch mTYPE
    case 0
        NofP=2; %the number of constellation points
    case 1
        NofP=4;
    case 2
        NofP=16;
    case 3
        NofP=64;
end
snr=snr_bit*log2(NofP);
L=length(snr_bit);
P=1;
sigma=sqrt(P./snr);

for m=1:L
    [s,c]=mod_symbols(Nsim,mTYPE);
    noise=(randn(1,Nsim)+1i*randn(1,Nsim))/sqrt(2)*sigma(m);
    y=s+noise;
    
    if mTYPE==0
        y=real(y);  %BPSK
    end
    
    Nerror(m)=0;
    
    %get the simbol error rate
    for n=1:Nsim
        [unused,ind]=min(abs(y(n)-c).^2);   %c is a symbol
        %[minimum number, the index of it]
        s_detect(n)=c(ind);
        if s(n)~=s_detect(n)
            Nerror(m)=Nerror(m)+1;
        end
    end
    ser(m)=Nerror(m)/Nsim
    ber(m)=1-(1-ser(m))^(1/log2(NofP));
    ber1(m)=ser(m)/log2(NofP); %another way of calculate 
    
    fprintf('SNR=%gdB and BER=%g done \n',snr_bit(m),ber(m));
end

semilogy(snrdB_bit,ber,'-o');
grid on
xlabel('SNR per bit in dB')
ylabel('BER')
switch mTYPE
    case 0
        title('BPSK')
    case 1
        title('QPSK')
    case 2
        title('16-QAM')
    case 3
        title('64-QAM')
end
end
