function classpractice26()
clc
clear all

disp('----------------------------------------------------------------------')
disp('-------------------------class practice 2.6----------------------------')
disp('-------------------by Amber-GaoQi on 23/8.12---------------------')
disp('-----------------------------------------------------------------------')


disp('0-BPSK')              %1 bit per symbol
disp('1-QPSK')              % 2 bits per symbol
disp('2-16-QAM')        % 4 bits per symbol
disp('3-64-QAM')        % 6 bits per symbol

fprintf('\n')

Mod=input('choose 0, 1, 2 or 3: ');
snrdB=input('average SNR in dB: '); %10,15,20,30
Nsim=input('the number of random transmissions: ')
Nch=input('the number of channel: '); 


snr_b=10.^(snrdB/10);


switch Mod
    case 0%BPSK
        NofP=2;
    case 1%QPSK
        NofP=4;
    case 2%16-QAM
        NofP=16;
    case 3%64-QAM
        NofP=64;
end
snr=snr_b*log2(NofP);
L=length(snr_b);

A=1;
for m=1:L
    sigma(m)=sqrt(1/2/snr(m));
end

for m=1:L
    for k=1:Nch
        [s,c]=mod_symbols(Nsim,Mod);
        noise=(randn(1,Nsim)+1i*randn(1,Nsim))/sqrt(2)*sigma(m);
        h=(randn+1i*randn)/sqrt(2);
        y=h*s+noise;

        if Mod==0
            y=real(y);  %BPSK
        end

        Nerror(m)=0;

        %get the simbol error rate
        for n=1:Nsim
            [unused,ind]=min(abs(y(n)/h-c).^2);   %c is a symbol
            %[minimum number, the index of it]
            s_detect(n)=c(ind);
            if s(n)~=s_detect(n)
                Nerror(m)=Nerror(m)+1;
            end
        end
    end 
    ser(m)=Nerror(m)/(Nsim*Nch);
    ber(m)=ser(m)/log2(NofP); %another way of calculate 
    
    fprintf('SNR=%g dB done \n',snrdB(m));
end

semilogy(snrdB,ber,'o-');
grid on
xlabel('SNR per bit in dB')
ylabel('BER')
switch Mod
    case 0
        title('BPSK with fading')
    case 1
        title('QPSK with fading')
    case 2
        title('16-QAM with fading')
    case 3
        title('64-QAM with fading')
end
end


        
        