function classpractice22()
clc
clear all

disp('----------------------------------------------------------------------')
disp('-------------------------class practice 2.2----------------------------')
disp('-------------------by Amber-GaoQi on 23/7/27---------------------')
disp('-----------------------------------------------------------------------')


disp('0-BPSK')              %1 bit per symbol
disp('1-QPSK')              % 2 bits per symbol
disp('2-16-QAM')        % 4 bits per symbol
disp('3-64-QAM')        % 6 bits per symbol

fprintf('\n')

Mod=input('choose 0, 1, 2 or 3: ');
N=input('the number of random samples: '); %200,1000
snrdB=input('average SNR in dB: '); %10,15,20,30
snr=10^(snrdB/10);
p=1;
sigma=sqrt(p/snr); %snr=p/sigma^2

switch Mod
    case 0%BPSK
        c=[-1 1];
        s=c(randi(2,1,N));
        noise = (randn(1,N)+1i*randn(1,N))/sqrt(2)*sigma;
        y=s+noise;
        plot(y,'x')
        grid on
    case 1%QPSK
        c=[-1-1i -1+1i 1-1i 1+1i];
        pav=mean(abs(c).^2);
        c=c/sqrt(pav);
        s=c(randi(4,1,N));
        noise=(randn(1,N)+1i*randn(1,N))/sqrt(2)*sigma;
        y=s+noise;
        plot(y,'x')
        grid on
    case 2%16-QAM
        sI=[-3 -1 1 3 ];
        sQ=[-3 -1 1 3];
        k=1;
        for m=1:4
            for n=1:4
                c(k)=sI(m)+1i*sQ(n);
                k=k+1;
            end
        end
        pav=mean(abs(c).^2);
        c=c/sqrt(pav);
        s=c(randi(16,1,N));
        noise=(randn(1,N)+1i*randn(1,N))/sqrt(2)*sigma;
        y=s+noise;
        plot(y,'x')
        grid on
    case 3%64-QAM
        sI=[-7 -5 -3 -1 1 3 5 7];
        sQ=[-7 -5 -3 -1 1 3 5 7];
        k=1;
        for m=1:8
            for n=1:8
                c(k)=sI(m)+1i*sQ(n);
                k=k+1;
            end
        end
        pav=mean(abs(c).^2);
        c=c/sqrt(pav);
        noise=(randn(1,N)+1i*randn(1,N))/sqrt(2)*sigma;
        s=c(randi(64,1,N))+noise;
        plot(s,'x')
        grid on
end
end

        
        