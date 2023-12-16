function classpractice41()
    clc
    clear all
    
    disp('----------------------------------------------------------------------')
    disp('-------------------------class practice 4.1----------------------------')
    disp('-------------------by Amber-GaoQi on 23/8/12---------------------')
    disp('-----------------------------------------------------------------------')
    
    Nf=input('the number of frequency channels: ');
    Nsim=input('the number of random channels for simulation: ');
    avSNRdB=input('average SNR for each frequency channel in dB: '); %0:10:
    avSNR=10.^(avSNRdB/10);
    sigma=sqrt(1./avSNR);
    P=1;
    
    for m=1:length(avSNR)
        for n=1:Nsim
            H=(randn(1,Nf)+1i*randn(1,Nf))/sqrt(2);
            water_level=fzero(@(x) fun_waterlevel(x,sigma(m),H,P),P/2);
            C(m,n)=0;
            C_equal(m,n)=0;
            for f=1:Nf
                p(f)=max([0 water_level-sigma(m)^2/abs(H(f))^2]);
                C(m,n)=C(m,n)+log2(1+p(f)*abs(H(f))^2/sigma(m)^2);
                C_equal(m,n)=C_equal(m,n)+log2(1+P/Nf*abs(H(f))^2/sigma(m)^2);
            end
        end
        Cav(m)=mean(C(m,:));
        Cav_e(m)=mean(C_equal(m,:));
        fprintf('SNR=%gdB is completed\n',avSNRdB(m));
    end
    
    plot(avSNRdB,Cav,'b-')
    hold on
    grid on
    plot(avSNRdB,Cav_e,'r--')
    xlabel('Average SNR per channel in dB')
    ylabel('Average capacity')
    legend('water-filling','equal power allocation')
    
end