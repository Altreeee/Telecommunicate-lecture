function classpractice42()
    clc
    clear all
    
    disp('----------------------------------------------------------------------')
    disp('-------------------------class practice 4.2----------------------------')
    disp('-------------------by Amber-GaoQi on 23/8/19---------------------')
    disp('-----------------------------------------------------------------------')
    
    U=input('the number of cells: ');
    Nf=input('the number of frequency channels: ');
    SNRdB=input('SNR in dB: '); %0:10:
    SNR=10.^(SNRdB/10);
    P=1;
    sigma=sqrt(1./SNR);
    Nsim=input('the number of random simulations');
    
    for m=1:length(SNR)
        for n=1:Nsim
            H=(randn(Nf,U,U)+1i*randn(Nf,U,U))/sqrt(2);
            p=sWFpa(H,sigma(m),P);
            for u=1:U%each user
                C(u)=0;
                for f=1:Nf%each frequency
                    c(f)=sigma(m)^2;
                    for u_=1:U
                        if u_~=u
                            c(f)=c(f)+p(f,u_)*abs(H(f,u_,u))^2;
                        end
                    end
                    sinr(f,u)=p(f,u)*abs(H(f,u,u))^2/c(f);
                    C(u)=C(u)+log2(1+sinr(f,u));
                end
            end
            Net_Cap(m,n)=sum(C);
        end
       avC(m)=mean(Net_Cap(m,:));
       fprintf('SNR=%gdB is %g\n',SNRdB(m),avC(m));
    end
    
    Yplot=input('Plot the results (Y=1/N=0)? ');
    if Yplot==1
        plot(SNRdB,avC);
        grid on
        xlabel('Average SNR in dB')
        ylabel('Average capacity')
    end
    hold on
    
    plot(SNRdB,Cav_e,'r--')
    
    legend('water-filling','equal power allocation')
    
end