function classpractice23()
clc
clear

disp('----------------------------------------------------------------------')
disp('-------------------------class practice 2.3----------------------------')
disp('-------------------by Amber-GaoQi on 23/7/26---------------------')
disp('-----------------------------------------------------------------------')

snrdB=input('the range of SNR in dB for simulations: ');%0:1:10
snr=10.^(snrdB/10);
A=1;
sigma = sqrt(A^2/2./snr);

Nsim=input('the number of random samples: ')

tic

c=[-1 1];   %A=1, BPSK

for m=1:length(snr)
    Nerror(m)=0;    %initial the error
    
    for n=1:Nsim
        s(n)=c(randi(2));
        noise(n)=randn*sigma(m);
        y(n)=s(n)+noise(n);
        
        s_det(n)=(y(n)>0)*1+(y(n)<=0)*(-1); %get a new s after receive s with noise
        %use matrices to avoid for loop
        %if y(n)>0
        %   s_det(n)=1;
        %else
        %   s_det(n)=-1;
        %end
        
        if s(n)~=s_det(n) %there's a mistake
            Nerror(m)=Nerror(m)+1;
        end
    end
    
    BER(m)=Nerror(m)/Nsim; %average rate
    fprintf('SMR=%g BER=%g done\n', snrdB(m),BER(m))
end

toc

figure(1)
semilogy(snrdB,BER); %draw the plot
xlabel('SNR in dB')
ylabel('simulated BER')
title('BPSK')
grid on;

BER_anal=analytical_ber(snrdB);

%campare two plots
figure(2);
semilogy(snrdB,BER,'o-');
hold on
semilogy(snrdB,BER_anal,'s-');
xlabel('SNR in dB')
ylabel('BER')
legend('simulated BER','analytical BER')
title('BPSK')
grid on;

end