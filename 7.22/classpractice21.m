function classpractice21()

clc
clear all

disp('----------------------------------------------------------------------')
disp('-------------------------class practice 2.1----------------------------')
disp('-------------------by Amber-GaoQi on 23/7/26---------------------')
disp('-----------------------------------------------------------------------')

D=input('distance:');   %1:1:100
sigma=input('standard deviation of shadowing in dB:');
N=input('the number of samples for each distance:');
alpha=input('the pathloss exponent,alpha:'); %2:0.5:4

for a=1:length(alpha)
    for m=1:length(D)
        for n=1:N
             PrDB(m,n)=-10*alpha(a)*log10(D(m))+randn*sigma;
        end
    end
    
    figure(a)

    for m=1:length(D)
        semilogx(D(m)*ones(1,N),PrDB(m,:),'x')
        hold on
    end

    xlabel('distance in log scale')
    ylabel('the normalised received power in dB')
    str=sprintf('shadow effect with alpha = %g and sigma =%g dB',alpha,sigma);
    title(str);
end    
end
