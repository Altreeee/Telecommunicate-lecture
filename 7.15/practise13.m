function classpractice13(void)
%%
clc
clear all


%% the introduction of this program
disp('------------------------------------');
disp('This is Class Practice 1.3');
disp('This is done by student GaoQi-Amber, 16/7/2023');
disp('given the K and L and we can get N for all combinations');
disp('------------------------------------')

%% get the value
K=input('what is the range of k?');
L=input('what is the range of l?');
alpha=input('what is the the pathloss exponent?');

fprintf('\n')
disp('the outputs are: ')
fprintf('\n')

%% Calculate the size of N and output
%prompt
fprintf('the outputs are：\n');

%now we add a new layer to loop all the alpha's situation
%each loop a would add 0.5
 for a=1:0.5:alpha
        for k = 1:K  
                for l = 1:L
                    N = k^2 + k*l + l^2;    %Calculate N
                    
                    SIRapp=(sqrt(3*N))^alpha/6; %Calculate the approxinated SIR
                    SIPapp_dB=10*log10(SIRapp); %Convert it to dB form
                    
                    Q=sqrt(3*N);
                    %Calculate the accurate SIR
                    SIRacc=1/(((2*(Q+1)^alpha+(Q-1)^alpha)/(Q^2-1)^alpha)+(((Q+0.5)^alpha+(Q-0.5)^alpha)/(Q^2-0.25)^alpha)+(1/Q^alpha));
                    %Convert it to dB form
                    SIRacc_dB=10*log10(SIRacc);
                    
                    fprintf('(K, L) = (%d, %d), alpha is %g, N = %d and the approxinated SIR is %g dB the accurate SIR is %g dB\n', k, l, a, N,SIPapp_dB,SIRacc_dB);    %output
                end
        end
 end
