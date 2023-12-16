function classpractice12(void)
%%
clc
clear all


%% the introduction of this program
disp('------------------------------------');
disp('This is Class Practice 1.2');
disp('This is done by student GaoQi-Amber, 16/7/2023');
disp('We want to know what is the SIR when N is 7 or 12');
disp('------------------------------------')

%% get the value
%get the value of alpha from user
alpha = input('please input the alpha：');


%% Calculate the SIR and output
disp('when the N is 7: ');
N=7;
Q=sqrt(3*N);
SIR=1/(((2*(Q+1)^alpha+(Q-1)^alpha)/(Q^2-1)^alpha)+(((Q+0.5)^alpha+(Q-0.5)^alpha)/(Q^2-0.25)^alpha)+(1/Q^alpha));
fprintf('SIR is %g\n',SIR);

disp('when the N is 12: ');
N=12;
Q=sqrt(3*N);
SIR=1/(((2*(Q+1)^alpha+(Q-1)^alpha)/(Q^2-1)^alpha)+(((Q+0.5)^alpha+(Q-0.5)^alpha)/(Q^2-0.25)^alpha)+(1/Q^alpha));
fprintf('SIR is %g\n',SIR);
