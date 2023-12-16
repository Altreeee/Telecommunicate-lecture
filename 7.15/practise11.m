function classpractice11(void)
%%
clc
clear all


%% the introduction of this program
disp('------------------------------------');
disp('This is Class Practice 1.1');
disp('This is done by student GaoQi-Amber, 16/7/2023');
disp('given the K and L and we can get N for all combinations');
disp('------------------------------------')

%% get the value
% The user is prompted for K and L
disp('please input the range of K：');
%get the value of K from user
K = input('K = ');
disp('please input the range of L：');
L = input('L = ');
disp('please input the range of L：');

%% Calculate the size of N and output
%prompt
fprintf('the outputs are：\n');
%In the first layer cycle, L grows from 1 to L, and K increases by 1
for k = 1:K
    %The second layer loops
    for l = 1:L
        N = k^2 + k*l + l^2;    %Calculate N
        fprintf('(K, L) = (%d, %d), N = %d \n', k, l, N);    %output
    end
end
