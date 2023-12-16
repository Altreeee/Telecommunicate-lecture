function classpractice21()

clc
clear all

disp('----------------------------------------------------------------------')
disp('-------------------------class practice 2.1----------------------------')
disp('-------------------by Amber-GaoQi on 23/7/26---------------------')
disp('-----------------------------------------------------------------------')

% Input parameters from the user
D = input('distance:');   % 1:1:100
sigma = input('standard deviation of shadowing in dB:');
N = input('the number of samples for each distance:');
alpha = input('the pathloss exponent,alpha:'); % 2:0.5:4

figure; % Create a single figure outside the loops

% Create subplots arrangement based on the number of alpha values
num_alphas = length(alpha);
num_rows = ceil(sqrt(num_alphas));
num_cols = ceil(num_alphas / num_rows);

for a = 1:num_alphas
    for m = 1:length(D)
        for n = 1:N
             % Calculate The ratio of PR to P0 in dB scale with shadowing for each distance
             PrDB(m, n) = -10 * alpha(a) * log10(D(m)) + randn * sigma; 
                                                        %RANDN produces a standard normal distribution
        end
    end
    
    % Create a subplot for each alpha value
    subplot(num_rows, num_cols, a);
    
    % Plot each set of data on the subplot
    for m = 1:length(D)
        semilogx(D(m) * ones(1, N), PrDB(m, :), 'x')
        hold on
    end
    
    xlabel('distance in log scale')
    ylabel('the normalized received power in dB')
    str = sprintf('shadow effect with alpha = %g and sigma = %g dB', alpha(a), sigma);
    title(str);
end

% Adjust the layout to make the subplots fit nicely
sgtitle('Combined Plots'); % Add a common title for the combined figure

end
