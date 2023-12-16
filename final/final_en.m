function final()

clc
clear all

disp('----------------------------------------------------------------------')
disp('-------------------------final presentation----------------------------')
disp('-------------------by Amber-GaoQi on 23/9/2---------------------')
disp('-----------------------------------------------------------------------')

%%
    M=input('the number of cells: ');
    N=input('reuse factor: ');
    Nc=input('the number of frequency channels: ');
    Um=input('average number of users in each cell: ');
    snrdB=input('average SNR for each user in dB: ');
    snr=10^(snrdB/10);
    P=1;
    sigma=sqrt(P/snr);

%% Create an adjacency matrix to represent neighboring relationships
min_distance_threshold = 0.9999999; % Minimum distance between cell centers

% Create an adjacency matrix to represent neighboring relationships, initially set to 0
adjacency_matrix = zeros(M, M);

% Initialize lists for center point coordinates
x_center_positions = zeros(1, M);
y_center_positions = zeros(1, M);

% Generate the first center point (0, 0)
x_positions(1) = 0;
y_positions(1) = 0;

% Define polar coordinate parameters
radii = 1; % Radius
angles =[0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3]; % Angles at 60-degree intervals

% Generate the remaining center points
for i = 2:M
    n=2;
    found_position = false; % Add a flag variable to control the outer loop
    while found_position == false
        % Generate x and y points, ensuring they are integers, around the origin in concentric circles
        h=1+floor((i-2)/6);
        if mod(n, 6) == 1
            x_positions(n) = x_center_positions(h)+(radii * sin(angles(1)));
            y_positions(n) = y_center_positions(h) + (radii * cos(angles(1)));
        elseif mod(n, 6) == 2
            x_positions(n) = x_center_positions(h)+(radii * sin(angles(2)));
            y_positions(n) = y_center_positions(h) + (radii * cos(angles(2)));
        elseif mod(n, 6) == 3
            x_positions(n) = x_center_positions(h)+(radii * sin(angles(3)));
            y_positions(n) = y_center_positions(h)+ (radii * cos(angles(3)));
        elseif mod(n, 6) == 4
            x_positions(n) = x_center_positions(h)+(radii * sin(angles(4)));
            y_positions(n) = y_center_positions(h) + (radii * cos(angles(4)));
        elseif mod(n, 6) == 5
            x_positions(n) = x_center_positions(h)+(radii * sin(angles(5)));
            y_positions(n) = y_center_positions(h)+ (radii * cos(angles(5)));
        elseif mod(n, 6) == 0
            x_positions(n) = x_center_positions(h)+(radii * sin(angles(6)));
            y_positions(n) = y_center_positions(h)+(radii * cos(angles(6)));
        end

        % Check the distance between this center point and existing center points
        too_close = false;
        for j = 1:i-1
            distance = sqrt((x_positions(n) -  x_center_positions(j))^2 + (y_positions(n) -  y_center_positions(j))^2);
            if distance < min_distance_threshold
                too_close = true;
                n=n+1;
                break;
            end
        end

        if ~too_close
            % If the distance is suitable, add the coordinates of this center point to the list
            x_center_positions(i) = x_positions(n);
            y_center_positions(i) = y_positions(n);
            n=n+1;
            found_position = true; % Set the flag variable to true
            break;
        end

        if n > 1000 % Add a maximum number of attempts to prevent an infinite loop
            error('Unable to generate suitable positions.');
        end
    end
    
    if ~found_position
        error('Unable to generate suitable positions.');
    end
end

% Output the positions of each cell
% disp('Cell positions:');
% disp([x_center_positions; y_center_positions]);

% Calculate neighboring relationships based on cell positions
for i = 1:M
    for j = 1:M
        if i ~= j
            % Calculate the distance between cell i and cell j
            distance = sqrt((x_center_positions(i) - x_center_positions(j))^2 + (y_center_positions(i) - y_center_positions(j))^2);

            % If the distance is less than the threshold, set the corresponding position in the adjacency matrix to 1
            if distance < 1.1
                adjacency_matrix(i, j) = 1;
            end
        end
    end
end

% Output the adjacency matrix
disp('Adjacency matrix:');
disp(adjacency_matrix);

%%
% Initialize a frequency list
fr = 1:N;
% Create an array to store the frequency used by each cell
grid_fr = zeros(1, M);
% Start the allocation process
for i = 1:M
    % Find the frequencies used by neighboring cells of the current cell
    neighbor_fr = grid_fr(adjacency_matrix(i, :) == 1);
    % Find available frequencies that have not been used
    available_fr = setdiff(fr, neighbor_fr);
    % Randomly select an available frequency and assign it to the current cell
    if ~isempty(available_fr)
        selected_fr = available_fr(randi(length(available_fr)));
        grid_fr(i) = selected_fr;
    else
        % If no available frequencies are found, select any frequency randomly
        selected_fr = fr(randi(length(fr)));
        grid_fr(i) = selected_fr;
    end
end
% Print the frequencies used by each cell
for i = 1:M
    disp(['Cell ' num2str(i) ' uses frequency f' num2str(grid_fr(i))]);
end

%%
% First, set all diagonal elements of the adjacency_matrix to 1.
% Then, check the positions j in columns of adjacency_matrix(i, :) == 1 representing cell j
% If they have the same frequency as cell i, change the 1 at position (i,j) in adjacency_matrix to r
r = 0.7;
% Set diagonal elements from 0 to 1
adjacency_matrix = adjacency_matrix + eye(M);

% Check positions in each row where adjacency_matrix(i, :) == 1
for i = 1:M
    neighbors = find(adjacency_matrix(i, :) == 1); % Find indices of neighboring cells
    current_frequency = grid_fr(i); % Get the current cell's frequency
    
    for j = neighbors
        if i ~= j % Ensure that i and j are not the same cell
            if grid_fr(j) == current_frequency % If the frequencies are the same
                adjacency_matrix(i, j) = r; % Change the 1 at (i, j) to r
            end
        end
    end
end

rho = adjacency_matrix;
% Traverse the adjacency matrix and change 1 to 0 except for the diagonal elements
for i = 1:M
    for j = 1:M
         if i ~= j && rho(i, j) == 1
            rho(i, j) = 0;
         end
     end
 end
 
%%    
    Nc_ = zeros(1, N);  % Initialize the Nc_ array
    for i = 1:N
        if i ~= N
            Nc_(i) = round(Nc / N);
        else
            % Calculate Nc_(N): Nc minus the sum of all previous Nc_(i)
            Nc_(N) = Nc - sum(Nc_(1:N-1));
        end

        % Generate a part of f(N, 1:Nc_(N))
        if i == 1
            g(i, 1:Nc_(i)) = 1:Nc_(i);
        else
            g(i, 1:Nc_(i)) = sum(Nc_(1:i-1)) + 1:sum(Nc_(1:i-1)) + Nc_(i);
        end
    end

    for i=1:M
        num=grid_fr(i);
        f(i, 1:Nc_(num))=g(num,1:Nc_(num));
    end

    % Setting the number of active users for all the cells
    for c=1:M
        U(c)=poissrnd(Um);
        Nchannel(c)=length( find( f(c , : )>0));%check
        if U(c)>Nchannel(c)
            fprintf('We cannot serve all the %d users in cell %d\n',Nchannel(c),c)
        end
    end
    
    % Generate fading channels for each serving user from all base stations
    for c1=1:M
        for c2=1:M
            for u=1:min(Nchannel(c2),U(c2))
                h(c1,c2,u)=(randn+1i*randn)/sqrt(2); % From cell c1 to user u in cell c2
            end
        end
    end
    
    % Compute the user capacity for each cell
    for c=1:M
        for u=1:min(Nchannel(c),U(c))
            I(c,u)=0;
            for c_=1:M
                if rho(c_,c)~=0&&c_~=c
                    if U(c_)>=u&&Nchannel(c_)>=u
                         I(c,u)=I(c,u)+rho(c_,c)*P*abs(h(c_,c,u))^2;
                    end
                end
            end
            sinr(c,u)=P*abs(h(c,c,u))^2/(I(c,u)+sigma^2);
            capacity(c,u)=log2(1+sinr(c,u));
        end
        cell_capacity(c)=sum(capacity(c, : ));
    end
    total_net=sum(cell_capacity);
    fprintf('Capacity is %g\n',total_net);
   
end
