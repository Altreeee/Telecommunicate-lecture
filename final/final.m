function final()

clc
clear all

disp('----------------------------------------------------------------------')
disp('-------------------------final presentation----------------------------')
disp('-------------------by Amber-GaoQi on 23/9/2---------------------')
disp('-----------------------------------------------------------------------')

%%
    M=input('the number of cells');
    N=input('reuse factor: ');
    Nc=input('the number of frequency channels: ');
    Um=input('average number of users in each cell: ');
    snrdB=input('average SNR for each user in dB: ');
    snr=10^(snrdB/10);
    P=1;
    sigma=sqrt(P/snr);



%% 创建表示相邻关系的邻接矩阵
min_distance_threshold = 0.9999999; % 每两个cell中心之间的最小距离

% 创建表示相邻关系的邻接矩阵，初始都设为0
adjacency_matrix = zeros(M, M);

% 初始化中心点坐标列表
x_center_positions = zeros(1, M);
y_center_positions = zeros(1, M);

% 生成第一个中心点（0，0）
x_positions(1) = 0;
y_positions(1) = 0;

% 定义极坐标生成参数
radii = 1; % 半径
angles =[0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3]; % 60度间隔的角度

% 生成剩余的中心点
for i = 2:M
    n=2;
    found_position = false; % 添加标志变量来控制外部循环
    while found_position == false
                    % 生成x和y点，确保都是整数，围绕原点一圈一圈生成
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

                    

                    % 检查该中心点与已有中心点的距离
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
                        % 如果距离合适，将该中心点的坐标添加到列表中
                        x_center_positions(i) = x_positions(n);
                        y_center_positions(i) = y_positions(n);
                        n=n+1;
                        found_position = true; % 设置标志变量为true
                        break;
                    end
                    
                    if n > 1000 % 添加一个最大尝试次数的限制，防止无限循环
                         error('无法生成合适的位置。');
                    end
    end
                        if ~found_position
                            error('无法生成合适的位置。');
                        end
end

% 输出每个cell的位置
%disp('cell的位置：');
%disp([x_center_positions; y_center_positions]);

% 根据cell位置计算相邻关系
for i = 1:M
    for j = 1:M
        if i ~= j
            % 计算cell i和cell j之间的距离
            distance = sqrt((x_center_positions(i) - x_center_positions(j))^2 + (y_center_positions(i) - y_center_positions(j))^2);

            % 如果距离小于阈值，将邻接矩阵中对应位置设为1
            if distance < 1.1
                adjacency_matrix(i, j) = 1;
            end
        end
    end
end

% 输出邻接矩阵
disp('adjacency_matrix：');
disp(adjacency_matrix);

%%
% 初始化频率列表
fr = 1:N;

% 创建一个数组来存储每个cell的频率
grid_fr = zeros(1, M);

% 开始分配过程
for i = 1:M
    % 找到与当前cell相邻的cell的频率
    neighbor_fr = grid_fr(adjacency_matrix(i, :) == 1);
    
    % 找到尚未使用的频率
    available_fr = setdiff(fr, neighbor_fr);
    
    % 随机选择一个可用频率并将其分配给当前cell
    if ~isempty(available_fr)
        selected_fr = available_fr(randi(length(available_fr)));
        grid_fr(i) = selected_fr;
    else
        % 如果没有可用频率，则任选一个频率分配
        selected_fr = fr(randi(length(fr)));
        grid_fr(i) = selected_fr;
    end
end

% 打印每个cell的频率
for i = 1:M
    disp(['cell ' num2str(i) ' used f' num2str(grid_fr(i))]);
end

%%
%首先将adjacency_matrix矩阵中的对角线数字全部从0变为1.
%然后检查每行adjacency_matrix(i, :) == 1的位置j处列所代表的cellj所对应的频率
%如果与celli所对应的频率相同，则将adjacency_matrix矩阵中(i,j)该处的1改为r
r=0.7;
% 将对角线上的 0 改为 1
adjacency_matrix = adjacency_matrix + eye(M);

% 检查每行 adjacency_matrix(i, :) == 1 的位置
for i = 1:M
    neighbors = find(adjacency_matrix(i, :) == 1); % 找到相邻cell的索引
    current_color = grid_fr(i); % 获取当前cell的频率
    
    for j = neighbors
        if i ~= j % 确保 i 和 j 不是同一个cell
            if grid_fr(j) == current_color % 如果频率相同
                adjacency_matrix(i, j) = r; % 将 (i, j) 处的 1 改为 r
            end
        end
    
    end
end

rho = adjacency_matrix;
% 遍历邻接矩阵，将除了对角线以外的 1 改为 0
for i = 1:M
    for j = 1:M
         if i ~= j && rho(i, j) == 1
            rho(i, j) = 0;
         end
     end
 end
 
    
    %%
    Nc_ = zeros(1, N);  % 初始化 Nc_ 数组
    for i = 1:N
        if i ~= N
            Nc_(i) = round(Nc / N);
        else
            % 计算 Nc_(N)：Nc 减去前面所有 Nc_(i) 的和
            Nc_(N) = Nc - sum(Nc_(1:N-1));
        end

        % 生成 f(N, 1:Nc_(N)) 的一部分
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

    %setting the number of active users for all the cells
    for c=1:M
        U(c)=poissrnd(Um);
        Nchannel(c)=length( find( f(c , : )>0));%check
        if U(c)>Nchannel(c)
            fprintf('we cannot server all the %d users in cell %d\n',Nchannel(c),c)
        end
    end
    
    %generate fading channels for each serving user from all base stations
    for c1=1:M
        for c2=1:M
            for u=1:min(Nchannel(c2),U(c2))
                h(c1,c2,u)=(randn+1i*randn)/sqrt(2); %from cell c1 to user u in cell c2
            end
        end
    end
    
    %compute the user capacity for each cell
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
    fprintf('capacity is %g\n',total_net);
   
end


    