function [transport1,transport2,transport3] = generateTransport(n)
% 构造交通方式的矩阵
% 三个矩阵分别代表货车、火车和直升机
% 使用函数
    min_n1 = 1000; 
    % 生成一个 n x n 的随机矩阵  
    upper_triangle = triu(randi([0, 1], n, n), 1);
    transport1 = upper_triangle + upper_triangle' + eye(n);
    % transport1 = randi([0, 1], n, n); % 生成全 0 和 1 矩阵  
    % 将对角线元素设为 0 
    transport1(1:n+1:end) = 0; 
    a = sum(transport1(:))/2;
    while a < min_n1
        row = randi(n);  
        col = randi(n);  
        if row ~= col && transport1(row, col) == 0  
        transport1(row, col) = 1;  
        transport1(col, row) = 1;  % 保持对称  
        min_n1 = min_n1 - 1;  % 减少需要填充的1的数量  
        end 
    end
%     transport2 = zeros(n,n);
%     transport3 = zeros(n,n);
    transport2 = randi([0, 1], n, n); % 生成全 0 和 1 矩阵     
    transport2(1:n+1:end) = 0; % 将对角线元素设为 0  
    transport3 = randi([0, 1], n, n); % 生成全 0 和 1 矩阵  
    transport3(1:n+1:end) = 0; % 将对角线元素设为 0 
end
