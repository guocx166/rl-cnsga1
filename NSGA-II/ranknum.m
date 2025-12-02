function [sortedDict,sorted_indices] = ranknum(filename,DemandPoint)
data1 = readtable(filename, 'ReadVariableNames', false);
data1 = table2array(data1); 
data = zeros(size(data1, 1), 3, 3);
for i = 1:size(data1, 1) % 遍历每个受灾点  
    data(i, :, :) = reshape(data1(i, :), [3, 3]); % 将每行数据重塑为 3x3 矩阵  
end
weights = [0.3, 0.4, 0.3]; 

% 计算受灾点数量
n_count = size(data, 1); % 受灾点数量
m_count = size(data, 2); % 指标数量 (这里为3个指标)

% 1. 计算每个三角模糊数的期望值
expected_values = zeros(n_count, m_count);
for i = 1:n_count
    for j = 1:m_count
        % 计算每个模糊三角数的期望值
        l = data(i,j,1);
        m = data(i,j,2);
        u = data(i,j,3);
        expected_values(i,j) = (l + 4 * m + u) / 6;  % 三角模糊数的期望值
    end
end

% 2. 归一化处理
norm_data = zeros(n_count, m_count);
for j = 1:m_count
    max_value = max(expected_values(:, j));  % 找到最大值
    min_value = min(expected_values(:, j));  % 找到最小值
    norm_data(:, j) = (expected_values(:, j) - min_value) / (max_value - min_value);  % 归一化
end

% 3. 加权归一化矩阵
weighted_data = norm_data .* weights; 

% 4. 确定理想解和反理想解
ideal_solution = max(weighted_data);   % 各指标的理想解
anti_ideal_solution = min(weighted_data); % 各指标的反理想解

% 5. 计算与理想解和反理想解的距离
D_plus = sqrt(sum((weighted_data - ideal_solution).^2, 2)); % 与理想解的距离
D_minus = sqrt(sum((weighted_data - anti_ideal_solution).^2, 2)); % 与反理想解的距离

% 6. 计算相对接近度
relative_closeness = D_minus ./ (D_plus + D_minus);

% 7. 排序
[~, sorted_indices] = sort(relative_closeness, 'descend');

% 8. 输出排序结果

% 创建字典（Map）  
sortedDict = containers.Map('KeyType', 'double', 'ValueType', 'double');  
% 将节点和排序结果存入字典  
for i = 1:length(DemandPoint)  
    sortedDict(DemandPoint(i)) = sorted_indices(i); % 将每个点和值对应加入字典  
end 
% disp('受灾点的排序结果:');
% for i = 1:length(sorted_indices)
%     fprintf('受灾点 %d: 相对接近度 = %.4f\n', sorted_indices(i), relative_closeness(sorted_indices(i)));
% end
end

