function diversity = calculateDiversity(objectives)
    % chromosome: 输入种群，元胞形式
    % V: 决策变量的数量
    % 默认假设目标数量为2，因为函数的描述是基于两个目标
    % 提取帕累托前沿解
    % 使用
%     isPareto = cell2mat(chromosome(:, V + M + 1)); % 假定M+V+1列为1表示帕累托前沿
%     minNumber = min(isPareto);
%     paretoSolutions = chromosome(isPareto == minNumber, :);
%     % paretoSolutions = chromosome(isPareto == 1, :);
%     
%     % 提取目标值
%     objectives = cell2mat(paretoSolutions(:, V+1:V+2));
    A = size(objectives, 1); % 帕累托前沿解的个数
    
    % 计算解之间的指定距离
    d = zeros(A, 1);
    for i = 1:A
        % 计算每个解与其他解的目标距离，求最小值
        distances = inf;
        for j = 1:A
            if i ~= j
                distance = sum(abs(objectives(i, :) - objectives(j, :)));
                if distance < distances
                    distances = distance;
                end
            end
        end
        d(i) = distances;
    end
    
    % 平均距离
    avg_d = sum(d) / (A-1);
    
    % 方差
    delta = sum((d - avg_d).^2) / (A-1);
    
    % 计算离散度
    diversity = delta / avg_d;
end

