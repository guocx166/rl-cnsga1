function action = epsilonGreedy(Q_table, state, epsilon)
% epsilon-greedy算法来选择动作
    if rand < epsilon
        % 探索：随机选择一个动作
        action = randi(size(Q_table, 2)); % 选择范围在 1 到 num_actions 之间
    else
        % 利用：选择当前Q表中最大Q值的动作
        [~, action] = max(Q_table(state, :)); % 返回最大值的索引
    end
end