function [path] = dijkstra1(w, start, terminal, a)  
    % Dijkstra算法实现  
    % w: 带权邻接矩阵  
    % start: 起点编号 (1到n)  
    % terminal: 终点编号 (1到n)  
    % a: 不可达节点编号  

    [n1, ~] = size(w);  
    
    % 将不可达节点的行和列设置为无穷大  
    if terminal ~= a  
        w(a, :) = inf;  
        w(:, a) = inf;  
    end  
    
    n = n1;  % 顶点数  
    label = inf(n, 1);  % 到各点的最短路径  
    label(start) = 0;   % 起始点到自身的距离是0  
    f = nan(n, 1);      % 前驱节点  
    S = false(n, 1);    % 已访问的节点集  

    % 主循环  
    for count = 1:n  
        % 找到当前未访问且最小的距离节点  
        u = -1;  
        minDistance = inf;  
        for i = 1:n  
            if ~S(i) && label(i) < minDistance  
                minDistance = label(i);  
                u = i;  
            end  
        end  
        
        % 如果找不到可达节点，退出循环  
        if u == -1  
            break;  
        end  
        
        % 标记为已访问  
        S(u) = true;  

        % 更新邻接节点的距离  
        for v = 1:n  
            if ~S(v) && w(u, v) < inf  
                if label(v) > label(u) + w(u, v)  
                    label(v) = label(u) + w(u, v);  
                    f(v) = u;  % 更新前驱节点  
                end  
            end  
        end  
    end  

    % 获取最短路径和总距离  
    minDistance = label(terminal);  
    if isinf(minDistance)  
        path = [];  % 没有路径  
        return;  
    end  

    % 重建路径  
    path = [];  
    while terminal ~= start  
        path = [terminal, path];  % 将终点加入路径  
        terminal = f(terminal);    % 更新终点  
        if isnan(terminal)  % 防止无前驱导致的死循环  
            path = [];  % 没有路径  
            minDistance = inf;  
            return;  
        end  
    end  
    path = [start, path];  % 加入起点  
end  
