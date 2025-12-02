function [PopulationPaths1,PopulationPaths2,PopulationTrans1,PopulationTrans2] = generate(Dis,trans,trans_all,L1, L2, transit_nodes,pop,SupplyPoint,TransPoint,DemandPoint)
% 随机生成初始种群，使用函数
L = L1 + L2;
PopulationPaths1 = cell(L, pop);
PopulationPaths2 = cell(L,pop);
random_numbers = randi([7, 12], L, pop);
% 种群1
for i = 1:L1
    s = SupplyPoint(randi([1, length(SupplyPoint)]));
    t = TransPoint(i);
    for j = 1:pop
        r = random_numbers(i,j);
        selected_nodes = datasample(transit_nodes, r, 'Replace', false); % 或者使用 randperm 选择
        path = [s,selected_nodes,t];
        ll = nnz(path); % 求其长度
        pathnew = [path(1)];
        pointend = path(ll);
        for k = 1:ll-1
            a = path(k); b = path(k+1);
            if trans(a,b) == 1
                pathnew = [pathnew b];
            elseif trans(a, b) == 0
                pathadd = dijkstra1(Dis, a, b, pointend);
                pathnew = [pathnew pathadd(2:end)];
            end
        end
        transnew = generateTrans(trans_all,pathnew);
        PopulationPaths1{i,j} = pathnew;
        PopulationTrans1{i,j} = transnew;
    end
end
for i = 1:L2
    s = TransPoint(randi([1, length(TransPoint)]));
    t = DemandPoint(i);
    for j = 1:pop
        r = random_numbers(i+L1,j);
        selected_nodes = datasample(transit_nodes, r, 'Replace', false); % 或者使用 randperm 选择
        path = [s,selected_nodes,t];
        ll = nnz(path); % 求其长度
        pathnew = [path(1)];
        pointend = path(ll);
        for k = 1:ll-1
            a = path(k); b = path(k+1);
            if trans(a,b) == 1
                pathnew = [pathnew b];
            elseif trans(a, b) == 0
                pathadd = dijkstra1(Dis, a, b, pointend);
                pathnew = [pathnew pathadd(2:end)];
            end
        end
        transnew = generateTrans(trans_all,pathnew);
        PopulationPaths1{i+L1,j} = pathnew;
        PopulationTrans1{i+L1,j} = transnew;
    end
end
% 种群2
for i = 1:L1
    s = SupplyPoint(randi([1, length(SupplyPoint)]));
    t = TransPoint(i);
    for j = 1:pop
        r = random_numbers(i,j);
        selected_nodes = datasample(transit_nodes, r, 'Replace', false); % 或者使用 randperm 选择
        path = [s,selected_nodes,t];
        ll = nnz(path); % 求其长度
        pathnew = [path(1)];
        pointend = path(ll);
        for k = 1:ll-1
            a = path(k); b = path(k+1);
            if trans(a,b) == 1
                pathnew = [pathnew b];
            elseif trans(a, b) == 0
                pathadd = dijkstra1(Dis, a, b, pointend);
                pathnew = [pathnew pathadd(2:end)];
            end
        end
        transnew = generateTrans(trans_all,pathnew);
        PopulationPaths2{i,j} = pathnew;
        PopulationTrans2{i,j} = transnew;
    end
end
for i = 1:L2
    s = TransPoint(randi([1, length(TransPoint)]));
    t = DemandPoint(i);
    for j = 1:pop
        r = random_numbers(i+L1,j);
        selected_nodes = datasample(transit_nodes, r, 'Replace', false); % 或者使用 randperm 选择
        path = [s,selected_nodes,t];
        ll = nnz(path); % 求其长度
        pathnew = [path(1)];
        pointend = path(ll);
        for k = 1:ll-1
            a = path(k); b = path(k+1);
            if trans(a,b) == 1
                pathnew = [pathnew b];
            elseif trans(a, b) == 0
                pathadd = dijkstra1(Dis, a, b, pointend);
                pathnew = [pathnew pathadd(2:end)];
            end
        end
        transnew = generateTrans(trans_all,pathnew);
        PopulationPaths2{i+L1,j} = pathnew;
        PopulationTrans2{i+L1,j} = transnew;
    end
end
PopulationPaths1 = PopulationPaths1';
PopulationTrans1 = PopulationTrans1';
PopulationPaths2 = PopulationPaths2';
PopulationTrans2 = PopulationTrans2';
end