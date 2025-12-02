function [migP1, migP2] = migratePopulations(individuVide,DemandC1,DemandC2,GAParameters,L,Dis,P1, P2, m, Capacity,CostTrans,SupplyPoint, TransPoint, DemandPoint,start_time)
    % 确保输入的种群是元胞数组
    % 种群的迁移
%     if ~iscell(population1) || ~iscell(population2)
%         error('种群必须是元胞数组');
%     end

    % 获取种群的大小
pop = GAParameters.PopSize;
N = GAParameters.PopSize;
migP1 = repmat(individuVide,N,1);
migP2 = repmat(individuVide,N,1);
L1 = length(TransPoint);
L2 = length(DemandPoint);
LL = length(SupplyPoint);
for i = 1:pop  
    new11{i} = P1(i).Val1; % 路径
    new12{i} = P1(i).Val2; % 交通工具
    population1(i,1:L) = new11{i};
    population1(i,1+L:L+L) = new12{i};
%     Paths1(i, :) = new11{i};  % 将每个单元格中的元胞内容放回原位 
%     Trans1(i, :) = new12{i};
    new21{i} = P2(i).Val1; % 路径
    new22{i} = P2(i).Val2; % 交通工具
    population2(i,1:L) = new21{i};
    population2(i,1+L:L+L) = new22{i};
%     Paths2(i, :) = new21{i};  % 将每个单元格中的元胞内容放回原位 
%     Trans2(i, :) = new22{i};
end

pop_size1 = size(population1, 1);
pop_size2 = size(population2, 1);

% 检查m是否不超过两个种群的大小
if m > pop_size1 || m > pop_size2
    error('m不能大于两个种群的大小');
end

    % 从种群1中随机选择m个不重复的个体
idx1 = randperm(pop_size1, m);

    % 从种群2中随机选择m个不重复的个体
idx2 = randperm(pop_size2, m);

    % 交换选中的个体
for i = 1:m
    temp = population1(idx1(i),:); % 保存种群1中选中的个体
    population1(idx1(i),:) = population2(idx2(i),:); % 用种群2的个体替换
    population2(idx2(i),:) = temp; % 用种群1的个体替换
end
Paths1 = population1(:,1:L);
Trans1 = population1(:,L+1:L+L);
Paths2 = population2(:,1:L);
Trans2 = population2(:,L+1:L+L);
[T1,C1]=fitnessPopulation(pop,Dis,Paths1,Trans1,start_time);%计算适应度函数
[T2,C2]=fitnessPopulation(pop,Dis,Paths1,Trans1,start_time);
[Demand11] = CaculateDemand(DemandC1, pop, Paths1, SupplyPoint,TransPoint,DemandPoint);
Demand111 = cell2mat(Demand11);
[Demand22] = CaculateDemand(DemandC2, pop, Paths2, SupplyPoint,TransPoint,DemandPoint);
Demand222 = cell2mat(Demand22);
[CostTotal1,PopulationTransNum1] = updateCost(Capacity,CostTrans,Trans1,pop,Demand111(:,1:LL),Demand111(:,LL+1:LL+L1),Demand111(:,LL+L1+1:LL+L1+L2));        
[CostTotal2,PopulationTransNum2] = updateCost(Capacity,CostTrans,Trans2,pop,Demand222(:,1:LL),Demand222(:,LL+1:LL+L1),Demand222(:,LL+L1+1:LL+L1+L2));
Cost1 = sum(CostTotal1,2)+sum(C1,2);
% Time1 = sum(T1,2);
Time1 = max(T1(:, 1:L1), [], 2) + max(T1(:, (L1 + 1):L1+L2), [], 2);
Time2 = max(T2(:, 1:L1), [], 2) + max(T2(:, (L1 + 1):L1+L2), [], 2);
Cost2 = sum(CostTotal2,2)+sum(C1,2);
% Time2 = sum(T2,2);  
for i =1:pop
    ArrayPath1{i,1} = Paths1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    ArrayTran1{i,1} = Trans1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    ArrayPath2{i,1} = Paths2(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    ArrayTran2{i,1} = Trans2(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    migP1(i).Val1 = ArrayPath1{i}; 
    migP1(i).Val2 = ArrayTran1{i};
    migP2(i).Val1 = ArrayPath2{i}; 
    migP2(i).Val2 = ArrayTran2{i};
    migP1(i).ValObjective = [Time1(i);Cost1(i)];
    migP2(i).ValObjective = [Time2(i);Cost2(i)];
end
end

