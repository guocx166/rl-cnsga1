%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   NSGA2算法
%   这个写双阶段的鲁棒优化模型   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NSGA2(UserInput,ProblemParameters,GAParameters)
%%   初始化参数
N = GAParameters.PopSize;
gen = GAParameters.Gmax;
pop = GAParameters.PopSize;
prob1 = GAParameters.Pc; %temp
prob2 = GAParameters.Pm; %temp
start_time = GAParameters.StartTime;
discount_factor = GAParameters.discount_factor; % 折扣因子
learning_rate = GAParameters.learning_rate;   % 学习率
epsilon = GAParameters.epsilon;         % 探索率
num_states = GAParameters.num_states;        % 状态数量
num_actions = GAParameters.num_actions;       % 动作数量 (+0.02, 0, -0.02)
Capacity = [50,200,20];   % 交通工具的容量
CostTrans = [0.02,0.2,2]; % 交通工具租借一辆的成本
Q_table = zeros(num_states, num_actions);
mt = 0;
% individuVide.Val = [];
individuVide.Val1 = [];
individuVide.Val2 = [];
individuVide.ValObjective = [];
individuVide.Rank = [];
individuVide.CrowdingDistance = [];
individuVide.DominationSet = [];
individuVide.DominatedCount = [];
%记录迭代
minTimes=zeros(gen,1);
minCosts=zeros(gen,1);
% 读取数据并处理数据
data = dlmread('r50.TXT');
[point] = readData(data);
[m,~] = size(point);
[transport1,transport2,transport3] = generateTransport(m);
[trans,trans_all] = readTransport(transport1,transport2,transport3); % tran表示节点是否联通，trans_all表示节点间可选择的交通工具
Dis = distance(point);
for i=1:m
    for j=1:m
        if trans(i,j)==0
            Dis(i,j)=inf;
        end
    end
end
% 供应点坐标
% a1 = 50; a2 = 46; a3 = 39;
% SupplyPoint=[a1,a2,a3];
a1 = 50; a2 = 47; a3 = 18; a4 = 39;
SupplyPoint=[a1,a2,a3,a4];
% 中转点
% b1 = 32; b2 = 14;
% TransPoint = [b1,b2];
b1 = 32; b2 = 1; b3 = 3;
TransPoint = [b1,b2,b3];
% 受灾点
% c1 = 35; c2 = 25; c3 = 40;
% DemandPoint = [c1,c2,c3];
c1 = 35; c2 = 25; c3 = 40; c4 = 39;
DemandPoint = [c1,c2,c3,c4];
L1 = length(TransPoint);
L2 = length(DemandPoint);
LL = length(SupplyPoint);
L = L1 + L2;
rng(0);
DemandC1 = RandomDemand(pop,L2); % 为受灾节点的数量
DemandC2 = RandomDemand(pop,L2);
filename = 'ftopsis.xls';
[sortedDict,sorted_indices] = ranknum(filename,DemandPoint);
%% 初始种群
P1 = repmat(individuVide,N,1);
P2 = repmat(individuVide,N,1);
transit_nodes = setdiff(1:m, [a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,c4]); % 可以选择的转运节点
% transit_nodes = setdiff(1:m, [a1,a2,a3,b1,b2,c1,c2,c3]);
[PopulationPaths1,PopulationPaths2,PopulationTrans1,PopulationTrans2] = generate(Dis,trans,trans_all,L1, L2, transit_nodes,pop,SupplyPoint,TransPoint,DemandPoint);
for i = 1:pop  
    newCellArrayPath1{i,1} = PopulationPaths1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    newCellArrayTran1{i,1} = PopulationTrans1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    newCellArrayPath2{i,1} = PopulationPaths2(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    newCellArrayTran2{i,1} = PopulationTrans2(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    P1(i).Val1 = newCellArrayPath1{i}; 
    P1(i).Val2 = newCellArrayTran1{i};
    P2(i).Val1 = newCellArrayPath2{i}; 
    P2(i).Val2 = newCellArrayTran2{i};
end 
G = 1;
%% 计算适应值
P1 = EvaluationNSGA2(UserInput,GAParameters,P1,Dis,DemandC1,SupplyPoint, TransPoint, DemandPoint);
P2 = EvaluationNSGA2(UserInput,GAParameters,P2,Dis,DemandC1,SupplyPoint, TransPoint, DemandPoint);
%% 非支配排序
[P1, F1] = NonDominatedSorting(P1);
[P2, F2] = NonDominatedSorting(P2);
%% 计算拥挤度
P1 = CrowdingDistance(P1,F1);
P2 = CrowdingDistance(P2,F2);
%% 筛选种群
[P1, ~] = SortPopulation(P1);
[P2, ~] = SortPopulation(P2);
tmp1 = cat(2,P1.ValObjective);
tmp2 = cat(2,P2.ValObjective);
phi_10 = calculateDiversity(tmp1');
phi_20 = calculateDiversity(tmp2');
%% 开始迭代
%-----------------------------------------
GlobalminTimes = zeros(pop,1);
GlobalminCosts = zeros(pop,1);
GlobalmaxSatisfactionRate = zeros(pop,1);
while (G < GAParameters.Gmax)
% NSGA2流程    
    MP1 = SelectionTournoi(individuVide, UserInput.Algorithme, P1);
    MP2 = SelectionTournoi(individuVide, UserInput.Algorithme, P2);
% 交叉变异
    for i = 1:pop  
        newCellArray11{i} = MP1(i).Val1; % 路径
        newCellArray12{i} = MP1(i).Val2; % 交通工具
        PopulationPaths1(i, :) = newCellArray11{i};  % 将每个单元格中的元胞内容放回原位 
        PopulationTrans1(i, :) = newCellArray12{i};
        newCellArray21{i} = MP2(i).Val1; % 路径
        newCellArray22{i} = MP2(i).Val2; % 交通工具
        PopulationPaths2(i, :) = newCellArray21{i};  % 将每个单元格中的元胞内容放回原位 
        PopulationTrans2(i, :) = newCellArray22{i};
    end 
    for p = 1:L
        [childPathOutput1,childTranOutput1] = cross(pop,prob1,prob2,PopulationPaths1(:,p),PopulationTrans1(:,p),trans,transport1,transport2,transport3,trans_all,Dis);
        [childPathOutput2,childTranOutput2] = cross(pop,prob1,prob2,PopulationPaths2(:,p),PopulationTrans2(:,p),trans,transport1,transport2,transport3,trans_all,Dis);
        childPath1(:,p) = childPathOutput1;
        childTran1(:,p) = childTranOutput1;
        childPath2(:,p) = childPathOutput2;
        childTran2(:,p) = childTranOutput2;
    end
    Enfants1 = repmat(individuVide, N, 1);
    Enfants2 = repmat(individuVide, N, 1);
    for i = 1:pop  
        newCellArrayPath1{i,1} = childPath1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
        newCellArrayTran1{i,1} = childTran1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
        newCellArrayPath2{i,1} = childPath2(i, :); % 将整行（7 列）放入新的元胞数组的每一行
        newCellArrayTran2{i,1} = childTran2(i, :); % 将整行（7 列）放入新的元胞数组的每一行
        Enfants1(i).Val1 = newCellArrayPath1{i}; 
        Enfants1(i).Val2 = newCellArrayTran1{i};
        Enfants2(i).Val1 = newCellArrayPath2{i}; 
        Enfants2(i).Val2 = newCellArrayTran2{i};
    end
    %%
    [T1,C1]=fitnessPopulation(pop,Dis,childPath1,childTran1,start_time);%计算适应度函数
    [T2,C2]=fitnessPopulation(pop,Dis,childPath2,childTran2,start_time);
    [BaseDemandC1,SatisfactionRate1] = Embedded_two_stage(GAParameters,sorted_indices,DemandC1,childPath1,SupplyPoint,TransPoint,DemandPoint);
    [BaseDemandC2,SatisfactionRate2] = Embedded_two_stage(GAParameters,sorted_indices,DemandC1,childPath2,SupplyPoint,TransPoint,DemandPoint);
%     DemandAll = sum(DemandC1,2); % 所有节点的需求
%     BaseDemandC1 = 0.2 * DemandC1; % 需求点的基础满足率 要加进去算需求
%     [BaseDemand1] = CaculateDemand(BaseDemandC1, pop, childPath1,SupplyPoint,TransPoint,DemandPoint);% 最低需求满足
%     BaseDemand = cell2mat(BaseDemand1);
%     [Demand11] = CaculateDemand(DemandC1, pop, childPath1,SupplyPoint,TransPoint,DemandPoint); % 正常情况下
%     Demand111 = cell2mat(Demand11);
%     BaseDemandB1 = BaseDemand(:,LL+1:LL+L1);
%     DemandB1 = Demand111(:,LL+1:LL+L1);
%     SupplyB1 = 0.8 * DemandB1 - BaseDemandB1; % B实际能供应满足分配的物资 
% %     DemandA1 = Demand111(:, 1:LL); % 各个供应点的需求量
% %     SupplyA1 = DemandA1*0.8; % 需求的供货量为0.8
%     RemainDemandC1 = DemandC1 - BaseDemandC1;
%     for t = 1:pop
%         for k = 1:L2 %遍历每一个需求点
%             p1 = find(sorted_indices == k); % 需求点索引
%             tmp1 = DemandPoint(p1); % tmp1为需求点
%             p2 = find(DemandPoint == tmp1); % 找到需求点在矩阵中对应的位置
%             tmp2 = childPath1{t, L1 + p2}(1,1); % 找到该需求点的运送路径中的起点,即供应点
%             p3 = find(TransPoint == tmp2); % 找到集散中心在矩阵中对应的位置
%             if SupplyB1(t,p3) > 0 && SupplyB1(t,p3) -  RemainDemandC1(t,p2) > 0
%                 BaseDemandC1(t,p2) = BaseDemandC1(t,p2) + RemainDemandC1(t,p2);
%                 SupplyB1(t,p3) = SupplyB1(t,p3) -  RemainDemandC1(t,p2);
%             elseif SupplyB1(t,p3) > 0 && SupplyB1(t,p3) -  RemainDemandC1(t,p2) < 0
%                 BaseDemandC1(t,p2) = BaseDemandC1(t,p2) + SupplyB1(t,p3);
%                 SupplyB1(t,p3) = 0;
%             end
%         end
%     end
%     SatisfactionRate = BaseDemandC1 ./ DemandC1;
    punishCost1 = L2 - sum(SatisfactionRate1,2) ;
    punishCost2 = L2 - sum(SatisfactionRate2,2) ;
    [ChangeDemand1] = CaculateDemand(BaseDemandC1, pop, childPath1,SupplyPoint,TransPoint,DemandPoint); 
    [ChangeDemand2] = CaculateDemand(BaseDemandC2, pop, childPath2,SupplyPoint,TransPoint,DemandPoint); 
    ChangeDemand1 = cell2mat(ChangeDemand1);
    ChangeDemand2 = cell2mat(ChangeDemand2);
    [CostTotal1,PopulationTransNum1] = updateCost(Capacity,CostTrans,childTran1,pop,ChangeDemand1(:,1:LL),ChangeDemand1(:,LL+1:LL+L1),ChangeDemand1(:,LL+L1+1:LL+L1+L2));        
    Cost1 = (1+punishCost1).*(sum(CostTotal1,2)+sum(C1,2));
%     Time1 = (1+punishCost1).*(max(T1(:, 1:L1), [], 2) + max(T1(:, (L1 + 1):L1+L2), [], 2));
    Time1 = max(T1(:, (L1 + 1):L1+L2), [], 2)+ max(T1(:, 1:L1), [], 2);
    [CostTotal2,PopulationTransNum2] = updateCost(Capacity,CostTrans,childTran2,pop,ChangeDemand2(:,1:LL),ChangeDemand2(:,LL+1:LL+L1),ChangeDemand2(:,LL+L1+1:LL+L1+L2));        
    Cost2 = (1+punishCost2).*(sum(CostTotal2,2)+sum(C2,2));
%     Time2 = (1+punishCost2).*(max(T2(:, 1:L1), [], 2) + max(T2(:, (L1 + 1):L1+L2), [], 2));
    Time2 = max(T2(:, (L1 + 1):L1+L2), [], 2)+max(T2(:, 1:L1), [], 2);
    %%
%     % 两阶段优化模型的成本
%     [T1,C1]=fitnessPopulation(pop,Dis,childPath1,childTran1,start_time);%计算适应度函数
%     [T2,C2]=fitnessPopulation(pop,Dis,childPath2,childTran2,start_time);
%     [Demand11] = CaculateDemand(DemandC1, pop, childPath1,SupplyPoint,TransPoint,DemandPoint);
%     Demand111 = cell2mat(Demand11);
%     [Demand22] = CaculateDemand(DemandC2, pop, childPath2,SupplyPoint,TransPoint,DemandPoint);
%     Demand222 = cell2mat(Demand22);
%     [CostTotal1,PopulationTransNum1] = updateCost(Capacity,CostTrans,childTran1,pop,Demand111(:,1:LL),Demand111(:,LL+1:LL+L1),Demand111(:,LL+L1+1:LL+L1+L2));        
%     [CostTotal2,PopulationTransNum2] = updateCost(Capacity,CostTrans,childTran2,pop,Demand222(:,1:LL),Demand222(:,LL+1:LL+L1),Demand222(:,LL+L1+1:LL+L1+L2));
%     Cost1 = sum(CostTotal1,2)+sum(C1,2);
%     Time1 = max(T1(:, 1:L1), [], 2) + max(T1(:, (L1 + 1):L1+L2), [], 2);
%     Time2 = max(T2(:, 1:L1), [], 2) + max(T2(:, (L1 + 1):L1+L2), [], 2);
% %     Time1 = sum(T1,2);
%     Cost2 = sum(CostTotal2,2)+sum(C1,2);
% %     Time2 = sum(T2,2);  
    for i =1:pop
        Enfants1(i).ValObjective = [Time1(i);Cost1(i)];
        Enfants2(i).ValObjective = [Time2(i);Cost2(i)];
    end
    Mutants1 = Enfants1; % 变异后的个体
    Mutants2 = Enfants2;
    %-----------------------------------------
    %   Parents + Enfants mut閟
    %-----------------------------------------
    R1 = [P1 ; Mutants1];
    [R1, F1] = NonDominatedSorting(R1);
    R1 = CrowdingDistance(R1,F1);
    [Q1, ~] = SortPopulation(R1);
    P1 = Q1(1:GAParameters.PopSize);
    [P1, F1] = NonDominatedSorting(P1);
    P1 = CrowdingDistance(P1,F1);
    [P1, F1] = SortPopulation(P1);
    F1 = P1(F1{1});
    %
    R2 = [P2 ; Mutants2];
    [R2, F2] = NonDominatedSorting(R2);
    R2 = CrowdingDistance(R2,F2);
    [Q2, ~] = SortPopulation(R2);
    P2 = Q2(1:GAParameters.PopSize);
    [P2, F2] = NonDominatedSorting(P2);
    P2 = CrowdingDistance(P2,F2);
    [P2, F2] = SortPopulation(P2);
    F2 = P2(F2{1});
    phi_1_current = calculateDiversity((cat(2,P1.ValObjective))');
    phi_2_current = calculateDiversity((cat(2,P2.ValObjective))');
    %% 灾变策略
    if ~mod(G,30)
        % 灾变策略时用于替换个体
        % 指定CSV文件名  
%         filename1 = 'path1.csv';
%         filename2 = 'path2.csv';
%         [P1,P2] = Catastrophe(P1,P2,GAParameters,filename1,filename2,trans_all,trans,Dis, L, Capacity,CostTrans,DemandC1,DemandC2,SupplyPoint,TransPoint,DemandPoint,start_time);
        [P1,P2] = Catastrophe(P1,P2,GAParameters,transit_nodes,trans_all,trans,Dis, Capacity,CostTrans,DemandC1,DemandC2,SupplyPoint,TransPoint,DemandPoint,start_time); 
        R1 = [P1 ; Mutants1];
        [R1, F1] = NonDominatedSorting(R1);
        R1 = CrowdingDistance(R1,F1);
        [Q1, ~] = SortPopulation(R1);
        P1 = Q1(1:GAParameters.PopSize);
        [P1, F1] = NonDominatedSorting(P1);
        P1 = CrowdingDistance(P1,F1);
        [P1, F1] = SortPopulation(P1);
        F1 = P1(F1{1});
        %
        R2 = [P2 ; Mutants2];
        [R2, F2] = NonDominatedSorting(R2);
        R2 = CrowdingDistance(R2,F2);
        [Q2, ~] = SortPopulation(R2);
        P2 = Q2(1:GAParameters.PopSize);
        [P2, F2] = NonDominatedSorting(P2);
        P2 = CrowdingDistance(P2,F2);
        [P2, F2] = SortPopulation(P2);
        F2 = P2(F2{1});
    end
    %% 种群迁移，强化学习
    % 强化学习循环，跟主体一起循环
    % 计算状态
    state = getState(phi_10, phi_20, phi_1_current, phi_2_current);
    % 选取动作 (此处简单实现随机选取)
    action = epsilonGreedy(Q_table, state, epsilon);
    A_values = [0.02, 0, -0.02];
    A = A_values(action);
    % 更新迁徙参数
    mt = mt + A * pop;
    
    % 计算奖赏
    reward1 = computeReward(phi_1_current, phi_10);
    reward2 = computeReward(phi_2_current, phi_20);
    reward = reward1 + reward2;
    if mt > 0
        mtt = round(mt);
        [migP1, migP2] = migratePopulations(individuVide,DemandC1,DemandC2,GAParameters,L,Dis,P1, P2, mtt, Capacity,CostTrans,SupplyPoint, TransPoint, DemandPoint,start_time);
    end
    % 重新计算种群的多样性
    phi_1_next = calculateDiversity((cat(2,migP1.ValObjective))');
    phi_2_next = calculateDiversity((cat(2,migP1.ValObjective))');
    % 计算下一个状态（此处用于演示，通常应在实际环境中获取）
    next_state = getState(phi_10, phi_20, phi_1_next, phi_2_next);
    best_next_action = max(Q_table(next_state, :));
    Q_table(state, action) = Q_table(state, action) + ...
    learning_rate * (reward + discount_factor * best_next_action - Q_table(state, action));
    % 找到最优m对应的动作
    [~, best_action] = max(Q_table(state, :));
    optimal_m = mt + A_values(best_action) * pop;
    % 种群迁移 
    if ~mod(G,20)
        P1 = migP1;
        P2 = migP2;
        R1 = [P1 ; Mutants1];
        [R1, F1] = NonDominatedSorting(R1);
        R1 = CrowdingDistance(R1,F1);
        [Q1, ~] = SortPopulation(R1);
        P1 = Q1(1:GAParameters.PopSize);
        [P1, F1] = NonDominatedSorting(P1);
        P1 = CrowdingDistance(P1,F1);
        [P1, F1] = SortPopulation(P1);
        F1 = P1(F1{1});
        %
        R2 = [P2 ; Mutants2];
        [R2, F2] = NonDominatedSorting(R2);
        R2 = CrowdingDistance(R2,F2);
        [Q2, ~] = SortPopulation(R2);
        P2 = Q2(1:GAParameters.PopSize);
        [P2, F2] = NonDominatedSorting(P2);
        P2 = CrowdingDistance(P2,F2);
        [P2, F2] = SortPopulation(P2);
        F2 = P2(F2{1});
    end
    %%  画图
    AffichageResultats(F1, G);
    AllObject = (cat(2,P1.ValObjective))';
    sumTime = AllObject(:,1);
    [sortTime, sortIndex1] = sort(sumTime);
    minTime = sortTime(1);
    minTimes(G,1) = minTime;
%     fprintf('代数:%d   最短时间:%.2fh \n',G,minTime);
    sumCost = AllObject(:,2);
    [sortCost, sortIndex2] = sort(sumCost);
    minCost = sortCost(1);
    minCosts(G,1) = minCost;
    SatisfactionRate = sort(SatisfactionRate1,1);
    minSatisfactionRate(G,:) = SatisfactionRate(1,:);
    maxS = SatisfactionRate(pop,:);
    totalmaxSatisfactionRate = sum(maxS,2)/L2;
%     fprintf('代数:%d   最小成本:%.2f万元 \n',G,minCost);
    if G > 1
        GlobalminTimes(G,1) = min(minTime,GlobalminTimes(G-1,1));
        GlobalminCosts(G,1) = min(minCost,GlobalminCosts(G-1,1));
        GlobalmaxSatisfactionRate(G,1) = max(totalmaxSatisfactionRate,GlobalmaxSatisfactionRate(G-1,1));
    else  
        GlobalminTimes(1,1) = minTime;
        GlobalminCosts(1,1) = minCost;
        GlobalmaxSatisfactionRate(G,1) = totalmaxSatisfactionRate;
    end
    fprintf('代数:%d   最短时间:%.2fh \n',G,GlobalminTimes(G,1));
    fprintf('代数:%d   最小成本:%.2f万元 \n',G, GlobalminCosts(G,1));

    
    G = G+1;
    
    
    
end
%-----------------------------------------
% 帕累托前沿
%-----------------------------------------
AffichageResultats(F1, G);
% 画迭代图
figure 
% plot(minTimes, 'MarkerFaceColor', 'red','LineWidth',1);
hold on;
plot(GlobalminTimes,'MarkerFaceColor', 'red','LineWidth',1)
title('收敛曲线图（每一代的最短时间）');
set(gca,'ytick',1:1:10); 
ylabel('时间');
xlabel('迭代次数');
grid on
figure 
% plot(minCosts, 'MarkerFaceColor', 'red','LineWidth',1);
hold on;
plot(GlobalminCosts,'MarkerFaceColor', 'red','LineWidth',1)
title('收敛曲线图（每一代的最小成本）');
set(gca,'ytick',1:1:10); 
ylabel('成本');
xlabel('迭代次数');
grid on
figure 
% plot(minCosts, 'MarkerFaceColor', 'red','LineWidth',1);
hold on;
plot(GlobalmaxSatisfactionRate,'MarkerFaceColor', 'red','LineWidth',1)
title('收敛曲线图（每一代的最大满足率）');
set(gca,'ytick',1:1:10); 
ylabel('物资的满足率');
xlabel('迭代次数');
grid on

k = size(minSatisfactionRate, 2);  

% 创建一个新的图形窗口  
figure;  

% 循环遍历每一列  
for i = 1:k  
    subplot(2, 2, i);  % 将画布划分为 1 行 k 列，激活第 i 个子图  
%     plot(minSatisfactionRate(:, i), 'o-'); % 绘制第 i 列的数据，使用图形标记和线  
    scatter(1:size(minSatisfactionRate, 1), minSatisfactionRate(:, i), 'filled'); % 绘制散点图 
    title(['Column ' num2str(i)]); % 为每个子图添加标题  
    xlabel('Row Index'); % x 轴标签  
    ylabel(['Value of Column ' num2str(i)]); % y 轴标签
    ylim([0, 1]); % 统一 Y 轴范围
end  

% 设置整体标题  
annotation('textbox', [0, 0.95, 1, 0.05], 'String', 'Plots of Each Column in a Matrix', ...  
           'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14); 