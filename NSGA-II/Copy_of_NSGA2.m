%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   NSGA2算法
%   这个是单阶段的鲁棒优化模型，直接是考虑的最差情况   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NSGA2(UserInput,ProblemParameters,GAParameters)
%%   初始化参数
N = GAParameters.PopSize;
gen = GAParameters.Gmax;
pop = GAParameters.PopSize;
prob1 = GAParameters.Pc; %temp
prob2 = GAParameters.Pm; %temp
start_time = GAParameters.StartTime;
Capacity = [50,200,20];   % 交通工具的容量
CostTrans = [0.02,0.2,2]; % 交通工具租借一辆的成本
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
data = dlmread('R1_2_1.TXT');
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
C1 = point(c1,:); C2 = point(c2,:); C3 = point(c3,:); C4 = point(c4,:);
L1 = length(TransPoint);
L2 = length(DemandPoint);
L = L1 + L2;
Lpop = L+1; % 把需求点也算到个体中后，个体的长度
rng(0);
DemandC1 = RandomDemand(pop,4); % 4可能需要修改，为受灾节点的数量

%% 初始种群
P = repmat(individuVide,N,1);
transit_nodes = setdiff(1:m, [a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,c4]); % 可以选择的转运节点
% transit_nodes = setdiff(1:m, [a1,a2,a3,b1,b2,c1,c2,c3]);
[PopulationPaths1,PopulationPaths2,PopulationTrans1,PopulationTrans2] = generate(Dis,trans,trans_all,L1, L2, transit_nodes,pop,SupplyPoint,TransPoint,DemandPoint);
for i = 1:pop  % 遍历每一行  
    newCellArrayPath{i,1} = PopulationPaths1(i, :); % 将整行（7 列）放入新的元胞数组的每一行  
end
for i = 1:pop  
    newCellArrayPath{i,1} = PopulationPaths1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    newCellArrayTran{i,1} = PopulationTrans1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    P(i).Val1 = newCellArrayPath{i}; 
    P(i).Val2 = newCellArrayTran{i};
end 
G = 1;
%% 计算适应值
P = EvaluationNSGA2(UserInput,GAParameters,P,Dis,DemandC1,SupplyPoint, TransPoint, DemandPoint);

%% 非支配排序
[P, F] = NonDominatedSorting(P);

%% 计算拥挤度
P = CrowdingDistance(P,F);

%% 筛选种群
[P, ~] = SortPopulation(P);

%% 开始迭代
%-----------------------------------------
while (G < GAParameters.Gmax)
% NSGA2流程    
    MP = SelectionTournoi(individuVide, UserInput.Algorithme, P);
% 交叉变异
    for i = 1:pop  
        newCellArray1{i} = MP(i).Val1; % 路径
        newCellArray2{i} = MP(i).Val2; % 交通工具
        PopulationPaths(i, :) = newCellArray1{i};  % 将每个单元格中的元胞内容放回原位 
        PopulationTrans(i, :) = newCellArray2{i};
    end 
    for p = 1:L
        [childPathOutput1,childTranOutput1] = cross(pop,prob1,prob2,PopulationPaths,PopulationTrans,trans,transport1,transport2,transport3,trans_all,Dis);
        childPath1(:,p) = childPathOutput1;
        childTran1(:,p) = childTranOutput1;
    end 
    Enfants = repmat(individuVide, N, 1);
    for i = 1:pop  
        newCellArrayPath{i,1} = childPath1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
        newCellArrayTran{i,1} = childTran1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
        Enfants(i).Val1 = newCellArrayPath{i}; 
        Enfants(i).Val2 = newCellArrayTran{i};
    end 
    % [t,c]=fitness(GAParameters,Dis,Enfants);
    % 计算适应值函数
    [T1,C1]=fitnessPopulation(pop,Dis,childPath1,childTran1,start_time);%计算适应度函数
    max_values1 = max(DemandC1);
    DemandRe1 = repmat(max_values1, size(C1, 1), 1);
    [tmp1] = CaculateDemand(DemandRe1, pop, childPath1,SupplyPoint,TransPoint,DemandPoint);
    DemandMax1 = cell2mat(tmp1);
    [CostTotal1,PopulationTransNum1] = updateCost(Capacity,CostTrans,childTran1,pop,DemandMax1(:,1:4),DemandMax1(:,5:8),DemandMax1(:,9:11));        
    Cost1 = sum(CostTotal1,2)+sum(C1,2);
    Time1 = max(T1,2);
    for i =1:pop
        Enfants(i).ValObjective = [Time1(i);Cost1(i)];
    end
    Mutants1 = Enfants; % 变异后的个体


    %-----------------------------------------
    %   Parents + Enfants mut閟
    %-----------------------------------------
    R = [P ; Mutants1];
    [R, F] = NonDominatedSorting(R);
    R = CrowdingDistance(R,F);
    [Q, ~] = SortPopulation(R);
    P = Q(1:GAParameters.PopSize);
    [P, F] = NonDominatedSorting(P);
    P = CrowdingDistance(P,F);
    %-----------------------------------------
    %   按 rank 和 distance 排序
    %-----------------------------------------
    [P, F] = SortPopulation(P);
    % Store F1
    F1 = P(F{1});
    %%  画图
    AffichageResultats(F1, G);
    sumTime=Time1;
    [sortTime, sortIndex1] = sort(sumTime);
    minTime = sortTime(1);
    minTimes(G,1) = minTime;
    fprintf('代数:%d   最短时间:%.2fh \n',G,minTime);
    sumCost=Cost1;
    [sortCost, sortIndex2] = sort(sumCost);
    minCost = sortCost(1);
    minCosts(G,1) = minCost;
    fprintf('代数:%d   最小成本:%.2f万元 \n',G,minCost);
    G = G+1;
end

%-----------------------------------------
% 帕累托前沿
%-----------------------------------------
AffichageResultats(F1, G);
% 画迭代图
figure 
plot(minTimes, 'MarkerFaceColor', 'red','LineWidth',1);
title('收敛曲线图（每一代的最短时间）');
set(gca,'ytick',1:1:10); 
ylabel('时间');
xlabel('迭代次数');
grid on
figure 
plot(minCosts, 'MarkerFaceColor', 'red','LineWidth',1);
title('收敛曲线图（每一代的最小成本）');
set(gca,'ytick',1:1:10); 
ylabel('成本');
xlabel('迭代次数');
grid on
