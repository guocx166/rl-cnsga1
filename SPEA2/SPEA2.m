%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SPEA2
%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SPEA2(UserInput,ProblemParameters,GAParameters)

%-----------------------------------------
%   Init individu vide
%-----------------------------------------
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
individuVide.S = [];
individuVide.R = [];
individuVide.sigma = [];
individuVide.sigmaK = [];
individuVide.D = [];
individuVide.F = [];

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
a1 = 190; a2 = 128; a3 = 192; a4 = 72;
SupplyPoint=[a1,a2,a3,a4];
A1 = point(a1,:); A2 = point(a2,:); A3 = point(a3,:); A4 = point(a4,:);
pointsSupply = {A1, A2, A3, A4};
% 中转点
b1 = 60; b2 = 10; b3 = 110;
TransPoint = [b1,b2,b3];
B1 = point(b1,:); B2 = point(b2,:);B3 = point(b3,:);
% 受灾点
c1 = 116; c2 = 124; c3 = 51; c4 = 87;
DemandPoint = [c1,c2,c3,c4];
C1 = point(c1,:); C2 = point(c2,:); C3 = point(c3,:); C4 = point(c4,:);
L1 = length(TransPoint);
L2 = length(DemandPoint);
L = L1 + L2;
rng(0);
DemandC1 = RandomDemand(pop,L2); % 为受灾节点的数量
DemandC2 = RandomDemand(pop,L2);
%-----------------------------------------
%  初始化种群
%-----------------------------------------
Pop1 = repmat(individuVide,N,1);
Pop2 = repmat(individuVide,N,1);
transit_nodes = setdiff(1:m, [a1,a2,a3,a4,b1,b2,b3,c1,c2,c3,c4]); % 可以选择的转运节点
[PopulationPaths1,PopulationPaths2,PopulationTrans1,PopulationTrans2] = generate(Dis,trans,trans_all,L1, L2, transit_nodes,pop,SupplyPoint,TransPoint,DemandPoint);
for i = 1:pop  
    newCellArrayPath1{i,1} = PopulationPaths1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    newCellArrayTran1{i,1} = PopulationTrans1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    newCellArrayPath2{i,1} = PopulationPaths2(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    newCellArrayTran2{i,1} = PopulationTrans2(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    Pop1(i).Val1 = newCellArrayPath1{i}; 
    Pop1(i).Val2 = newCellArrayTran1{i};
    Pop2(i).Val1 = newCellArrayPath2{i}; 
    Pop2(i).Val2 = newCellArrayTran2{i};
end
Pop1 = EvaluationNSGA2(UserInput,GAParameters,Pop1,Dis,DemandC1,SupplyPoint, TransPoint, DemandPoint);
Pop2 = EvaluationNSGA2(UserInput,GAParameters,Pop2,Dis,DemandC1,SupplyPoint, TransPoint, DemandPoint);
archive1 = [];
archive2 = [];
G = 1;

%-----------------------------------------
%   It閞ation pendant Gmax g閚閞ations
%-----------------------------------------
while (G < GAParameters.Gmax) 
    
    P1 = [Pop1; archive1];
    P2 = [Pop2; archive2];
    
    %-----------------------------------------
    %   Evaluation
    %-----------------------------------------
    P1 = EvaluationSPEA2(UserInput, GAParameters,P1,SupplyPoint,TransPoint,DemandPoint,Capacity,CostTrans,Dis,start_time,DemandC1);
    P2 = EvaluationSPEA2(UserInput, GAParameters,P2,SupplyPoint,TransPoint,DemandPoint,Capacity,CostTrans,Dis,start_time,DemandC2);
    %-----------------------------------------
    %   S閘ection de l'environnement
    %-----------------------------------------
    archive1 = SelectionEnvironnement(P1, GAParameters);
    archive2 = SelectionEnvironnement(P2, GAParameters);
    
    %-----------------------------------------
    %   Traitements
    %-----------------------------------------
    %Approximation du front de pareto
    paretoFront1 = archive1([archive1.F] < 1);
    paretoFront2 = archive2([archive2.F] < 1);
    
    %Plot R閟ultat
    AffichageResultats(paretoFront1, G);
    
    %V閞ification condition de fin
    if G >= GAParameters.Gmax
        break;
    end 
    
    %-----------------------------------------
    %   锦标赛选择
    %-----------------------------------------
    MP1 = SelectionTournoi(individuVide, UserInput.Algorithme, archive1);
    MP2 = SelectionTournoi(individuVide, UserInput.Algorithme, archive2);
    
    %-----------------------------------------
    %   交叉变异
    %-----------------------------------------
    % Enfants = SimulatedBinaryCrossover(individuVide, GAParameters, ProblemParameters, MP);
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
        [childPathOutput1,childTranOutput1] = cross(pop,prob1,prob2,PopulationPaths1,PopulationTrans1,trans,transport1,transport2,transport3,trans_all,Dis);
        [childPathOutput2,childTranOutput2] = cross(pop,prob1,prob2,PopulationPaths2,PopulationTrans2,trans,transport1,transport2,transport3,trans_all,Dis);
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
    % 两阶段优化模型的成本
    [T1,C1]=fitnessPopulation(pop,Dis,childPath1,childTran1,start_time);%计算适应度函数
    [T2,C2]=fitnessPopulation(pop,Dis,childPath2,childTran2,start_time);
    [Demand11] = CaculateDemand(DemandC1, pop, childPath1,SupplyPoint,TransPoint,DemandPoint);
    Demand111 = cell2mat(Demand11);
    [Demand22] = CaculateDemand(DemandC2, pop, childPath2,SupplyPoint,TransPoint,DemandPoint);
    Demand222 = cell2mat(Demand22);
    [CostTotal1,PopulationTransNum1] = updateCost(Capacity,CostTrans,childTran1,pop,Demand111(:,1:4),Demand111(:,5:8),Demand111(:,9:11));        
    [CostTotal2,PopulationTransNum2] = updateCost(Capacity,CostTrans,childTran2,pop,Demand222(:,1:4),Demand222(:,5:8),Demand222(:,9:11));
    Cost1 = sum(CostTotal1,2)+sum(C1,2);
    Time1 = sum(T1,2);
    Cost2 = sum(CostTotal2,2)+sum(C1,2);
    Time2 = sum(T2,2);  
    for i =1:pop
        Enfants1(i).ValObjective = [Time1(i);Cost1(i)];
        Enfants2(i).ValObjective = [Time2(i);Cost2(i)];
    end

    %-----------------------------------------
    %   Mutation (polynomial mutation)
    %-----------------------------------------
    % Pop = PolynomialMutation(UserInput.Algorithme,UserInput.Probleme, GAParameters, ProblemParameters, Enfants);
    Pop1 = Enfants1; % 变异后的个体
    Pop2 = Enfants2;
    G = G+1;
    AllObject = (cat(2,P1.ValObjective))';
    sumTime = AllObject(:,1);
    [sortTime, sortIndex1] = sort(sumTime);
    minTime = sortTime(1);
    minTimes(G,1) = minTime;
    fprintf('代数:%d   最短时间:%.2fh \n',G,minTime);
    sumCost = AllObject(:,2);
    [sortCost, sortIndex2] = sort(sumCost);
    minCost = sortCost(1);
    minCosts(G,1) = minCost;
    fprintf('代数:%d   最小成本:%.2f万元 \n',G,minCost);
end

%-----------------------------------------
%   Affichage r閟ultats finaux
%-----------------------------------------
AffichageResultats(paretoFront1, G);





