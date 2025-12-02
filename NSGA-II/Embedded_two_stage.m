function [BaseDemandC1,SatisfactionRate] = Embedded_two_stage(GAParameters,sorted_indices,DemandC1,childPath1,SupplyPoint,TransPoint,DemandPoint)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
pop = GAParameters.PopSize;
LL = length(SupplyPoint);
L1 = length(TransPoint);
L2 = length(DemandPoint);
DemandAll = sum(DemandC1,2); % 所有节点的需求
BaseDemandC1 = 0.2 * DemandC1; % 需求点的基础满足率 要加进去算需求
[BaseDemand1] = CaculateDemand(BaseDemandC1, pop, childPath1,SupplyPoint,TransPoint,DemandPoint);% 最低需求满足
BaseDemand = cell2mat(BaseDemand1);
[Demand11] = CaculateDemand(DemandC1, pop, childPath1,SupplyPoint,TransPoint,DemandPoint); % 正常情况下
Demand111 = cell2mat(Demand11);
BaseDemandB1 = BaseDemand(:,LL+1:LL+L1);
DemandB1 = Demand111(:,LL+1:LL+L1);
SupplyB1 = 0.8 * DemandB1 - BaseDemandB1; % B实际能供应满足分配的物资 
%     DemandA1 = Demand111(:, 1:LL); % 各个供应点的需求量
%     SupplyA1 = DemandA1*0.8; % 需求的供货量为0.8
RemainDemandC1 = DemandC1 - BaseDemandC1;
for t = 1:pop
    for k = 1:L2 %遍历每一个需求点
        p1 = find(sorted_indices == k); % 需求点索引
        tmp1 = DemandPoint(p1); % tmp1为需求点
        p2 = find(DemandPoint == tmp1); % 找到需求点在矩阵中对应的位置
        tmp2 = childPath1{t, L1 + p2}(1,1); % 找到该需求点的运送路径中的起点,即供应点
        p3 = find(TransPoint == tmp2); % 找到集散中心在矩阵中对应的位置
        if SupplyB1(t,p3) > 0 && SupplyB1(t,p3) -  RemainDemandC1(t,p2) > 0
            BaseDemandC1(t,p2) = BaseDemandC1(t,p2) + RemainDemandC1(t,p2);
            SupplyB1(t,p3) = SupplyB1(t,p3) -  RemainDemandC1(t,p2);
        elseif SupplyB1(t,p3) > 0 && SupplyB1(t,p3) -  RemainDemandC1(t,p2) < 0
            BaseDemandC1(t,p2) = BaseDemandC1(t,p2) + SupplyB1(t,p3);
            SupplyB1(t,p3) = 0;
        end
    end
end
SatisfactionRate = BaseDemandC1 ./ DemandC1;
end

