function [CostTotal,PopulationTransNum] = updateCost(Capacity,CostTrans,childTran1,pop,DemandA,DemandB,DemandC)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
% Capacity = [50,200,20]; % 交通工具的容量
% CostTrans = [0.02,0.2,2]; % 交通工具租借一辆的成本
[~,L1] = size(DemandB); % 3
[~,L2] = size(DemandC); % 4
[~,LL] = size(DemandA); % 4
PopulationTransNum = {};
CostTotal = [];
for i = 1:pop
    for j1 = 1:L1 % 3是集散中心的数量
        Trans1 = childTran1{i,j1};
        lt1 = nnz(Trans1);
        Cost1 = 0;
        TransNum1 = [];
        for k =1:lt1
            tmp = Trans1(k);
            Cost1 =  Cost1 + CostTrans(tmp) * ceil(DemandB(j1)/Capacity(tmp)); % 该种运输路线的成本
            TransNum1(k) =  ceil(DemandB(j1)/Capacity(tmp)); % 交通工具的使用数量
        end
        CostTotal(i,j1) = Cost1;
        PopulationTransNum{i,j1} = TransNum1;
    end
    for j2 = L1+1:L1+L2 % 4是受灾点的数量
        Trans2 = childTran1{i,j2};
        lt2 = nnz(Trans2);
        Cost2 = 0;
        TransNum2 = [];
        for k =1:lt2
            tmp = Trans2(k);
            Cost2 =  Cost2 + CostTrans(tmp) * ceil(DemandC(j2)/Capacity(tmp)); % 该种运输路线的成本
            TransNum2(k) =  ceil(DemandC(j2)/Capacity(tmp)); % 交通工具的使用数量
        end
        CostTotal(i,j2) = Cost2;
        PopulationTransNum{i,j2} = TransNum2;
    end 
end
end

