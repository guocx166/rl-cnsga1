function [Demand11] = CaculateDemand(DemandC1, pop, chromosome1_tmp,SupplyPoint,TransPoint,DemandPoint)
% 各个节点的物资需求量
% 使用函数
L1 = length(TransPoint); % 3
L2 = length(DemandPoint); % 4
LL = length(SupplyPoint); % 4
DemandB1 = zeros(pop,L1);% 3为集散中心的个数
% b1,b2,b3为三个中转节点
for i = 1:pop
    for j = L1+1:L1+L2 % 4:7这里后面要改成符号比较好一些
        tmp1 = chromosome1_tmp{i,j}(1,1); % 找到集散中心
        tmp2 = chromosome1_tmp{i,j}(1,end); % 找到需求点
        tmp3 = DemandC1(i, DemandPoint == tmp2);
        DemandB1(i,TransPoint==tmp1) = DemandB1(i,TransPoint==tmp1) + tmp3;
    end
end
DemandA1 = zeros(pop,LL); % 4为供应点的个数
for i = 1:pop
    for j = 1:L1 % 这里后面要改成符号比较好一些
        tmp1 = chromosome1_tmp{i,j}(1,1); % 找到供应点
        tmp2 = chromosome1_tmp{i,j}(1,end); % 找到集散中心
        tmp3 = DemandB1(i, TransPoint == tmp2);
        DemandA1(i,SupplyPoint==tmp1) = DemandA1(i,SupplyPoint==tmp1)+tmp3;
    end
end
Demand1 = [DemandA1, DemandB1, DemandC1];
Demand11 = mat2cell(Demand1, ones(size(Demand1, 1), 1), size(Demand1, 2));
end

