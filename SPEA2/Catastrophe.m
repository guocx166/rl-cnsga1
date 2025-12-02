function [P1,P2] = Catastrophe(P1,P2,GAParameters,filename1,filename2,trans_all,trans,Dis, L, Capacity,CostTrans,DemandC1,DemandC2,SupplyPoint,TransPoint,DemandPoint,start_time)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
pop = GAParameters.PopSize;
for i = 1:pop  
    new11{i} = P1(i).Val1; % 路径
    new12{i} = P1(i).Val2; % 交通工具
    Paths1(i, :) = new11{i};  % 将每个单元格中的元胞内容放回原位 
    Trans1(i, :) = new12{i};
    new21{i} = P2(i).Val1; % 路径
    new22{i} = P2(i).Val2; % 交通工具
    Paths2(i, :) = new21{i};  % 将每个单元格中的元胞内容放回原位 
    Trans2(i, :) = new22{i};
end
[ChangePopulationPaths1,ChangePopulationPaths2,ChangePopulationTrans1,ChangePopulationTrans2] = generatePop(filename1,filename2,trans_all,trans,Dis, L, pop);
    % 替换原种群中的个体  
    % 随机选择n个个体的索引  
selectedIndices = randperm(pop, pop/2);  
for si = 1:pop/2  
    Paths1(selectedIndices(si),:) = ChangePopulationPaths1(si,:); % 使用{}操作符替换单元格内容 
    Trans1(selectedIndices(si),:) = ChangePopulationTrans1(si,:);
    Paths2(selectedIndices(si),:) = ChangePopulationPaths2(si,:); % 使用{}操作符替换单元格内容 
    Trans2(selectedIndices(si),:) = ChangePopulationTrans2(si,:);
end
[T1,C1]=fitnessPopulation(pop,Dis,Paths1,Trans1,start_time);%计算适应度函数
[T2,C2]=fitnessPopulation(pop,Dis,Paths2,Trans2,start_time);
[Demand11] = CaculateDemand(DemandC1, pop, Paths1,SupplyPoint,TransPoint,DemandPoint);
Demand111 = cell2mat(Demand11);
[Demand22] = CaculateDemand(DemandC2, pop, Paths2,SupplyPoint,TransPoint,DemandPoint);
Demand222 = cell2mat(Demand22);
[CostTotal1,PopulationTransNum1] = updateCost(Capacity,CostTrans,Trans1,pop,Demand111(:,1:4),Demand111(:,5:8),Demand111(:,9:11));        
[CostTotal2,PopulationTransNum2] = updateCost(Capacity,CostTrans,Trans2,pop,Demand222(:,1:4),Demand222(:,5:8),Demand222(:,9:11));
Cost1 = sum(CostTotal1,2)+sum(C1,2);
Time1 = sum(T1,2);
Cost2 = sum(CostTotal2,2)+sum(C1,2);
Time2 = sum(T2,2);  
for i =1:pop
    ArrayPath1{i,1} = Paths1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    ArrayTran1{i,1} = Trans1(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    ArrayPath2{i,1} = Paths2(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    ArrayTran2{i,1} = Trans2(i, :); % 将整行（7 列）放入新的元胞数组的每一行
    P1(i).Val1 = ArrayPath1{i}; 
    P1(i).Val2 = ArrayTran1{i};
    P2(i).Val1 = ArrayPath2{i}; 
    P2(i).Val2 = ArrayTran2{i};    
    P1(i).ValObjective = [Time1(i);Cost1(i)];
    P2(i).ValObjective = [Time2(i);Cost2(i)];
end
end

