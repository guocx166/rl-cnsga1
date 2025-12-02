%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   计算适应值函数
%
%      
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = EvaluationNSGA2(UserInput,GAParameters,Po,Dis,DemandC1,SupplyPoint, TransPoint, DemandPoint)

%-----------------------------------------
%	Init variables
%-----------------------------------------
N = size(Po,1);
P = Po;
pop = GAParameters.PopSize;
start_time = GAParameters.StartTime;
L1 = length(TransPoint);
L2 = length(DemandPoint);
for i = 1:pop  
    newCellArray1{i} = P(i).Val1; % 路径
    newCellArray2{i} = P(i).Val2; % 交通工具
    PopulationPaths(i, :) = newCellArray1{i};  % 将每个单元格中的元胞内容放回原位 
    PopulationTrans(i, :) = newCellArray2{i};
end 
%-----------------------------------------
%	Affectation des valeurs de fitness
%-----------------------------------------
[T1,C1]=fitnessPopulation(pop,Dis,PopulationPaths,PopulationTrans,start_time);%计算适应度函数
% sumsT1 = sum(T1,2);
sumsT1 = max(T1(:, 1:L1), [], 2) + max(T1(:, (L1 + 1):L1+L2), [], 2);
sumsC1 = sum(C1,2);
max_values1 = max(DemandC1); 
DemandRe1 = repmat(max_values1, size(sumsC1, 1), 1);
[tmp1] = CaculateDemand(DemandRe1, pop, PopulationPaths, SupplyPoint, TransPoint, DemandPoint);
DemandMax1 = cell2mat(tmp1);
sumsC1 = sumsC1 + DemandMax1;

for i =1:pop
    P(i).ValObjective = [sumsT1(i);sumsC1(i)];
end

end