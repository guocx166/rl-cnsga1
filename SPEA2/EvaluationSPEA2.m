%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   蓆ape d'関aluation
%
%       Affectation d'une 'qualit?' ? chaque
%       individu, une valeur de fitness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = EvaluationSPEA2(UserInput, GAParameters,P,SupplyPoint,TransPoint,DemandPoint,Capacity,CostTrans,Dis,start_time,DemandC1)

%-----------------------------------------
%	Init variables
%-----------------------------------------
N = GAParameters.PopSize;
NA = GAParameters.ArchiveSize;
NP = size(P,1);
dominance = false(NP,NP);
k = round(sqrt(N+NA));

%On initialise S ? 0
for i = 1:NP
    P(i).S = 0;
end

%-----------------------------------------
%	Test de domination => donne S
%-----------------------------------------
for i = 1:NP
   
    for j = i+1:NP
        
        %Si les variables de d閏ision de i dominent j
        if all(P(i).ValObjective <= P(j).ValObjective) && any(P(i).ValObjective < P(j).ValObjective)
            P(i).S = P(i).S + 1;
            dominance(i,j) = true;
            
        %Si les variables de d閏ision j dominent i
        elseif all(P(j).ValObjective <= P(i).ValObjective) && any(P(j).ValObjective < P(i).ValObjective)
            P(j).S = P(j).S + 1;
            dominance(j,i) = true;
        end
    end
end

%-----------------------------------------
%	Calcul de R = Raw Fitness
%-----------------------------------------
S = [P.S];
for i = 1:NP
    P(i).R = sum(S(dominance(:,i)));
end

%-----------------------------------------
%	Calcul de densit? D
%-----------------------------------------
LL = length(SupplyPoint);
L1 = length(TransPoint);
L2 = length(DemandPoint);
for i = 1:N  
    newCellArray11{i} = P(i).Val1; % 路径
    newCellArray12{i} = P(i).Val2; % 交通工具
    PopulationPaths1(i, :) = newCellArray11{i};  % 将每个单元格中的元胞内容放回原位 
    PopulationTrans1(i, :) = newCellArray12{i};
end
[T1,C1]=fitnessPopulation(N,Dis,PopulationPaths1,PopulationTrans1,start_time);%计算适应度函数
[Demand11] = CaculateDemand(DemandC1, N, PopulationPaths1,SupplyPoint,TransPoint,DemandPoint);
Demand111 = cell2mat(Demand11);
[CostTotal1,PopulationTransNum1] = updateCost(Capacity,CostTrans,PopulationTrans1,N,Demand111(:,1:LL),Demand111(:,LL+1:LL+L1),Demand111(:,LL+L1+1:LL+L1+L2));        
Cost1 = sum(CostTotal1,2)+sum(C1,2);
Time1 = sum(T1,2);
z = [];
for i = 1:N
    z = [z; Time1(i); Cost1(i)];
end
if NP > N
    for i = 1:N  
        newCellArray11{i} = P(i+N).Val1; % 路径
        newCellArray12{i} = P(i+N).Val2; % 交通工具
        PopulationPaths1(i, :) = newCellArray11{i};  % 将每个单元格中的元胞内容放回原位 
        PopulationTrans1(i, :) = newCellArray12{i};
    end
    [T1,C1]=fitnessPopulation(N,Dis,PopulationPaths1,PopulationTrans1,start_time);%计算适应度函数
    [Demand11] = CaculateDemand(DemandC1, N, PopulationPaths1,SupplyPoint,TransPoint,DemandPoint);
    Demand111 = cell2mat(Demand11);
    [CostTotal1,PopulationTransNum1] = updateCost(Capacity,CostTrans,PopulationTrans1,N,Demand111(:,1:LL),Demand111(:,LL+1:LL+L1),Demand111(:,LL+L1+1:LL+L1+L2));        
    Cost1 = sum(CostTotal1,2)+sum(C1,2);
    Time1 = sum(T1,2);
    for i = 1:N
        z = [z; Time1(i); Cost1(i)];
    end
end

%tableau de distance entre les individus, puis tri
sigma = pdist2(abs(z),abs(z)); 
sigma = sort(sigma);

for i = 1:NP
    P(i).sigma = sigma(:,i);
    P(i).sigmaK = P(i).sigma(k);
    
    %densit? D
    P(i).D = 1/(P(i).sigmaK+2);
end

%-----------------------------------------
%	Calcul de F
%-----------------------------------------
for i = 1:NP
    P(i).F = P(i).R + P(i).D;
end







