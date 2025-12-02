function [T,C] = fitnessPopulation(pop,Dis,PopulationPaths1,PopulationTrans1,start_time)
% fitnessPopulation适应度函数，包括时间和成本两项
% 成本包括运输成本和转运成本；时间包括运输时间、转运时间以及等待时间,chromosome1为路径，chromosome2为交通工具染色体，
% 此函数中有交通工具运行速度、转运时间、火车发车的时间；百公里运输成本、转运成本变量。
% 1为公路，2为火车，3为飞机
% 使用函数
%% 输入各种参数
v=[90,60,180];%速度,分别为汽车、火车和飞机
%start_time=9.00;%运输开始时间，作为输入变量
ex_time=[0,1,1.5;
    1,0,2;
    1.5,2,0];
T=0;%交通工具的运行时间
t1=[8.00,10.00,12.00,14.00,16.00,18.00];%火车的班次
%t2=13.00; %直升机允许起飞时间点
%t3=18.00;%直升机结束起飞时间点
ex_cost=[0,1,1.5;
    1,0,2;
    1.5,2,0];%节点处的转换成本
% tr_cost=[0.5,0.01,5];%各种交通工具运输每百公里的成本
tr_cost=[0.1,0.01,1];%各种交通工具运输每百公里的成本
[M,N] = size(PopulationPaths1);% 行数为种群个数，列数为每个个体中包含的方案
C = zeros(pop,N); % 总成本
T = zeros(pop,N); % 总时间
%% 函数主体
for i = 1:M % 循环每一个个体
    for j = 1:N % 循环个体中的每一个元素
        chromosome1 = PopulationPaths1{i,j};
        chromosome2 = PopulationTrans1{i,j};
        for k=1:length(chromosome1)-1
            if chromosome1(k)~=0 && chromosome1(k+1)~=0
                a=chromosome2(k);dis=Dis(chromosome1(k),chromosome1(k+1));
                % DisDisp(k) = dis;
                % 这里好像可以简化，直接用a索引
                C(i,j)=C(i,j)+dis/100*tr_cost(a);
                T(i,j)=T(i,j)+dis/v(a);
            else
            C(i,j)=C(i,j);
            T(i,j)=T(i,j);
            end
        end
        for p = 1:length(chromosome2)-1
            a=chromosome2(p);b=chromosome2(p+1);
            if a~=0 && b~=0 && a~=b
                C(i,j) = C(i,j) + ex_cost(a,b); % 成本
                T(i,j) = T(i,j) + ex_time(a,b); % 时间
            else
                C(i,j)=C(i,j);
                T(i,j) = T(i,j);
            end
        end
    end
end
end