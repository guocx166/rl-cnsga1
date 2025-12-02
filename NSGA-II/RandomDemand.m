function DemandReshape = RandomDemand(pop,n)
%UNTITLED 此处显示有关此函数的摘要
%   生成各个点的需求
num = pop*n;
Beta = 0.4; % 控制预算水平
% Beta = 0.4; % 控制预算水平
Tao = Beta * num; % 不确定需求预算的取值范围
Q1 = 100; % 假设平均需求量为100
Epsilon = 0.5; 
Q2 = Epsilon * Q1; % 需求量的偏差
random_numbers = rand(1, num);
total_sum = sum(random_numbers);
% 如果和大于Tao，进行调整, 控制其预算范围
while total_sum > Tao  
        scalingFactor = Tao / total_sum;  
        random_numbers = random_numbers * scalingFactor;
        total_sum = sum(random_numbers);
end
for i = 1:num
        Eta = randi([0, 1]) * 2 - 1;
        Demand(i) =  Q1 + Q2 * random_numbers(i)*Eta;
end
DemandReshape = reshape(Demand, [pop, n]);
