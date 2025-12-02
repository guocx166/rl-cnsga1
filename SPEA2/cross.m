function[childPath,childTran] = cross(N,prob1,prob2,parent_chromosome1,parent_chromosome2,trans,transport1,transport2,transport3,trans_all,Dis)
% 染色体的交叉,生成的为元胞形式
% N为种群大小，parent_chromosome1和2都是N*1列的元胞
% 使用函数
childPath = cell(N,1);
childTran = cell(N,1);
max_iterations = 10; % 设置最大循环上限
for i=1:2:N
    random1 = rand();%随机产生一个随机数
    %随机选择两个父代个体，此时的parent_1和parent_2代表选的位置，对于路径和交通工具是公用的。
    parent_1=round(N*rand(1));
        if parent_1<1
            parent_1=1;
        end
    parent_2=round(N*rand(1));
    if parent_2<1
        parent_2=1;
    end
    iteration = 0;
    while isequal(parent_chromosome1(parent_1),parent_chromosome1(parent_2)) && iteration < max_iterations
        parent_2 = round(N*rand(1));
        if parent_2 < 1
            parent_2 = 1;
        end
        iteration = iteration + 1;
    end
%     if isequal(parent_chromosome1(parent_1), parent_chromosome1(parent_2))
%     % 根据需求给parent_chromosome1(parent_2)赋新的值
%     parent_chromosome1(parent_2) = 
%     [parent_chromosome1(parent_2),parent_chromosome2] = generatePaths_dis(trans,trans_all, A(1), B(j), 1); % 请根据情况设定newValue
%     end
    % 转化为矩阵的格式
    parent_path_1 = cell2mat(parent_chromosome1(parent_1));%第一条路径,转化为矩阵的格式
    parent_tran_1 = cell2mat(parent_chromosome2(parent_1));%第一条交通
    parent_path_2 = cell2mat(parent_chromosome1(parent_2));%第二条路径
    parent_tran_2 = cell2mat(parent_chromosome2(parent_2));%第二条交通
%%%%%% 至此选好了父代的染色体 
    %% 交叉  
    lp1 = nnz(parent_path_1);
    lp2 = nnz(parent_path_2);
    lt1 = nnz(parent_tran_1);
    lt2 = nnz(parent_tran_2);
    l=min(lp1,lp2);%路径不一样长的时候确保交叉点在短染色体上
    if prob1 >= random1 && l>2 % 如果可以直接到达，则不进行交叉
        % 计算非0长度
        startIndex = randi(l-2)+1;%随机产生一个染色体长度内的数
        %采用单点交叉的模式进行交叉,交叉后两个染色体的长度交换了
        child_path_1(1:startIndex-1) = parent_path_1(1:startIndex-1);
        child_path_1(startIndex:lp2) = parent_path_2(startIndex:lp2);
        child_path_2(1:startIndex-1) = parent_path_2(1:startIndex-1);
        child_path_2(startIndex:lp1) = parent_path_1(startIndex:lp1);
        child_tran_1(1:startIndex-1) = parent_tran_1(1:startIndex-1);%交通工具的交叉节点在路径交叉节点前1位
        child_tran_1(startIndex:lt2) = parent_tran_2(startIndex:lt2);
        child_tran_2(1:startIndex-1) = parent_tran_2(1:startIndex-1);
        child_tran_2(startIndex:lt1) = parent_tran_1(startIndex:lt1);
        %% 检查路径是否连通并且进行变异
        [childPath{i},childTran{i}] = check_mutation(prob2,child_path_1,child_tran_1,startIndex,trans,transport1,transport2,transport3,trans_all,Dis);
        [childPath{i+1},childTran{i+1}] = check_mutation(prob2,child_path_2,child_tran_2,startIndex,trans,transport1,transport2,transport3,trans_all,Dis);
    else 
        childPath{i}=parent_path_1(parent_path_1~=0);
        childTran{i}=parent_tran_1(parent_tran_1~=0);
        childPath{i+1}=parent_path_2(parent_path_2~=0);
        childTran{i+1}=parent_tran_2(parent_tran_2~=0);
    end
end