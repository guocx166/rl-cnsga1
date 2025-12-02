function[child_path_output,child_tran_output] = check_mutation(prob2,child_path_1,child_tran_1,startIndex,trans,transport1,transport2,transport3,trans_all,Dis)
% 检查路径是否连通并进行交叉 输出的是一个数组
% 使用函数
%检查路径是否连通
    m = length(child_path_1);
    n = length(child_tran_1);
    a=child_path_1(startIndex-1);b=child_path_1(startIndex);c=child_tran_1(startIndex-1);
    if trans(a,b)~=0 % 交叉后的路径可以连通       
            % child_path(i+j-1,1:m)=child_path_merge(j,1:m);
        tr=eval(['transport',num2str(c)]);%索引变量，不知道对不对
        if tr(a,b)==1 % 交叉后的交通工具可以连通两点
            % child_path(i+j-1) = child_path_merge(j); % 可以连通直接赋值给元胞
            child_path=child_path_1;
            child_tran=child_tran_1;
        else % 交叉后的交通工具不能连通两点
            child_path = child_path_1;
            child_tran_1(startIndex-1)=get_trans(a,b,trans_all);
            child_tran=child_tran_1;
        end
    elseif trans(a,b)==0 % 如果交叉后的两点之间不可以连通
        if a~=b%交换路径后两个节点不重复
            % 连通路径
            start=a;terminal=b;
            [~,path1]=dijkstra(Dis,start,terminal);
            l3=length(path1);
            child_path(1:startIndex-1)=child_path_1(1:startIndex-1);
            child_path(startIndex-1:startIndex+l3-2)=path1(1:l3);%有重复覆盖一下，防止索引出现0的情况
            child_path(startIndex+l3-1:l3+m-2)=child_path_1(startIndex+1:m);
            % 交通方式
            child_tran(1:startIndex-1)=child_tran_1(1:startIndex-1);
            %child_tran(startIndex-1) = get_trans(a,path1(2),trans_all);
            for i3=1:l3-1
                child_tran(startIndex+i3-2)=get_trans(path1(i3),path1(i3+1),trans_all);%找出新增路径上的交通工具
            end
            child_tran(startIndex+l3-2:n+l3-2)=child_tran_1(startIndex:n);
        else % 当交换的路径中出现重复的情况时
            child_path(1:startIndex-1)=child_path_1(1:startIndex-1);
            child_path(startIndex:m-1)=child_path_1(startIndex+1:m);
            child_path(m)=0;
            child_tran(1:startIndex-1)=child_tran_1(1:startIndex-1);
            child_tran(startIndex-1:n-1)=child_tran_1(startIndex:n);%重复赋值一个数防止索引出现0的情况
            child_tran(n)=0;
        end
    end
    child_path_output = child_path(child_path~=0);
    child_tran_output = child_tran(child_tran~=0);
    %% 变异
    random2 = rand();
    if prob2 >= random2%变异
        l4=nnz(child_path);%记录非0路径的长度
        mutateIndex = randi(l4-2)+1;%在除去首末的其余点位寻找变异点
        e=child_path(mutateIndex-1);f=child_path(mutateIndex);h = child_path(mutateIndex+1);g=child_path(mutateIndex+1);%找到变异点及其前后的路径节点
        trans_muta=Dis;
        trans_muta(e,f) = inf; trans_muta(f,h) = inf; % 去除可行路径，保证变异后跟原路径不同，只打断一条可行路径
        [~,path2]=dijkstra(trans_muta,e,h);
        l5=length(path2);
        if l5 == 2 % 如果可以直接连通
            child_path_muta(1:mutateIndex-1)=child_path(1:mutateIndex-1);
            % child_path_muta(mutateIndex:mutateIndex+l5-3)=path2(2:l5-1);
            child_path_muta(mutateIndex+l5-2:l4+l5-3)=child_path(mutateIndex+1:l4);
            child_tran_muta(1:mutateIndex-2)=child_tran(1:mutateIndex-2);
            for i5=1:l5-1
                child_tran_muta(mutateIndex-2+i5) = get_trans(path2(i5),path2(i5+1),trans_all);
            end
            child_tran_muta(mutateIndex+l5-2:l4-2)=child_tran(mutateIndex+1:l4-1);
            child_path_output = child_path_muta(child_path_muta~=0);
            child_tran_output = child_tran_muta(child_tran_muta~=0);
        else % 不能直接连通需要重新连接时    
            child_path_muta(1:mutateIndex-1)=child_path(1:mutateIndex-1);
            child_path_muta(mutateIndex:mutateIndex+l5-3)=path2(2:l5-1);
            child_path_muta(mutateIndex+l5-2:l4+l5-3)=child_path(mutateIndex+1:l4);
            child_tran_muta(1:mutateIndex-2)=child_tran(1:mutateIndex-2);
            for i5=1:l5-1
                child_tran_muta(mutateIndex-2+i5) = get_trans(path2(i5),path2(i5+1),trans_all);
            end
            child_tran_muta(mutateIndex+l5-2:l4+l5-2-1-1)=child_tran(mutateIndex+1:l4-1);
            child_path_output = child_path_muta(child_path_muta~=0);
            child_tran_output = child_tran_muta(child_tran_muta~=0);
    end
    %% 判断是否存在非法个体，并进行筛选
    l6=nnz(child_path_output);%记录非0路径的长度
    l7=nnz(child_tran_output);%记录非0交通工具的长度
    if l7+1~=l6%如果存在非法个体，则重新生成路径
        child_tran_new = generateTrans(trans_all,child_path_output);
        child_tran_output = child_tran_new(child_tran_new~=0);
    end
end