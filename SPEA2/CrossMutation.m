function[child_path,child_tran]=CrossMutation(GAParameters,parent_chromosome1,parent_chromosome2,trans,trans1,trans2,trans3,trans_all)
prob1 = GAParameters.Pc; % 交叉
prob2 = GAParameters.Pm; % 变异
[~,m] = size(parent_chromosome1);
[~,n] = size(parent_chromosome2);

% parent_chromosome1=parent_chromosome(:,1:m);%父代路径染色体
% parent_chromosome2=parent_chromosome(:,m+1:m+n);%父代交通方式染色体
[N,~]=size(parent_chromosome1);
for i=1:2:N
    random1 = rand();%随机产生一个随机数
%     %随机选择两个父代个体，此时的parent_1和parent_2代表选的位置，对于路径和交通工具是公用的。
%     parent_1=round(N*rand(1));
%         if parent_1<1
%             parent_1=1;
%         end
%     parent_2=round(N*rand(1));
%     if parent_2<1
%         parent_2=1;
%     end
%     while parent_1 == parent_2
%         parent_2 = round(N*rand(1));
%         if parent_2 < 1
%             parent_2 = 1;
%         end
%     end
    random_number = randperm(N, 2);
    parent_1 = random_number(1);
    parent_2 = random_number(2);
    parent_path_1 = parent_chromosome1(parent_1,:);%第一条路径
    parent_tran_1 = parent_chromosome2(parent_1,:);%第一条交通
    parent_path_2 = parent_chromosome1(parent_2,:);%第二条路径
    parent_tran_2 = parent_chromosome2(parent_2,:);%第二条交通
    if prob1 >= random1
%%%%%%至此选好了父代的染色体
        l1=0;%计算路径非0部分的长度
        for i1=1:m
            if parent_path_1(i1)~=0
                l1=l1+1;
            end
        end
        l2=0;
        for i2=1:m
            if parent_path_2(i2)~=0
                l2=l2+1;
            end
        end
        l=min(l1,l2);%路径不一样长的时候确保交叉点在短染色体上
        startIndex = randi(l-2)+1;%随机产生一个2到25前的一个数
        %采用单点交叉的模式进行交叉
        child_path_1(1:startIndex-1) = parent_path_1(1:startIndex-1);
        child_path_1(startIndex:m) = parent_path_2(startIndex:m);
        child_path_2(1:startIndex-1) = parent_path_2(1:startIndex-1);
        child_path_2(startIndex:m) = parent_path_1(startIndex:m);
        child_tran_1(1:startIndex-1) = parent_tran_1(1:startIndex-1);%交通工具的交叉节点在路径交叉节点前1位
        child_tran_1(startIndex-1:n) = parent_tran_2(startIndex-1:n);
        child_tran_2(1:startIndex-1) = parent_tran_2(1:startIndex-1);
        child_tran_2(startIndex-1:n) = parent_tran_1(startIndex-1:n);
        child_path_merge=[child_path_1;child_path_2];%合起来进行循环
        child_tran_merge=[child_tran_1;child_tran_2];%合起来进行循环
        for j =1:2%对新产生的自带染色体都进行循环
            a=child_path_merge(j,startIndex-1);b=child_path_merge(j,startIndex);c=child_tran_merge(j,startIndex-1);
            if trans(a,b)~=0%交叉后的路径可以连通       
                    child_path(i+j-1,1:m)=child_path_merge(j,1:m);
                    tr=eval(['trans',num2str(c)]);%索引变量，不知道对不对
                    if tr(a,b)==1%交叉后的交通工具可以连通两点
                        child_tran(i+j-1,1:n)=child_tran_merge(j,1:n);
                    else%交叉后的交通工具不能连通两点
                        child_tran_merge(j,startIndex-1)=get_trans(a,b,trans_all);
                        child_tran(i+j-1,1:n)=child_tran_merge(j,1:n);
                    end
            elseif trans(a,b)==0 %如果交叉后的两点之间不可以连通
                if a~=b%交换路径后两个节点不重复
                    start=a;terminal=b;
                    [~,path1]=dijkstra(trans,start,terminal);
                    l3=length(path1);
                    child_path(i+j-1,1:startIndex-1)=child_path_merge(j,1:startIndex-1);
                    child_path(i+j-1,startIndex-1:startIndex+l3-2)=path1(1:l3);%有重复覆盖一下，防止索引出现0的情况
                    child_path(i+j-1,startIndex+l3-1:l3+m-2)=child_path_merge(j,startIndex+1:m);
                    child_tran(i+j-1,1:startIndex-1)=child_tran_merge(j,1:startIndex-1);
                    for i3=1:l3-1
                        child_tran(i+j-1,startIndex-1+i3-1)=get_trans(path1(i3),path1(i3+1),trans_all);%找出新增路径上的交通工具
                    end
                    child_tran(i+j-1,startIndex+l3-2:n+l3-2)=child_tran_merge(j,startIndex:n);
                else%当交换的路径中出现重复的情况时
                    child_path(i+j-1,1:startIndex-1)=child_path_merge(j,1:startIndex-1);
                    child_path(i+j-1,startIndex:m-1)=child_path_merge(j,startIndex+1:m);
                    child_path(i+j-1,m)=0;
                    child_tran(i+j-1,1:startIndex-1)=child_tran_merge(j,1:startIndex-1);
                    child_tran(i+j-1,startIndex-1:n-1)=child_tran_merge(j,startIndex:n);%重复赋值一个数防止索引出现0的情况
                    child_tran(i+j-1,n)=0;
                end
            end       
        end
        %prob2=0.2;
    else
        child_path(i,1:m)=parent_path_1;child_tran(i,1:n)=parent_tran_1;
        child_path(i+1,1:m)=parent_path_2;child_tran(i+1,1:n)=parent_tran_2;
        %prob2=0.8;%提高变异概率
    end
    random2 = rand();
    if prob2 >= random2%变异
        lm=length(child_path(i,:));
        l4=0;%记录非0路径的长度
        for i4=1:lm%计算非0路径长度
            if child_path(i,i4)~=0
                l4=l4+1;
            end
        end
        mutateIndex = randi(l4-2)+1;%在除去首末的其余点位寻找变异点
        child_path_muta=child_path(i,:);
        child_tran_muta=child_tran(i,:);
        e=child_path(i,mutateIndex-1);f=child_path(i,mutateIndex);g=child_path(i,mutateIndex+1);%找到变异点及其前后的路径节点
        trans_muta=trans;
        trans_muta(e,f)=0;%去除可行路径，保证变异后跟原路径不同，只打断一条可行路径
        [~,path2]=dijkstra(trans_muta,e,g);
        l5=length(path2);
        child_path(i,1:mutateIndex-1)=child_path_muta(1:mutateIndex-1);
        child_path(i,mutateIndex:mutateIndex+l5-3)=path2(2:l5-1);
        child_path(i,mutateIndex+l5-2:lm+l5-3)=child_path_muta(mutateIndex+1:lm);
        child_tran(i,1:mutateIndex-1)=child_tran_muta(1:mutateIndex-1);
        for i5=1:l5-1
            child_tran(i,mutateIndex-1+i5-1)=get_trans(path2(i5),path2(i5+1),trans_all);
        end
        child_tran(i,mutateIndex-1+l5-1-1+1:l4+l5-2-1-1)=child_tran_muta(mutateIndex+1:l4-1);
    end
    %[~,m]=size(child_path);
    %[~,n]=size(child_tran);
    %判断是否存在非法个体，并进行筛选
    t=0;
    for k=1:2
        ll=length(child_path(i+t,:));
        l6=0;%记录非0路径的长度
        for i6=1:ll
           if child_path(i+t,i6)~=0
               l6=l6+1;
           end
        end 
        ln=length(child_tran(i+t,:));
        l7=0;%记录非0交通工具的长度
        for i7=1:ln  
            if child_tran(i+t,i7)~=0
                l7=l7+1;
            end
        end
        if l7+1~=l6%如果存在非法个体，则重新生成路径
            transport=change_transport(child_path(i+t,:),ll,trans_all);
            for i8=1:l7
                child_tran(i+t,i8)=transport(i8);
            end
            child_tran(i+t,l7+1:ln)
        end 
        t=t+1;
    end
end
end
