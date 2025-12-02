function [PopulationPaths1,PopulationPaths2,PopulationTrans1,PopulationTrans2] = generatePop(filename1,filename2,trans_all,trans,Dis,L,pop)
%UNTITLED 此处显示有关此函数的摘要
%   根据文件输入函数
% 使用函数
dataTable1 = readtable(filename1);
dataTable1 = table2array(dataTable1);
dataTable2 = readtable(filename2);
dataTable2 = table2array(dataTable2);
[md1,~] = size(dataTable1);
% start = 1;
% terminal = 10;
% C = [SupplyPoint];
% [shortestPath] = dijkstra1(Dis, start, terminal, C);

PopulationPaths1 = cell(L, pop);
PopulationPaths2 = cell(L,pop);
%%
for t = 1:L
    for i = 1:pop
        pathtmp = dataTable1(i+(t-1)*pop, :); % 对于每一行元素进行连通性处理
        pathtmp = pathtmp(pathtmp~=0); % 去掉0元素
        ll = nnz(pathtmp); % 求其长度
        pathnew = [pathtmp(1)];
        pointend = pathtmp(ll);
        for j = 1:ll-1
            a = pathtmp(j); b = pathtmp(j+1);
            if trans(a,b) == 1
                pathnew = [pathnew b];
            elseif trans(a, b) == 0
                pathadd = dijkstra1(Dis, a, b, pointend);
                pathnew = [pathnew pathadd(2:end)];
            end
        end
        transnew = generateTrans(trans_all,pathnew);
        PopulationPaths1{t,i} = pathnew;
        PopulationTrans1{t,i} = transnew;
    end
end

for t = 1:L
    for i = 1:pop
        pathtmp = dataTable2(i+(t-1)*pop, :); % 对于每一行元素进行连通性处理
        pathtmp = pathtmp(pathtmp~=0); % 去掉0元素
        ll = nnz(pathtmp); % 求其长度
        pathnew = [pathtmp(1)];
        pointend = pathtmp(ll);
        for j = 1:ll-1
            a = pathtmp(j); b = pathtmp(j+1);
            if trans(a,b) == 1
                pathnew = [pathnew b];
            elseif trans(a, b) == 0
                pathadd = dijkstra1(Dis, a, b, pointend);
                pathnew = [pathnew pathadd(2:end)];
            end
            transnew = generateTrans(trans_all,pathnew);
        PopulationPaths2{t,i} = pathnew;
        PopulationTrans2{t,i} = transnew;
        end
    end
end
PopulationPaths1 = PopulationPaths1';
PopulationTrans1 = PopulationTrans1';
PopulationPaths2 = PopulationPaths2';
PopulationTrans2 = PopulationTrans2';
end

