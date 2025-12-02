function [min,path]=dijkstra(w,start,terminal)
%输入变量w为所求图的带权邻接矩阵，start、terminal分别为路径的起点和终点的编号。
%返回path为从start到termial的最短路径以及长度min
%先判断是否是直达终点的路径，如果不是，就把到终点的路径去挑，防止去到终点后又折返
if terminal~=25
    for i=1:25
        w(i,25)=0;
        for j=1:25
            w(25,j)=0;
            if w(i,j)~=1
                w(i,j)=inf;
            end
        end
    end
else
    for i=1:25
        for j=1:25
            if w(i,j)~=1
                w(i,j)=inf;
            end
        end
    end
end
n=size(w,1); label(start)=0; f(start)=start;
%n为所求图的顶点个数，label存放到各点的最短路径，f(v)表示v的父顶点用来还原路径

%初始化将除了start以外的顶点label均设置为无穷大
for i=1:n
    if i~=start
       label(i)=inf;
    end
end

%s数组存放已经搜好的顶点集，初始化只有start
s(1)=start; u=start;
while length(s)<n
    %遍历一遍顶点，将不在顶点集中的顶点选出来进行下面的if判定
    for i=1:n
        ins=0;
        for j=1:length(s)
            if i==s(j)
                ins=1;
            end  
        end
        %判断是否有中继顶点使得它们之间的距离更短，如果有的话更新距离并更新前驱结点
        if ins==0
            v=i;
            if label(v)>(label(u)+w(u,v))
                label(v)=(label(u)+w(u,v)); f(v)=u;
            end 
        end
    end   
    v1=0;
    k=inf;
    %同上再次进行遍历，找到目前最短的路径顶点v1，放入顶点集并改变u的值
    for i=1:n
        ins=0;
        for j=1:length(s)
            if i==s(j)
                ins=1;
            end
        end
        if ins==0
            v=i;
            if k>label(v)
                k=label(v);  v1=v;
            end  
        end
    end
    s(length(s)+1)=v1;  
    u=v1;
end

min=label(terminal); path(1)=terminal;
i=1; 

%按倒序结果推出最短路径
while path(i)~=start
    path(i+1)=f(path(i));
    i=i+1 ;
end
path(i)=start;
L=length(path);
%翻转得到最短路径
path=path(L:-1:1);


