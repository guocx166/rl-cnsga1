function [transport] = change_transport(path,l,trans_all)
%根据路径生成交通工具
%   此处显示详细说明

for n =1:l-1%生成的交通工具染色体比路径少1
    a=path(n);b=path(n+1);
    if a~=0 && b~=0
        if trans_all(a,b)==1
            transport(n)=1;
        elseif trans_all(a,b)==2
            transport(n)=2;
        elseif trans_all(a,b)==3
            transport(n)=3;
        elseif trans_all(a,b)==4
            transport(n)=unidrnd(2);%随机选择1、2
        elseif trans_all(a,b)==5
            A=[1,3];
            transport(n)=A(randi(numel(A),1,1));
            transport(n)=sort(transport(n));%随机选择1、3
        elseif trans_all(a,b)==6
            A=[2,3];
            transport(n)=A(randi(numel(A),1,1));
            transport(n)=sort(transport(n));%随机选择2、3
        elseif trans_all(a,b)==7
            transport(n)=unidrnd(3);%随机选择1、2、3
        end
    else
        transport(n)=0;
    end
end
end

