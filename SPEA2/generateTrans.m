function [f2] = generateTrans(trans_all,path)
% 初始化种群中的交通工具，采用实数编码方式
% f2为交通方式
% 使用函数
[x,y]=size(path);
%根据路径去选择对应的交通工具
for m=1:x
    for n =1:length(path(m,:))-1%生成的交通工具染色体比路径少1
        a=path(m,n);b=path(m,n+1);
        if a~=0 && b~=0
            if trans_all(a,b)==1
                transport(m,n)=1;
            elseif trans_all(a,b)==2
                transport(m,n)=2;
            elseif trans_all(a,b)==3
                transport(m,n)=3;
            elseif trans_all(a,b)==4
                transport(m,n)=unidrnd(2);%随机选择1、2
            elseif trans_all(a,b)==5
                A=[1,3];
                transport(m,n)=A(randi(numel(A),1,1));
                transport(m,n)=sort(transport(m,n));%随机选择1、3
            elseif trans_all(a,b)==6
                A=[2,3];
                transport(m,n)=A(randi(numel(A),1,1));
                transport(m,n)=sort(transport(m,n));%随机选择2、3
            elseif trans_all(a,b)==7
                transport(m,n)=unidrnd(3);%随机选择1、2、3
            end
        else
            transport(m,n)=0;
        end
    end
end
f2=transport;
