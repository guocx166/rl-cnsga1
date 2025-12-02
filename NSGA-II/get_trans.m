function z = get_trans(x,y,trans_all)
%UNTITLED3 用于判断交通方式，输出为交通方式序号，即1,2或3
%   x,y为两个节点，z为节点间交通工具，

if trans_all(x,y)==1
    z=1;
elseif trans_all(x,y)==2
    z=2;
elseif trans_all(x,y)==3
    z=3;
elseif trans_all(x,y)==4
    z=unidrnd(2);%随机选择1、2
elseif trans_all(x,y)==5
    A=[1,3];
    z=A(randi(numel(A),1,1));
    z=sort(z);%随机选择1、3
elseif trans_all(x,y)==6
    A=[2,3];
    z=A(randi(numel(A),1,1));
    z=sort(z);%随机选择2、3
elseif trans_all(x,y)==7
    z=unidrnd(3);%随机选择1、2、3
end
end



