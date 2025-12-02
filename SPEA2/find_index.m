function [a] = find_index(b,c)
%UNTITLED2 此处显示有关此函数的摘要
% %  , 找到矩阵中比元素大的数的位置,
% c=[8,10,12,14,16,18,20,22];
% b=14;
for i=1:length(c)
    if b<=c(i)
        a=i;
    break
    end
end

