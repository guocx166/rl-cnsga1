function [dis] = distance(point)
%distance 这里计算各点之间的距离
%   距离为欧式距离，计算结果为一个对称矩阵
for i =1:max(size(point))
    for j=1:max(size(point))
         rij = sqrt( (point(i,1)-point(j,1))^2 + (point(i,2)-point(j,2))^2 );
                tij = round(rij);
                dis(i,j)= tij;
    end
end

