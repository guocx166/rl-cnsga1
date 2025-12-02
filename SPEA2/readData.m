function [point] = readData(data)
% 读取文件并画图，将数据放大50倍，选取供应点，中转节点和需求点并通过蒙特卡洛模拟生成最初各个需求点的需求量
% 分三部分，读取数据;扩大距离;指定供应、中转节点
% 读取txt文件 修改后  
% data = dlmread(file_name); % 读取 TXT 文件
% 写入CSV文件
outputFile = 'output_data.csv'; % 输出的 CSV 文件名  
csvwrite(outputFile, data); % 保存数据为 CSV 格式 
% 读取CSV文件  
dataTable = readtable(outputFile);
data = table2array(dataTable(:, 1:3)); % 提取数据的前3列并转换为数组  
% 截取前两列：序号和横坐标（第三列实际上在下面的plot语句中也会被使用）  
[x,~] = size(data);
ids = data(:,1);         % 序号  
X = data(:,2);           % 横坐标  
Y = data(:,3);           % 纵坐标  
city1 = data(:,2:3);
% 扩大距离
for i = 1:x-1 %将欧式距离扩大40倍
    for j =1:x
        X(j)=20*(city1(j,1)-city1(i,1))-city1(i,1);
        
        Y(j)=20*(city1(j,2)-city1(i,2))-city1(i,2);
    end
end
% 绘制散点图  
figure;  
scatter(X, Y, 'filled');  
hold on;  
% 为每个点添加序号标注  
for i = 1:length(ids)  
    text(X(i), Y(i), num2str(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');  
end  
point = [X,Y];
xlabel('X轴');  
ylabel('Y轴');  
title('坐标散点图');  
grid on; % 添加网格  
hold off; % 释放图形  
axis equal; % 轴比例相同
end

