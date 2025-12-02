function [trans,trans1,trans2,trans3,trans_all] = read_transport(filename1,filename2,filename3)
%得到可以连通的节点，其余输出为各自交通方式连接节点,trans表示节点是否可以连通
%trans_all说明节点处可选用的交通工具，其中1为公路，2为铁路，3为直升机，4为公路、铁路，5为公路、直升机；6为铁路、直升机；7为所有交通工具均可
% 此处显示详细说明
trans1=xlsread(filename1);
trans2=xlsread(filename2);
trans3=xlsread(filename3);
trans1=trans1(2:26,2:26);%公路运输
trans2=trans2(2:26,2:26);%铁路运输
trans3=trans3(2:26,2:26);%飞机运输
[m,n]=size(trans1);
for i=1:m
    for j =1:m
        if trans1(i,j)==1
            trans(i,j)=1;
        elseif trans2(i,j)==1 
            trans(i,j)=1;
        elseif trans3(i,j)==1
            trans(i,j)=1;
        else
            trans(i,j)=0;
        end
    end
end
trans_all=[];
%trans_all=zeros(m,n);
for i=1:m
    for j =1:m
        if trans1(i,j)==1 && trans2(i,j)==1 && trans3(i,j)==1
            trans_all(i,j)=7;
        elseif trans1(i,j)==1 &&trans2(i,j)==1 && trans3(i,j)~=1
            trans_all(i,j)=4;
        elseif trans1(i,j)==1 &&trans2(i,j)~=1 && trans3(i,j)~=1
            trans_all(i,j)=1;
        elseif trans1(i,j)~=1 && trans2(i,j)==1 && trans3(i,j)==1
            trans_all(i,j)=5;
        elseif trans1(i,j)~=1 && trans2(i,j)==1 && trans3(i,j)~=1
            trans_all(i,j)=2;
        elseif trans1(i,j)~=1 && trans2(i,j)~=1 && trans3(i,j)==1
            trans_all(i,j)=3;
        elseif trans1(i,j)==1 && trans2(i,j)~=1 &&trans3(i,j)==1
            trans_all(i,j)=6;
        end
    end
end

