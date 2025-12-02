function[chromosome]=mutation(filename,m,n,trans,trans1,trans2,trans3,trans_all)
%灾变操作
%部分初始化种群，并重新开始进化，当拥挤度为0的个体超过种群大小的1/2时，开始进行灾变。
transport=[];
f1=xlsread(filename);
[N1,M1]=size(f1);
for i=1:N1
    l=0;
    for j=1:M1
        if f1(i,j)~=0
            l=l+1;
        end
    end
    transport(i,1:l-1)=change_transport(f1(i,:),l,trans_all);
end
[N2,M2]=size(transport);
f1(:,M1+1:m)=0;transport(:,M2+1:n)=0;
chromosome=[f1,transport];
end
