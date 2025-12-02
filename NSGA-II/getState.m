function state = getState(phi_10, phi_20, phi_1, phi_2)
    % 根据种群多样性比率获得状态
    % 使用函数
    condition1 = (phi_1 > phi_10) && (phi_2 > phi_20);
    condition2 = (phi_1 > phi_10) && (phi_2 == phi_20);
    condition3 = (phi_1 > phi_10) && (phi_2 < phi_20);
    condition4 = (phi_1 == phi_10) && (phi_2 > phi_20);
    condition5 = (phi_1 == phi_10) && (phi_2 == phi_20);
    condition6 = (phi_1 == phi_10) && (phi_2 < phi_20);
    condition7 = (phi_1 < phi_10) && (phi_2 > phi_20);
    condition8 = (phi_1 < phi_10) && (phi_2 == phi_20);
    condition9 = (phi_1 < phi_10) && (phi_2 < phi_20); 
    % 状态编号生成
% 根据条件生成状态编号
    % 这里用一个简单的if-elseif链来给每个状态分配唯一的编号
    if condition1
        state = 1;
    elseif condition2
        state = 2;
    elseif condition3
        state = 3;
    elseif condition4
        state = 4;
    elseif condition5
        state = 5;
    elseif condition6
        state = 6;
    elseif condition7
        state = 7;
    elseif condition8
        state = 8;
    elseif condition9
        state = 9;
    else
        state = 1; % 如果所有的解的值完全一致的话，那么就认为其多样性都在减少（这个时候可以引入灾变机制）
    end
end

