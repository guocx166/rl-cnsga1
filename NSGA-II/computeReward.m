function reward = computeReward(phi_current, phi_initial)
    % º∆À„Ω±…Õ
if isnan(phi_current) || isnan(phi_initial)  
    reward = -1;
else
    ratio = phi_current / phi_initial;
    if ratio == 1
        reward = 0;
    elseif ratio < 1
        reward = 0.5;
    elseif ratio > 1
        reward = -1;
    end
end
end

