function [SortedPopulation] = SortedPopobj(Population)
   
    PopObj = Population.objs;
    [N, ~] = size(PopObj);

    % 归一化目标值
    minObj = min(PopObj,[],1);
    maxObj = max(PopObj,[],1);
    rangeObj = maxObj - minObj;
    rangeObj(rangeObj==0) = 1; % 处理常值目标
    normObj = (PopObj - repmat(minObj,N,1)) ./ (repmat(rangeObj,N,1));
    
    % 计算两两角度矩阵（上三角部分）
    Angle = acos(1-pdist2(normObj, normObj, 'cosine'));
    Angle = triu(Angle, 1); % 取上三角避免重复和自比较
    
    % 找出最大角度对应的解对
    [~, max_idx] = max(Angle(:));
    [sol1_idx, sol2_idx] = ind2sub([N,N], max_idx);
    
    % 选择x坐标较小的作为种子点
    if PopObj(sol1_idx,1) < PopObj(sol2_idx,1)
        seed_idx = sol1_idx;
    else
        seed_idx = sol2_idx;
    end
    seed_point = normObj(seed_idx, :);
    
    % 计算所有解与种子点的距离
    distances = acos(1 - pdist2(normObj, seed_point, 'cosine'));
    
    % 创建排序索引（种子点排第一，其余按距离排序）
    [~, sorted_idx] = sort(distances);
    
    % 确保种子点在最前面
    if sorted_idx(1) ~= seed_idx
        % 如果种子点不在第一位，找到它并移到第一位
        seed_pos = find(sorted_idx == seed_idx);
        sorted_idx = [seed_idx; sorted_idx(1:seed_pos-1); sorted_idx(seed_pos+1:end)];
    end
    
    % 创建重新排序后的种群
    SortedPopulation = Population(sorted_idx);
end