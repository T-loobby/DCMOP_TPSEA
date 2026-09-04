function [selected_points,V] = select_points(points, num_sectors, obj,V)
    
    % 检查obj是否与points基本重合（obj的点都能在points中找到近邻）
    use_uniform = false;
    if ~isempty(obj) && size(obj,1) >= 2
        % 计算obj中每个点到points的最小距离
        min_dists = min(pdist2(obj, points), [], 2);
        % 如果obj中所有点都能在points中找到足够近的点（阈值0.015）
        if all(min_dists < 0.015)
            use_uniform = true;
            V=0;
        end
    end
    
    % 计算抽取间隔
    M = size(points, 1);
    step = max(floor(M / num_sectors), 1);  % 确保步长至少为1

    % 生成等间隔索引
    idx = 1:step:min(step*num_sectors, M);  % 确保不超出范围

    % 如果抽取的点数不足N，从末尾补充
    if length(idx) < num_sectors
        additional_idx = M - (num_sectors - length(idx)) + 1:M;
        idx = [idx, additional_idx];
    end

    % 抽取对应的点
    selected_points = points(idx, :);
end