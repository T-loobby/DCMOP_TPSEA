function [Nobj, Fobj] = obtain_matrix_mean(Pop_1, Pop_2, N)
    % 获取两个种群的目标值
    Nobj = Pop_1.objs; 
    original_Fobj = Pop_2.objs;
    
    %% 划界
    %{
    dotProducts = sum(Nobj(1:end-1, :) .* Nobj(2:end, :), 2);
    norms_i = sqrt(sum(Nobj(1:end-1, :).^2, 2));
    norms_i_plus_1 = sqrt(sum(Nobj(2:end, :).^2, 2));
    cosTheta = dotProducts ./ (norms_i .* norms_i_plus_1);
    cosTheta = max(min(cosTheta, 1), -1);
    theta = acos(cosTheta);
    naber_angle = rad2deg(theta);

    % 设定阈值并筛选符合条件的相邻点对
    threshold = 7;
    is_large_angle = naber_angle > threshold;
    [m,~]=size(Nobj(is_large_angle, :));
    Devided_ps=[];
    Devided_ps(1:2:2*m-1,:)=Nobj(is_large_angle, :);
    Devided_ps(2:2:2*m,:)=Nobj([false; is_large_angle], :);
    Rooms = dividePoints(Nobj, original_Fobj, Devided_ps);
    %}
    %% 判断皱褶，匹配，预测
    dot_products = original_Fobj * Nobj'; 
    norms_Fobj = sqrt(sum(original_Fobj.^2, 2)); 
    norms_Nobj = sqrt(sum(Nobj.^2, 2)); 
    norm_products = norms_Fobj * norms_Nobj'; 
    cos_theta = dot_products ./ norm_products;
    cos_theta(norm_products == 0) = 0; 
    angle_matrix = acos(cos_theta) * (180 / pi);

    % 初始匹配：为每个Pop_1个体找到Pop_2中最相似的个体
    [~, min_angle_indices] = min(angle_matrix, [], 1);
    
    % 找出Pop_2中未被匹配的个体索引
    un_match_F = setdiff(1:N, min_angle_indices);
    
    % 初始化Fobj（匹配后的Pop_2目标值）
    Fobj = original_Fobj(min_angle_indices, :); % 初始匹配部分
    
    % 计算所有匹配点对的平均位移向量
    matched_pairs = [min_angle_indices' (1:N)']; % [Pop_2索引, Pop_1索引]
    displacement_vectors = Nobj(matched_pairs(:,2),:) - original_Fobj(matched_pairs(:,1),:);
    avg_displacement = mean(displacement_vectors, 1);
    
    % 对未匹配的点应用平均位移向量
    for i = 1:length(un_match_F)
        a_idx = un_match_F(i); % Pop_2中未匹配的个体索引
        Fobj(a_idx, :) = original_Fobj(a_idx, :) + avg_displacement;
    end
    
    % 确保Fobj顺序正确（与Nobj对应）
    Fobj = Fobj(1:N, :); % 保持与Nobj相同大小
    
    % Nobj保持不变（即Pop_1.objs）
end