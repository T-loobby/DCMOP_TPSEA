function [curve_B_pred, tps_model] = predict_tps_deformation(A, B, lambda, eval_density)
    % 专为有序曲线设计的TPS预测，确保预测曲线的空间跨度不小于原始曲线
    % 输入:
    %   A - 原始曲线点(100×2)，按顺序连接构成曲线
    %   B - 目标曲线点(100×2)，与A点对点对应
    %   lambda - 正则化参数
    %   eval_density - 每两个原始点之间的插值点数
    
    % 参数验证
    if nargin < 4
        eval_density = 5; % 默认每段插值5个点
    end

    % 计算原始曲线的空间跨度（最远点距离）
    original_span = max(pdist(A));

    % 计算TPS系数（保持有序性）
    [w, v] = find_tps_coefficients_ordered(A, B, lambda);
    
    % 创建密集评估点（保持曲线顺序）
    t = linspace(0, 1, size(A,1))';
    t_eval = linspace(0, 1, eval_density*size(A,1))';
    
    % 参数化插值（保持顺序）
    A_eval = [interp1(t, A(:,1), t_eval, 'spline'), ...
              interp1(t, A(:,2), t_eval, 'spline')];
    
    % 预测形变（保持顺序）
    [fX, fY] = deform_curve_tps_ordered(A_eval, A, w, v);
    curve_B_pred = [fX, fY];
    
    % 确保预测曲线的空间跨度不小于原始曲线
    current_span = max(pdist(curve_B_pred));
    if current_span < original_span
        % 计算缩放因子
        scale_factor = original_span / current_span;
        % 计算曲线中心
        center = mean(curve_B_pred);
        % 缩放曲线
        curve_B_pred = (curve_B_pred - center) * scale_factor + center;
    end
    
    % 保存模型参数
    tps_model.control_points = A;
    tps_model.mapping_coeffs = w;
    tps_model.poly_coeffs = v;
    tps_model.lambda = lambda;
    
    % 顺序可视化
    %visualize_ordered_match(A, B, curve_B_pred);
end

%% 核心改进子函数
function [w, v] = find_tps_coefficients_ordered(ctrl_pts, target_pts, lambda)
    % 考虑曲线顺序的TPS系数计算
    n = size(ctrl_pts,1);
    
    % 计算径向基函数矩阵（添加顺序约束）
    K = zeros(n);
    for i = 1:n
        for j = 1:n
            r = norm(ctrl_pts(i,:) - ctrl_pts(j,:));
            if r > 0
                K(i,j) = r^2 * log(r);
            end
        end
    end
    
    % 添加邻近点约束权重
    for i = 2:n-1
        K(i,i) = K(i,i) + lambda * (1 + norm(ctrl_pts(i,:)-ctrl_pts(i-1,:)) + norm(ctrl_pts(i,:)-ctrl_pts(i+1,:)));
    end
    
    P = [ones(n,1), ctrl_pts];
    M = [K, P; P', zeros(3)];
    Y = [target_pts; zeros(3,2)];
    
    coeffs = pinv(M) * Y; % 使用伪逆提高稳定性
    w = coeffs(1:n,:);
    v = coeffs(n+1:end,:);
end

function [fX, fY] = deform_curve_tps_ordered(points, ctrl_pts, w, v)
    % 保持顺序的形变计算
    n_ctrl = size(ctrl_pts,1);
    n_points = size(points,1);
    U = zeros(n_points, n_ctrl);
    
    for i = 1:n_points
        for j = 1:n_ctrl
            r = norm(points(i,:) - ctrl_pts(j,:));
            if r > 0
                U(i,j) = r^2 * log(r);
            end
        end
    end
    
    V = [ones(n_points,1), points];
    deformed = U * w + V * v;
    fX = deformed(:,1);
    fY = deformed(:,2);
end

%% 顺序可视化
function visualize_ordered_match(A, B, B_pred)
    figure('Position', [100 100 1200 500]);
    
    % 子图1: 顺序对应关系
    subplot(1,2,1);
    plot(A(:,1), A(:,2), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'b'); 
    hold on;
    plot(B(:,1), B(:,2), 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    title('有序对应关系', 'FontSize', 12);
    legend('曲线A（有序）', '曲线B（有序）', 'Location', 'best');
    axis equal tight;
    grid on;
    
    % 子图2: 预测结果（带顺序）
    subplot(1,2,2);
    plot(A(:,1), A(:,2), 'b-', 'LineWidth', 2); 
    hold on;
    plot(B(:,1), B(:,2), 'r-', 'LineWidth', 2);
    plot(B_pred(:,1), B_pred(:,2), 'g-', 'LineWidth', 2);
    title('有序TPS预测（确保空间跨度）', 'FontSize', 12);
    legend('原始曲线A', '目标曲线B', '预测曲线', 'Location', 'best');
    axis equal tight;
    grid on;
    set(gcf, 'Color', 'w');
end