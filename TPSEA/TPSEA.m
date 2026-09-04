function Population = TPSEA(Mode, Problem, Population, boundary, N, maxgen, ft, nt, PF, preEvolution, D)

    %% 算法参数
    tau_t = ft;  % 环境变化频率
    n_t = nt;    % 环境变化幅度
    Lower = boundary.lower;
    Upper = boundary.upper;
    Dvalue=[]; IniDvalue= [];
    %% 存档和历史记录
    AllPop = [];      % 历史种群存档
    AllPopU = [];     % 无约束历史种群
    prepop = {};      % 预测种群
    
    %% 环境计数器
    T = 1;            % 环境种类计数器
    gen = 0;          % 总代数计数器
    cc = 0;           % 环境变化后的代数计数器
    VAR0 = 1;         % 初始约束违反值
    UPopulation=Population;
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N);
    %% 迁移参数
    migration_interval = 5;  
    max_migration_size = 10;  % 最大迁移数量
    %% 主循环
    while gen < maxgen
        %% 种群进化
        cc = cc + 1;
        if T >= 3
            % 双种群进化模式
            period = cc / tau_t;
            cp = 1;
            
            % 计算当前种群可行性比例
            feasible_ratio = sum(all([Population1.cons; Population2.cons] <= 0, 2)) / (N);

            % 调整VAR：可行性高时快速收紧，可行性低时保持宽松
            if feasible_ratio > 0.7
                 % 可行性高，快速衰减
                VAR = VAR0 * (1-period)^4;
            elseif feasible_ratio < 0.3
                 % 可行性低，缓慢衰减
                VAR = VAR0 * (1-period)^0.25;
            else
                 % 中等可行性，正常衰减
                VAR = VAR0 * (1-period)^cp*0.3;
            end

            %主种群1的进化
            [~,FrontNo1,CrowdDis1] = EnvironmentalSelection(Population1,1/2*N);
            MatingPool1 = TournamentSelection(2,1/2*N,FrontNo1,-CrowdDis1);
            if period <= 0.2
                Offdec1  = TPSEA_HybridOperator(Population1(MatingPool1),boundary,true);
            else
                Offdec1  = GA(Population1(MatingPool1).decs,boundary);
            end
            % % === 新增：预测生成子代 ===
            % if cc > 0 && ~isempty(V) && any(V)  % 仅在环境变化后启用且V不为0
            %     p_success =max(0, 0.1 - period);  % 线性衰减
            %     [Offdec1, ~] = generateX_V(Offdec1, Population1, p_success, size(Offdec1,1), V, Rn_STD);
            % end
            [Offspring1,~] = Individual(Problem,ft,nt,gen,[],preEvolution,Offdec1);

            %辅助种群的进化
            [~,FrontNo2,CrowdDis2] = EnvironmentalSelectionnocons(Population2,1/2*N);
            MatingPool2 = TournamentSelection(2,1/2*N,FrontNo2,-CrowdDis2);
            if period <= 0.2
                Offdec2  = TPSEA_HybridOperator(Population2(MatingPool2),boundary,false);
            else
                Offdec2  = GA(Population2(MatingPool2).decs,boundary);
            end
            % === 新增：预测生成子代 ===
            % if cc > 0 && ~isempty(V) && any(V)  % 仅在环境变化后启用且V不为0
            %     p_success =max(0, 0.5 - period);  % 线性衰减
            %     [Offdec2, ~] = generateX_V(Offdec2, Population2, p_success, size(Offdec2,1), V, Rn_STD);
            % end
            [Offspring2,~] = Individual(Problem,ft,nt,gen,[],preEvolution,Offdec2);

            %精英保留策略
            [Population1] = EnvironmentalSelection([Population1,Offspring1,Offspring2],1/2*N);
                    
            % 动态融合双排名策略
            if T>=3
                Population2 = EnvironmentalSelection_Hispro([Population2,Offspring2,Offspring1],AllPop((CPF_similar_env-1)*100+1:(CPF_similar_env)*100).objs,period,1/2*N, VAR);
            else
                Population2 = Auxiliray_task_EnvironmentalSelection2([Population2,Offspring2,Offspring1],1/2*N, VAR);
            end
            
            % 合并种群
            [Population(T,:),FrontNo,CrowdDis] = EnvironmentalSelection([Population(T,:),Offspring1,Offspring2],N);
            UPopulation = EnvironmentalSelectionnocons([UPopulation,Population(T,:),Offspring1,Offspring2],N);
        else
            MatingPool = TournamentSelection(2,N,CrowdDis,FrontNo);
            period = cc / tau_t;
            if period <= 0.2
                Offdec  = TPSEA_HybridOperator(Population(T,MatingPool),boundary,true);
            else
                Offdec  = GA(Population(T,MatingPool).decs,boundary);
            end
            [Offspring,~] = Individual(Problem,ft,nt,gen,[],preEvolution,Offdec);
            [Population(T,:),FrontNo,CrowdDis] = EnvironmentalSelection([Population(T,:),Offspring],N);
            UPopulation = EnvironmentalSelectionnocons([UPopulation,Population(T,:),Offspring],N);
        end
        %% 绘图
        % display population and decision vector
        if Mode == 1
            DrawPop(Population(all(Population.cons<=0,2)).objs,Problem,ft,nt);
%             DrawDec(Population(number,:).decs,Problem);
            if size(PF,2) == 1
                Pf = repmat(PF,1,number);
            else
                Pf = PF;
            end
%             Dvalue = DrawMetric(gen,Population(number,:).objs,Pf(number).PF,Dvalue);
%             Dvalue = DrawMetric(gen,Population(all(Population.cons<=0,2)).objs,Pf(number).PF,Dvalue);
        end
        gen = gen + 1;
        %% 环境变化检测
        change = DCchange(Population(T,:),Problem,ft,nt,gen,preEvolution);
        if change > 1e-6
            % 保存当前种群
            SortedPopulation = SortedPopobj(Population(T,:));
            SortedUPopulation = SortedPopobj(UPopulation);
            AllPop = [AllPop, SortedPopulation];
            AllPopU = [AllPopU, SortedUPopulation];
            
            % 反应环境变化
            cc = 0;
            T = T + 1;
            
            if T >= 3
                % 使用RMTT预测新环境种群
                [Dec, UDec, predictedPF,V] = peizhun( Lower, Upper, AllPop, AllPopU, N);
                prepop{T} = predictedPF;
                
                % %=== 新增：计算预测方向V和扰动Rn_STD ===
                % if V==1
                %     C_new = mean(Dec, 1);
                %     C_old = mean(AllPop(end-N+1:end).decs, 1);
                %     V = (C_new - C_old) / norm(C_new - C_old + eps);
                %     Rn_STD = norm(C_new - C_old) / (2*sqrt(D));                        % ======================================
                % end

                % 初始化两个子种群
                randsize = randperm(N);
                randa = randsize(1:0.5*N);  % 保持多样性
                
                Population(T,:) = Individual(Problem,ft,nt,gen,[],preEvolution,Dec);
                Population1 = EnvironmentalSelection(Population(T,:),1/2*N);
                Population2 = Individual(Problem,ft,nt,gen,[],preEvolution,UDec(randa,:));
                Population2 = EnvironmentalSelectionnocons([Population(T,:),Population2],1/2*N);
                UPopulation=[Population(T,:),Population2];
                
                % 计算约束违反值
                cons = [Population.cons; Population2.cons];
                cons(cons<0) = 0;
                VAR0 = max(sum(cons,2));
                if VAR0 == 0, VAR0 = 1; end
                
                % 寻找最相似的历史环境
                Similarity_CPF = 10000*ones(1,T-1);
                for env_idx = max([3,T-10]):T-1
                    if T >= 3
                        Similarity_CPF(env_idx) = sum(min(pdist2(prepop{end}, prepop{env_idx}),[],2));
                    else
                        Similarity_CPF(env_idx) = 10000;
                        prepop{env_idx} = AllPop((env_idx-1)*N+1:(env_idx)*N).objs;
                    end
                end
                [~, CPF_similar_env] = min(Similarity_CPF);
            else
                % 前两个环境随机初始化
                IniDec = unifrnd(repmat(boundary.lower, N, 1), ...
                                repmat(boundary.upper, N, 1));
                Population(T,:) = Individual(Problem,ft,nt,gen,[],preEvolution,IniDec);
            end
        end
        %% preEvaluation
        if gen-preEvolution+1 > 0
            if ~mod(gen-preEvolution,ft)
                T = floor((gen-preEvolution)/ft)+2;
                Population(T,:) = Population(T,:);
                if gen < maxgen
                    if size(PF,2) == 1
                        Pf = repmat(PF,1,T);
                    else
                        Pf = PF;
                    end
%                     IniDvalue = DrawIniMetric(T-1,newPopulation.objs,Pf(T).PF,IniDvalue,Problem);
                end
            end
        end
    end
end
