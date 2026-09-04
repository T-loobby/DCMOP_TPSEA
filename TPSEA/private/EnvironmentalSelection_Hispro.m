function Population = EnvironmentalSelection_Hispro(Population,preobj,period,N,VAR)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
if nargin==4
    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,length(Population));
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    %% For Historical interact（本质就是一个个体索引列表）
    idx=[];
    for insno=1:MaxFNo
        cc=find(FrontNo==insno);
        nowCrowdDis=CrowdDis(cc);
        [~,nownumb]=sort(nowCrowdDis,'descend');
        idx=[idx,cc(nownumb)];
    end
elseif nargin==5
    %约束处理分支：实现了改进的ε约束方法
    [~,~,idx] = Auxiliray_task_EnvironmentalSelection2(Population,length(Population),VAR);
end
    %基于历史Perato前沿对当前种群排序
    idx=idx';[~,idx2] = prehelpsort(Population,preobj,length(Population));
    indexforpop=[];indexforpop2=[];
    for zad=1:length(idx)
        indexforpop(idx(zad))=zad;
        indexforpop2(idx2(zad))=zad;
    end
    next1=zeros(1,length(Population));  
    next2=zeros(1,length(Population));
    next1(indexforpop<=N)=10;   %标记约束排名前N的个体
    next2(indexforpop2<=N)=1;   %标记历史排名前N的个体
    newnext=next1+next2;
    maybeworse=find(newnext==10);   %找出可能被淘汰的个体
%     period=0;
    %动态调整淘汰比例（随着period的增加而减少）
    badlength=floor(length(maybeworse)*(1-period));
    %淘汰约束排序中排名靠后的部分
    [~,zxc]=sort(indexforpop(maybeworse),'descend');
    next1(maybeworse(zxc(1:badlength)))=0;
    %找出可能被误解的个体（仅仅是被历史解选中）
    maybebetter=find(newnext==1);
    %用约束排序中的优秀个体替换
    [~,zxc2]=sort(indexforpop2(maybebetter));
    next1(maybebetter(zxc2(1:badlength)))=10;
    %% Population for next generation
    Population = Population(next1==10);
end