function [Population,idx] = prehelpsort(Population,hisccpopulation,N)

    hiscpopulation=hisccpopulation;
    idx=[];tt=1;
        %对历史前沿hisccpopulation进行非支配排序，由于属于同一前沿，故为1
        [FrontNo,MaxFNo] = NDSort(hisccpopulation,1);
        %计算历史前沿中解的拥挤度
        CrowdDis = CrowdingDistance(hisccpopulation,FrontNo);
        indidx=[];
        %仅仅处理第一前沿的解，分布均匀的排序在前
        for insno=1:1
            cc=find(FrontNo==insno);
            nowCrowdDis=CrowdDis(cc);
            [~,nownumb]=sort(nowCrowdDis,'descend');
            indidx=[indidx,cc(nownumb)];
        end 
        ndhispopulation=hisccpopulation(indidx,:);

        [~,sortdis]=sort(pdist2(ndhispopulation,Population.objs),2);
        %构建选择索引
        while length(idx)<N
            idx=[idx,sortdis(:,tt)'];
            idx=unique(idx,'stable');
            tt=tt+1;
        end
        %返回结果
        Population=Population(idx(1:N));
end