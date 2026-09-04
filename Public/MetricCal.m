function metric = MetricCal(Population,PF,runtime)
%METRICCAL Calculate mean IGD and HV over dynamic environments.

    environmentCount = size(Population,1);
    if size(PF,2) == 1
        PF = repmat(PF,1,environmentCount);
    end

    scoreIGD = zeros(1,environmentCount);
    scoreHV = zeros(1,environmentCount);
    for i = 1:environmentCount
        current = Population(i,:);
        current(any(current.cons > 0,2)) = [];
        if isempty(current)
            scoreIGD(i) = 10;
            scoreHV(i) = 0;
            continue;
        end

        popObj = current(NDSort(current.objs,1) == 1).objs;
        scoreIGD(i) = IGD(popObj,PF(i).PF);
        scoreHV(i) = HV(popObj,PF(i).PF);
    end

    metric.runtime = runtime;
    metric.IGD = mean(scoreIGD);
    metric.HV = mean(scoreHV);
end
