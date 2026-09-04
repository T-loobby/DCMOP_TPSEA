function OffDec = TPSEA_HybridOperator(Population,boundary,useCons)
% TPSEA_HybridOperator - Generate half offspring by DE and half by GA.

    if nargin < 3
        useCons = true;
    end

    Parent = Population.decs;
    N = size(Parent,1);

    DEDec = TPSEA_DEOperator(Population,boundary,useCons);
    GADec = GA(Parent,boundary);

    if size(GADec,1) < N
        GADec = [GADec; DEDec(size(GADec,1)+1:N,:)];
    elseif size(GADec,1) > N
        GADec = GADec(1:N,:);
    end

    OffDec = DEDec;
    nGA = floor(N/2);
    if nGA > 0
        gaSite = randperm(N,nGA);
        OffDec(gaSite,:) = GADec(gaSite,:);
    end
end
