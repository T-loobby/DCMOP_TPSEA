function OffDec = TPSEA_DEOperator(Population,boundary,useCons)
% TPSEA_DEOperator - DE/rand-to-best/1 with polynomial mutation.

    if nargin < 3
        useCons = true;
    end

    CR   = 0.9;
    F    = 0.5;
    proM = 1;
    disM = 20;

    Parent = Population.decs;
    [N,D]  = size(Parent);
    OffDec = Parent;

    if N <= 1 || D == 0
        OffDec = PolynomialMutation(OffDec,boundary,proM,disM);
        return;
    end

    if useCons && isprop(Population(1),'cons')
        [FrontNo,~] = NDSort(Population.objs,Population.cons,1);
    else
        [FrontNo,~] = NDSort(Population.objs,1);
    end
    Elite = find(FrontNo == 1);
    if isempty(Elite)
        Elite = 1:N;
    end

    for i = 1 : N
        xelite = Parent(Elite(randi(length(Elite))),:);
        donors = randperm(N);
        donors(donors == i) = [];
        while length(donors) < 3
            donors = [donors,randi(N)]; %#ok<AGROW>
        end
        r1 = donors(1);
        r2 = donors(2);
        r3 = donors(3);

        mutant = Parent(r1,:) + F*(xelite-Parent(r1,:)) + F*(Parent(r2,:)-Parent(r3,:));
        site = rand(1,D) < CR;
        site(randi(D)) = true;
        OffDec(i,site) = mutant(site);
    end

    OffDec = PolynomialMutation(OffDec,boundary,proM,disM);
end

function OffDec = PolynomialMutation(OffDec,boundary,proM,disM)

    [N,D] = size(OffDec);
    Lower = repmat(boundary.lower,N,1);
    Upper = repmat(boundary.upper,N,1);

    OffDec = min(max(OffDec,Lower),Upper);
    Site = rand(N,D) < proM/D;
    mu   = rand(N,D);

    temp = Site & mu <= 0.5;
    OffDec(temp) = OffDec(temp) + (Upper(temp)-Lower(temp)).* ...
        ((2.*mu(temp)+(1-2.*mu(temp)).*(1-(OffDec(temp)-Lower(temp))./ ...
        (Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);

    temp = Site & mu > 0.5;
    OffDec(temp) = OffDec(temp) + (Upper(temp)-Lower(temp)).* ...
        (1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*(1-(Upper(temp)-OffDec(temp))./ ...
        (Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));

    OffDec(OffDec < Lower) = Lower(OffDec < Lower);
    OffDec(OffDec > Upper) = Upper(OffDec > Upper);
end
