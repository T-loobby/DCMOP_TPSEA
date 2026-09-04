function [Population,boundary,V,N] = Individual(Problem,ft,nt,gen,N,preEvolution,Dec,AddProper)
%INDIVIDUAL Create or evaluate individuals for the DCF benchmark suite.

    validProblems = arrayfun(@(i) sprintf('DCF%d',i),1:10,'UniformOutput',false);
    if ~ismember(Problem,validProblems)
        error('Only DCF1--DCF10 are included in this package.');
    end

    M = 2;
    D = 10;
    if isempty(N)
        V = [];
    else
        [V,N] = UniformPoint(N,M);
    end

    if strcmp(Problem,'DCF2')
        boundary.lower = zeros(1,D);
    else
        boundary.lower = [0,-ones(1,D-1)];
    end
    boundary.upper = ones(1,D);

    problemFunction = str2func(Problem);
    if nargin == 7
        Population = problemFunction(ft,nt,gen,Dec,preEvolution);
    elseif nargin == 8
        Population = problemFunction(ft,nt,gen,Dec,preEvolution,AddProper);
    else
        Dec = unifrnd(repmat(boundary.lower,N,1),repmat(boundary.upper,N,1));
        Population = problemFunction(ft,nt,gen,Dec,preEvolution);
    end
end
