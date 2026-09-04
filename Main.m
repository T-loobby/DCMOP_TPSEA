function Main(smokeTest)
%MAIN Run TPSEA on the DCF1--DCF10 benchmark suite.
%   Main() runs every DCF problem once with the standard configuration.
%   Main(true) runs a short DCF1 smoke test for installation checking.

    if nargin < 1
        smokeTest = false;
    end

    rootDir = fileparts(mfilename('fullpath'));
    addpath(genpath(rootDir));

    problems = arrayfun(@(i) sprintf('DCF%d',i),1:10,'UniformOutput',false);
    ft = 20;
    nt = 5;
    preEvolution = 60;
    maxgen = preEvolution + 30*ft;
    populationSize = 100; % TPSEA's history model assumes N = 100.
    runTimes = 1;

    if smokeTest
        problems = {'DCF1'};
        ft = 5;
        preEvolution = 5;
        maxgen = preEvolution + 3*ft;
    end

    resultDir = fullfile(rootDir,'Results','TPSEA');
    if ~exist(resultDir,'dir')
        mkdir(resultDir);
    end

    for p = 1:numel(problems)
        problem = problems{p};
        PF = GeneratePF(problem,ft,nt,maxgen,preEvolution);

        for runNumber = 1:runTimes
            rng(runNumber,'twister');
            [initialPopulation,boundary,~,N] = Individual( ...
                problem,ft,nt,runNumber,populationSize,preEvolution);

            fprintf('TPSEA | %s | run %d/%d\n',problem,runNumber,runTimes);
            startTime = tic;
            finalPopulation = TPSEA(0,problem,initialPopulation,boundary, ...
                N,maxgen,ft,nt,PF,preEvolution);
            finalPopulation(end,:) = [];
            runtime = toc(startTime);
            metric = MetricCal(finalPopulation,PF,runtime);

            outputFile = fullfile(resultDir,sprintf('%s_run%02d.mat',problem,runNumber));
            save(outputFile,'problem','ft','nt','preEvolution','maxgen', ...
                'N','runtime','metric','finalPopulation');
            fprintf('Saved %s | IGD %.6g | HV %.6g | %.2f s\n', ...
                outputFile,metric.IGD,metric.HV,runtime);
        end
    end
end
