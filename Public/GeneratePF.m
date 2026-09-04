function PF = GeneratePF(Problem,ft,nt,maxgen,preEvolution)
%GENERATEPF Generate reference Pareto fronts for DCF1--DCF10.

    switch Problem
        case 'DCF1',  PF = DCF1.PF(ft,nt,maxgen,preEvolution);
        case 'DCF2',  PF = DCF2.PF(ft,nt,maxgen,preEvolution);
        case 'DCF3',  PF = DCF3.PF(ft,nt,maxgen,preEvolution);
        case 'DCF4',  PF = DCF4.PF(ft,nt,maxgen,preEvolution);
        case 'DCF5',  PF = DCF5.PF(ft,nt,maxgen,preEvolution);
        case 'DCF6',  PF = DCF6.PF(ft,nt,maxgen,preEvolution);
        case 'DCF7',  PF = DCF7.PF(ft,nt,maxgen,preEvolution);
        case 'DCF8',  PF = DCF8.PF(ft,nt,maxgen,preEvolution);
        case 'DCF9',  PF = DCF9.PF(ft,nt,maxgen,preEvolution);
        case 'DCF10', PF = DCF10.PF(ft,nt,maxgen,preEvolution);
        otherwise
            error('Only DCF1--DCF10 are included in this package.');
    end
end
