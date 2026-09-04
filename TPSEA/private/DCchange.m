function change = DCchange(Population,Problem,ft,nt,gen,preEvolution)
% Dectect the change
    N = size(Population,2);
    X = randperm(N);
    Current = Population(X(1:round(0.1*N)));

    Dec = Current.decs;
    [NewCurrent,~] = Individual(Problem,ft,nt,gen,[],preEvolution,Dec);
    obj1 = Current.objs;
    obj2 = NewCurrent.objs;
    change1 = sum(sum(abs(obj1-obj2),2));
    
    con1 = Current.cons;
    con2 = NewCurrent.cons;
    change2 = sum(sum(abs(con1-con2),2));
    change = change1 + change2;
end