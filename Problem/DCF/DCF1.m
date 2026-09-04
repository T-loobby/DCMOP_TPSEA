classdef DCF1 < handle
% <problem> <DCF>
% Dynamic MOP benchmark
% ft --- --- frequency of change
% nt --- --- severity of change
    properties(SetAccess = private)
        obj;  % the objective value
        dec;  % the decision vector
        con;  % the constraint violations
        add;  % Additional properties of the individual
    end
    methods
        %% Initialization
        function obj = DCF1(ft,nt,gen,Dec,preEvolution,AddProper)
            if nargin > 0
                [N,D] = size(Dec);
                obj(1,N) = DCF1;
              %% calculate objective value
              if gen < preEvolution
                t = 0;
              else
                t = floor((gen-preEvolution)/ft+1)* 1/nt;
              end
                g = 1 + sum((Dec(:,2:end)-0.5).^2,2); 
                Obj(:,1) = g.*Dec(:,1) ;
                Obj(:,2) = g.*(1-Dec(:,1));

                %% calculate constraint violations
                G = abs(sin(0.5*pi*t));
                c12 = (Obj(:,1)).^2 + (Obj(:,2)).^2 - (0.7+G).^2;
                Con = -c12;
                Obj = Obj + 2*t;
                %% generate population
                for i = 1 : length(obj)
                    obj(i).dec = Dec(i,:);
                    obj(i).obj = Obj(i,:);
                    obj(i).con = Con(i,:);
                end
                if nargin > 5
                    for i = 1 : length(obj)
                        obj(i).add = AddProper(i,:);
                    end
                end
            end
        end
          %% Get the matrix of decision variables of the population
        function value = decs(obj)
        %decs - Get the matrix of decision variables of the population.
            value = cat(1,obj.dec);
        end
        %% Get the matrix of objective values of the population
        function value = objs(obj)
        %objs - Get the matrix of objective values of the population.
            value = cat(1,obj.obj);
        end
        %% Get the matrix of constraint violations of the population
        function value = cons(obj)
        %cons - Get the matrix of constraint violations of the population.
            value = cat(1,obj.con);
        end
        function value = adds(obj,AddProper)
        %adds - Get the matrix of additional properties of the population.
            for i = 1 : length(obj)
                if isempty(obj(i).add)
                    obj(i).add = AddProper(i,:);
                end
            end
            value = cat(1,obj.add);
        end
    end
    methods (Static)   
       %% Sample reference points on Pareto front
        function P = PF(ft,nt,maxgen,preEvolution)
            V(:,1) = (0:1/(500-1):1)';
            V(:,2) = 1 - V(:,1) ;
            V1 = V./repmat(sqrt(sum(V.^2,2)),1,2);
            for i = 1 : ceil((maxgen-preEvolution)/ft+1)
                pf=[]; pf = V;
                pf2 = [];
                t  = (i-1) / nt;
                G = abs(sin(0.5*pi*t));
                Con1 = pf(:,1).^2 + pf(:,2).^2 - (0.7+G).^2;
                pf(Con1<0,:) = [];
                pf2 = (0.7+G)*V1;
                Con2 = pf2(:,1) + pf2(:,2) - 1;
                pf2(Con2<0,:) = [];
                pf = [pf;pf2];
                pf(NDSort(pf,1)~=1,:) = [];
%                 P(i) = struct('PF',pf);
                P(i) = struct('PF',pf+2*t);
            end 
        end
    end
end