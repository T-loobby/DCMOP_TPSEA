classdef DCF6 < handle
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
        function obj = DCF6(ft,nt,gen,Dec,preEvolution,AddProper)
            if nargin > 0
                [N,D] = size(Dec);
                obj(1,N) = DCF6;
              %% calculate objective value
              if gen < preEvolution
                t = 0;
              else
                t = floor((gen-preEvolution)/ft+1)* 1/nt;
              end
                G = sin(0.5*pi*t); 
                g = 1 + sum((Dec(:,2:end)-G).^2 ,2); 
                Obj(:,1) = g.*Dec(:,1);
                Obj(:,2) = g.*sqrt(1-Dec(:,1).^2);

                %% calculate constraint violations
                H = 1;
                if G < 0
                    H = -1;
                end
                ll  = cos(5*atan((Obj(:,2)./Obj(:,1)).^H).^4).^6;
                c12 = 1.1 - (Obj(:,1)./(0.9 + (0.1+0.7*abs(G))*ll)).^2 - (Obj(:,2)./(0.9 + (0.8-0.7*abs(G))*ll)).^2 ;
                Con = - c12;
                Obj = Obj + t;
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
            V = V./repmat(sqrt(sum(V.^2,2)),1,2);
            for i = 1 : ceil((maxgen-preEvolution)/ft+1)
                pf=[]; 
                pf = V;
                t  = (i-1) / nt;
                G = sin(0.5*pi*t);
                H = 1;
                if G < 0
                    H = -1;
                end
                ll  = cos(5*atan((pf(:,2)./pf(:,1)).^H).^4).^6;
                Con = 1.1 - (pf(:,1)./(0.9 + (0.1+0.7*abs(G))*ll)).^2 - (pf(:,2)./(0.9 + (0.8-0.7*abs(G))*ll)).^2;
                pf(Con<0,:) = [];
%                 P(i) = struct('PF',pf);
                P(i) = struct('PF',pf+t);
            end 
        end
    end
end