classdef DCF4 < handle
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
        function obj = DCF4(ft,nt,gen,Dec,preEvolution,AddProper)
            if nargin > 0
                [N,D] = size(Dec);
                obj(1,N) = DCF4;
              %% calculate objective value
              if gen < preEvolution
                t = 0;
              else
                t = floor((gen-preEvolution)/ft+1)* 1/nt;
              end
                G = cos(pi*t);
                g = 1 + sum(Dec(:,2:end).^2 ,2); 
                Obj(:,1) = g.*Dec(:,1);
                Obj(:,2) = g.*sqrt(1-Dec(:,1).^2);

                %% calculate constraint violations
                c11 = 3 - G - exp(Obj(:,1)) - 0.3*sin(4*pi*Obj(:,1)) - Obj(:,2);
                c12 = 4.1 - (1 + 0.3*Obj(:,1).^2 + Obj(:,1)) - 0.3*sin(4*pi*Obj(:,1)) - Obj(:,2);
                Con = c11.*c12;
                Obj = Obj + 3*t;
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
                pf=[]; pf2=[]; 
                pf = V;
                t  = (i-1) / nt;
                G = cos(pi*t);
                c11 = 3 - G - exp(pf(:,1)) - 0.3*sin(4*pi*pf(:,1)) - pf(:,2);
                c12 = 4.1 - (1 + 0.3*pf(:,1).^2 + pf(:,1)) - 0.3*sin(4*pi*pf(:,1)) - pf(:,2);
                Con = c11.*c12;
                pf(Con>0,:) = [];
                syms x y
                s=solve(3 - G - exp(x) - 0.3*sin(4*pi*x) - y==0,y==0,x,y);
                X=min(double(s.x));
                pf2(:,1) = (0:1/(500-1):X)';
                pf2(:,2) = 3 - G - exp(pf2(:,1)) - 0.3*sin(4*pi*pf2(:,1));
                OC = pf2(:,1).^2 + pf2(:,2).^2 -1;
                pf2(OC<0,:) = [];
                pf = [pf;pf2];
                pf(NDSort(pf,1)~=1,:) = [];
%                 P(i) = struct('PF',pf);
                P(i) = struct('PF',pf+3*t);
            end 
        end
    end
end