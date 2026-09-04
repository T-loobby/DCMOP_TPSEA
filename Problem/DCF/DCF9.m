classdef DCF9 < handle
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
        function obj = DCF9(ft,nt,gen,Dec,preEvolution,AddProper)
            if nargin > 0
                [N,D] = size(Dec);
                obj(1,N) = DCF9;
              %% calculate objective value
              if gen < preEvolution
                t = 0;
              else
                t = floor((gen-preEvolution)/ft+1)* 1/nt;
              end
                G = sin(0.5*pi*t);
                g = 1 + sum((Dec(:,2:end)-G).^2,2); 
                Obj(:,1) = g.*Dec(:,1);
                Obj(:,2) = g.*(1-Dec(:,1));

                %% calculate constraint violations
                W = 0.5*pi*t;                
                c11 = (sin(W)-cos(W))/(sin(W)+cos(W)).* Obj(:,1) - Obj(:,2) + (-2.2*(1-cos(W))+1.3)./(cos(W)+sin(W));
                c12 = (sin(W)-cos(W))/(sin(W)+cos(W)).* Obj(:,1) - Obj(:,2) + (-2.2*(1-cos(W))+1.8)./(cos(W)+sin(W));
                c1 = c11 .* c12;
                c31 = (sin(W)-cos(W))/(sin(W)+cos(W)).* Obj(:,1) - Obj(:,2) + (-2.2*(1-cos(W))+2.6)./(cos(W)+sin(W));
                c32 = (sin(W)-cos(W))/(sin(W)+cos(W)).* Obj(:,1) - Obj(:,2) + (-2.2*(1-cos(W))+3.1)./(cos(W)+sin(W));
                c3 = c31 .* c32;
                Con = -[c1 c3];
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
            for i = 1 : ceil((maxgen-preEvolution)/ft+1)
                pf=[];  pf2=[];  
                t  = (i-1) / nt;
                W = 0.5*pi*t;
                pf = V;
                c11 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf(:,1) - pf(:,2) + (-2.2*(1-cos(W))+1.3)./(cos(W)+sin(W));
                c12 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf(:,1) - pf(:,2) + (-2.2*(1-cos(W))+1.8)./(cos(W)+sin(W));
                c1a = c11 .* c12;
                c31 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf(:,1) - pf(:,2) + (-2.2*(1-cos(W))+2.6)./(cos(W)+sin(W));
                c32 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf(:,1) - pf(:,2) + (-2.2*(1-cos(W))+3.1)./(cos(W)+sin(W));
                c3a = c31 .* c32;
                pf(c1a<0 | c3a<0,:) = [];
                if 0.5*t-floor(0.5*t) >= 0.25 && 0.5*t-floor(0.5*t) <= 0.75
                    pf2(:,1) = [zeros(500,1); 1+V(:,1)];
                    pf2(:,2) = [1+V(:,1); zeros(500,1)];
                    c11 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf2(:,1) - pf2(:,2) + (-2.2*(1-cos(W))+1.3)./(cos(W)+sin(W));
                    c12 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf2(:,1) - pf2(:,2) + (-2.2*(1-cos(W))+1.8)./(cos(W)+sin(W));
                    c1a = c11 .* c12;
                    c31 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf2(:,1) - pf2(:,2) + (-2.2*(1-cos(W))+2.6)./(cos(W)+sin(W));
                    c32 = (sin(W)-cos(W))/(sin(W)+cos(W)).* pf2(:,1) - pf2(:,2) + (-2.2*(1-cos(W))+3.1)./(cos(W)+sin(W));
                    c3a = c31 .* c32;
                    pf2(c1a<0 | c3a<0,:) = [];
                elseif 0.5*t-floor(0.5*t) < 0.25
                    f12 = -(2.2*(1-cos(W))-1.3)/(cos(W)+sin(W));
                    f11 = (2.2*(1-cos(W))-1.3)/(sin(W)-cos(W));
                    f22 = -(2.2*(1-cos(W))-1.8)/(cos(W)+sin(W));
                    f21 = (2.2*(1-cos(W))-1.8)/(sin(W)-cos(W));
                    pf2 = [V.*repmat([f11 f12],size(V,1),1);V.*repmat([f21 f22],size(V,1),1)];
                    c11 = pf2(:,1) + pf2(:,2) -1;
                    pf2(c11<0,:) = [];
                elseif 0.5*t-floor(0.5*t) > 0.75   
                    f12 = -(2.2*(1-cos(W))-3.1)/(cos(W)+sin(W));
                    f11 = (2.2*(1-cos(W))-3.1)/(sin(W)-cos(W));
                    f22 = -(2.2*(1-cos(W))-2.6)/(cos(W)+sin(W));
                    f21 = (2.2*(1-cos(W))-2.6)/(sin(W)-cos(W));
                    pf2 = [V.*repmat([f11 f12],size(V,1),1);V.*repmat([f21 f22],size(V,1),1)];
                    c11 = pf2(:,1) + pf2(:,2) -1;
                    pf2(c11<0,:) = [];
                end
                pf = [pf; pf2];
                pf(NDSort(pf,1)~=1,:) = [];
%                 P(i) = struct('PF',pf);
                P(i) = struct('PF',pf+2*t);
            end 
        end
    end
end