function [POP,POPU,POF,V]=peizhun(Lower,Upper,AllPop,AllPopU,N)
Pop_1=AllPop(end-99:end);           %最近的一个环境的种群
Pop_2=AllPop(end-199:end-100);      %前一个环境的种群
Pop_1U=AllPopU(end-99:end);         %最近的一个环境的无约束种群
Pop_2U=AllPopU(end-199:end-100);    %前一个环境的无约束种群
Lower = repmat(Lower, N,1);
Upper = repmat(Upper, N, 1);
V=1;            %判断无约束Perato是否改变

POP=Pop_1.decs+(Pop_1.decs-Pop_2.decs);
POPU=Pop_1U.decs+(Pop_1U.decs-Pop_2U.decs);

%增加随机扰动（多样性：一部分我们也不预测）
randsize=randperm(N);
rand20=randsize(1:0.2*N);
randa=rand20(1:1/2*end);

POP(randa,:)=Pop_1(randa).decs;
POPU(randa,:)=Pop_1U(randa).decs;

%边界约束处理，让所有跨界的，重新赋予随机值
Dec = unifrnd(repmat(Lower,round(N),1),repmat(Upper,round(N),1));
POP(find(POP<Lower|POP>Upper))=Dec(find(POP<Lower|POP>Upper));
POPU(find(POPU<Lower|POPU>Upper))=Dec(find(POPU<Lower|POPU>Upper));

%预测Perato前沿
Pop_1=SortedPopobj(Pop_1);
Pop_2=SortedPopobj(Pop_2);
[Nobj, Fobj] = obtain_matrix_mean(Pop_1, Pop_2, N);
[curve_pred, ~] = predict_tps_deformation(Fobj,Nobj, 0.1);
% [curve_pred, ~] = predict_tps_deformation(Pop_2.objs,Pop_1.objs, 0.1);
[POF,V] = select_points(curve_pred, 100,Pop_1.objs,V);%Pop_1.objs

