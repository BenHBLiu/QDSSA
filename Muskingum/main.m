clear
close all
clc
format short



N = 100;
T = 500;
D = 5;
fobj = @(x) fitness(x); 

[f_QDSSA, x_QDSSA, B_QDSSA] = QDSSA(N, T, lb, ub, D, fobj);
[f_QSSA, x_QSSA, B_QSSA] = QSSA(N, T, lb, ub, D, fobj);
[f_DE, x_DE, B_DE] = LSHADE(N, T, lb, ub, D, fobj);
[f_DSSA, x_DSSA, B_DSSA] = DSSA(N, T, lb, ub, D, fobj);
[f_BKA, x_BKA, B_BKA] = BKA(N, T, lb, ub, D, fobj);
[f_SAO, x_SAO, B_SAO] = SAO(N, T, lb, ub, D, fobj);
[f_WSO, x_WSO, B_WSO] = WSO(N, T, lb, ub, D, fobj);
[f_SSA, x_SSA, B_SSA] = SSA(N, T, lb, ub, D, fobj);
[f_GA, x_GA, B_GA] = GA(N, T, lb, ub, D, fobj);
[f_PSO, x_PSO, B_PSO] = PSO(N, T, lb, ub, D, fobj);

[SSQ_PSO, SAD_PSO, MARE_PSO] = calSSM(x_PSO);
[SSQ_GA, SAD_GA, MARE_GA] = calSSM(x_GA);
[SSQ_SSA, SAD_SSA, MARE_SSA] = calSSM(x_SSA);
[SSQ_WSO, SAD_WSO, MARE_WSO] = calSSM(x_WSO);
[SSQ_SAO, SAD_SAO, MARE_SAO] = calSSM(x_SAO);
[SSQ_BKA, SAD_BKA, MARE_BKA] = calSSM(x_BKA);
[SSQ_DSSA, SAD_DSSA, MARE_DSSA] = calSSM(x_DSSA);
[SSQ_DE, SAD_DE, MARE_DE] = calSSM(x_DE);
[SSQ_QSSA, SAD_QSSA, MARE_QSSA] = calSSM(x_QSSA);
[SSQ_QDSSA, SAD_QDSSA, MARE_QDSSA] = calSSM(x_QDSSA);


parameter = [x_PSO; x_GA; x_SSA; x_WSO; x_SAO;...
             x_BKA; x_DSSA; x_DE; x_QSSA; x_QDSSA;];
         
statistic = [SSQ_PSO, SAD_PSO, MARE_PSO;...
             SSQ_GA, SAD_GA, MARE_GA;...
             SSQ_SSA, SAD_SSA, MARE_SSA;...
             SSQ_WSO, SAD_WSO, MARE_WSO;...
             SSQ_SAO, SAD_SAO, MARE_SAO;...
             SSQ_BKA, SAD_BKA, MARE_BKA;...
             SSQ_DSSA, SAD_DSSA, MARE_DSSA;...
             SSQ_DE, SAD_DE, MARE_DE;...
             SSQ_QSSA, SAD_QSSA, MARE_QSSA;...
             SSQ_QDSSA, SAD_QDSSA, MARE_QDSSA];

value = [parameter statistic];

xlswrite('value3.xlsx', value);


