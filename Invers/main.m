% warning off
clear
close all
clc
format short e
N=50;
MaxIT=500;

Function_name='F29';
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

runtimes = 3;

for j=1:runtimes
    j
    disp('CQPSO');
    [f_gbest(j),best(j,:),gbest(j,:)]=CQPSO(N,MaxIT,lb,ub,dim,fobj);
    disp('QDSSA');
    [QDSSA_Fbest(j),QDSSA_Best(j,:),QDSSA_xbest(j,:)]=QDSSA(N,MaxIT,lb,ub,dim,fobj);
    disp('QSSA');
    [QSSA_Fbest(j),QSSA_Best(j,:),QSSA_xbest(j,:)]=QSSA(N,MaxIT,lb,ub,dim,fobj);
    disp('DE');
    [DE_Fbest(j),DE_Best(j,:),DE_xbest(j,:)]=LSHADE(N,MaxIT,lb,ub,dim,fobj);
    disp('DSSA');
    [DSSA_Fbest(j),DSSA_Best(j,:),DSSA_xbest(j,:)]=DSSA(N,MaxIT,lb,ub,dim,fobj);
    disp('BKA');
    [BKA_Fbest(j),BKA_Best(j,:),BKA_xbest(j,:)]=BKA(N,MaxIT,lb,ub,dim,fobj);
    disp('SAO');
    [SAO_Fbest(j),SAO_Best(j,:),SAO_xbest(j,:)]=SAO(N,MaxIT,lb,ub,dim,fobj);
    disp('WSO');
    [WSO_Fbest(j),WSO_Best(j,:),WSO_xbest(j,:)]=WSO(N,MaxIT,lb,ub,dim,fobj);
    disp('SSA');
    [SSA_Fbest(j),SSA_Best(j,:),SSA_xbest(j,:)]=SSA(N,MaxIT,lb,ub,dim,fobj);
    disp('GA');
    [GA_Fbest(j),GA_Best(j,:),GA_xbest(j,:)]=GA(N,MaxIT,lb,ub,dim,fobj);
    disp('PSO');
    [PSO_Fbest(j),PSO_Best(j,:),PSO_xbest(j,:)]=PSO(N,MaxIT,lb,ub,dim,fobj);
end

[~,index_CQPSO]=min(f_gbest);
[~,index_QDSSA]=min(QDSSA_Fbest);
[~,index_QSSA]=min(QSSA_Fbest);
[~,index_DE]=min(DE_Fbest);
[~,index_DSSA]=min(DSSA_Fbest);
[~, index_BKA] = min(BKA_Fbest);
[~, index_SAO] = min(SAO_Fbest);
[~, index_WSO] = min(WSO_Fbest);
[~, index_SSA] = min(SSA_Fbest);
[~, index_GA] = min(GA_Fbest);
[~, index_PSO] = min(PSO_Fbest);



disp('----------------------------------------')
display(['The best solution obtained by CQPSO is : ', num2str(gbest(index_CQPSO,:))]);
display(['The best optimal value of the objective function found by CQPSO is : ', num2str(f_gbest(index_CQPSO))]);
disp('----------------------------------------')

display(['The best solution obtained by QDSSA is : ', num2str(QDSSA_xbest(index_QDSSA,:))]);
display(['The average value of the best solution obtained by QDSSA is : ', num2str(mean(QDSSA_xbest(index_QDSSA,:)))]);
display(['The best optimal value of the objective function found by QDSSA is : ', num2str(QDSSA_Fbest(index_QDSSA))]);
disp('----------------------------------------')

display(['The best solution obtained by QSSA is : ', num2str(QSSA_xbest(index_QSSA,:))]);
display(['The average value of the best solution obtained by QSSA is : ', num2str(mean(QSSA_xbest(index_QSSA,:)))]);
display(['The best optimal value of the objective function found by QSSA is : ', num2str(QSSA_Fbest(index_QSSA))]);
disp('----------------------------------------')

display(['The best solution obtained by DE is : ', num2str(DE_xbest(index_DE,:))]);
display(['The average value of the best solution obtained by DE is : ', num2str(mean(DE_xbest(index_DE,:)))]);
display(['The best optimal value of the objective function found by DE is : ', num2str(DE_Fbest(index_DE))]);
disp('----------------------------------------')

display(['The best solution obtained by DSSA is : ', num2str(DSSA_xbest(index_DSSA,:))]);
display(['The average value of the best solution obtained by DSSA is : ', num2str(mean(DSSA_xbest(index_DSSA,:)))]);
display(['The best optimal value of the objective function found by DSSA is : ', num2str(DSSA_Fbest(index_DSSA))]);
disp('----------------------------------------')

display(['The best solution obtained by BKA is : ', num2str(BKA_xbest(index_BKA,:))]);
display(['The average value of the best solution obtained by BKA is : ', num2str(mean(BKA_xbest(index_BKA,:)))]);
display(['The best optimal value of the objective function found by BKA is : ', num2str(BKA_Fbest(index_BKA))]);
disp('----------------------------------------')

display(['The best solution obtained by SAO is : ', num2str(SAO_xbest(index_SAO,:))]);
display(['The average value of the best solution obtained by SAO is : ', num2str(mean(SAO_xbest(index_SAO,:)))]);
display(['The best optimal value of the objective function found by SAO is : ', num2str(SAO_Fbest(index_SAO))]);
disp('----------------------------------------')

display(['The best solution obtained by WSO is : ', num2str(WSO_xbest(index_WSO,:))]);
display(['The average value of the best solution obtained by WSO is : ', num2str(mean(WSO_xbest(index_WSO,:)))]);
display(['The best optimal value of the objective function found by WSO is : ', num2str(WSO_Fbest(index_WSO))]);
disp('----------------------------------------')

display(['The best solution obtained by SSA is : ', num2str(SSA_xbest(index_SSA,:))]);
display(['The average value of the best solution obtained by SSA is : ', num2str(mean(SSA_xbest(index_SSA,:)))]);
display(['The best optimal value of the objective function found by SSA is : ', num2str(SSA_Fbest(index_SSA))]);
disp('----------------------------------------')

display(['The best solution obtained by GA is : ', num2str(GA_xbest(index_GA,:))]);
display(['The average value of the best solution obtained by GA is : ', num2str(mean(GA_xbest(index_GA,:)))]);
display(['The best optimal value of the objective function found by GA is : ', num2str(GA_Fbest(index_GA))]);
disp('----------------------------------------')

display(['The best solution obtained by PSO is : ', num2str(PSO_xbest(index_PSO,:))]);
display(['The average value of the best solution obtained by PSO is : ', num2str(mean(PSO_xbest(index_PSO,:)))]);
display(['The best optimal value of the objective function found by PSO is : ', num2str(PSO_Fbest(index_PSO))]);
disp('----------------------------------------')

hold on
semilogy(best(index_CQPSO, :), 'k-', 'LineWidth', 2); % CQPSO - 黑色实线
semilogy(QDSSA_Best(index_QDSSA, :), 'b:', 'LineWidth', 2); % QDSSA - 蓝色点线
semilogy(QSSA_Best(index_QSSA, :), 'm--', 'LineWidth', 2); % QSSA - 紫色虚线
semilogy(DE_Best(index_DE, :), 'r-.', 'LineWidth', 2); % DE - 红色点划线
semilogy(DSSA_Best(index_DSSA, :), 'c-', 'LineWidth', 2); % DSSA - 青色实线
semilogy(BKA_Best(index_BKA, :), 'y--', 'LineWidth', 2); % BKA - 黄色虚线
semilogy(SAO_Best(index_SAO, :), 'g:', 'LineWidth', 2); % SAO - 绿色点线
semilogy(WSO_Best(index_WSO, :), 'k-.', 'LineWidth', 2); % WSO - 黑色点划线
semilogy(SSA_Best(index_SSA, :), 'b-', 'LineWidth', 2); % SSA - 蓝色实线
semilogy(GA_Best(index_GA, :), 'r--', 'LineWidth', 2); % GA - 红色虚线
semilogy(PSO_Best(index_PSO, :), 'm:', 'LineWidth', 2); % PSO - 紫色点线

xlabel('Iteration')
ylabel('fitness')
legend('C-QPSO', 'QDSSA', 'QSSA', 'DE', 'DSSA', 'BKA', 'SAO', 'WSO', 'SSA', 'GA', 'PSO');
