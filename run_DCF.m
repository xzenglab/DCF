
fname = 'GeneGene_SDAE.mat';
%  indicator ('microarray', 'OMIM', 'orthologous', 'genes', 'phenotypes'), correspoding to 16,8,4,2,1

indicator = 0;
notes = sprintf('start evaluating novel disease by using DCF');
disp(notes);
ind = 1;
repeatTime = 20;
% load test file; this indicate disease appeared after 2014
%     load('emerge2014.mat');
load('splitsUniform.mat')
%     load('novel_ScoreMatrix_alpha=0.96_lambda=0.01.mat');
test = splits{1};
x = 1:100;
labels = [];
figure;
hold on;

% run DCF and evaluate 
% last bit control whether to add STRING
%     scores = DCF_mul_vieww(ind,200,100,100,2e-6,2e-2,fname,indicator, 0,1);
%     rate_DCF = recall(full(test1), scores, 100);
%     disp(rate_DCF(100));
%     plot(x, rate_DCF, '--', 'LineWidth', 2);
%     
% run DCF with STRINg and evaluate 
% 2e-6, 2e-2
% best 0.01, 0.05
for lambda = [1e-4 1e-3 1e-2 1e-6 1e-5 ]
    for alpha = [1e-2 2e-2 5e-2 1e-1 ] 
        scores = DCF_mul_views(ind,200,100,100,lambda,alpha,fname,indicator);
        % recall
        recall_DCF = recall(full(test), scores, 100);

        plot(x, recall_DCF, '--', 'LineWidth', 2);

        % cdf 
        cdf_DCF = cdf(full(test), scores, 100);

        plot(x, cdf_DCF, '--', 'LineWidth', 2);

        legend('recall','cdf');
        % res str
        res_str = sprintf('cdf: %.4f, recall: %.4f', cdf_DCF(100), recall_DCF(100));
        disp(res_str);
        send_mail_upon_finished('DGCCA evaluation finished', res_str, '18850544602@163.com');
    end 
end
    
    