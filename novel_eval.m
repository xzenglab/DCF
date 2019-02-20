% eval on katz and catapult 
% eval cdf of katz
try
    load('novel_testing2014.mat');
    load('novel_training.mat');
    test1 = novel_testing2014 - novel_training;
    test1(test1 == -1) = 0;
%   load('emerge2014.mat');
%     load('emerge2014.mat');
%     catapult_novel(1,0,9,10,10);
%     catapult(1,0,9,10,30);
	% katz_novel(1);
% 	load('ScoreMatrix_Katz_novel.mat')
% 	novelrate_katz = cdf(full(emerge2014), ScoreMatrix_KA, 100);
% 	save 'novelrate_katz.mat' 'novelrate_katz'
% 	send_mail_upon_finished('evaluation of katz novel finished', 'finished','18850544602@163.com'); 

    % eval cdf of catapult
	% catapult_novel(1,0,9,10,30);
%     IMC_novel(200,100,100,10,0.01);
    
    % plot recall 
    x = 1:100;
    figure1 = figure;
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    box on;
    set(gcf, 'position', [0 0 1200 900]);
    lineWidth = 2;
    
    % plot Katz 
    load('ScoreMatrix_Katz_STRING_novel.mat');
    novelrate_katz = recall(full(test1), ScoreMatrix_KA, 100) .* 100;
    plot(x, novelrate_katz(1:100), '--','LineWidth',lineWidth,'Color',[0.537254929542542 0.0705882385373116 0.537254929542542],'linesmoothing','on');
    
    % plot CATAPULT 
    load('ScoreMatrix_CATAPULT_STRING_novel.mat');
    novelrate_catapult = recall(full(test1), ScoreMatrix, 100) .* 100;
    plot(x, novelrate_catapult(1:100), '--','LineWidth',lineWidth,'Color',[0.392156862745098 0.788235294117647 0.372549019607843], 'linesmoothing','on');
    
    
    % plot IMC
    load('ScoreMatrix_IMC_with_STRING_novel.mat');
    novelrate_IMC = recall(full(test1), ScoreMatrix, 100) .* 100;
    plot(x, novelrate_IMC(1:100), '-','LineWidth',lineWidth,'Color',[0.23921568627451 0.43921568627451 0.8],'linesmoothing','on');
    
    % plot DCF
%     load('ScoreMatrix_DCF_novel_alpha.mat');
%     novelrate_DCF = recall(full(test1), ScoreMatrix, 100);
%     plot(x, novelrate_DCF(1:100), 'k-','LineWidth',lineWidth);
    
    % plot DCF-UB
%     load('ScoreMatrix_DCF-UB_novel.mat')
%     novelrate_DCF_UB = recall(full(test1), ScoreMatrix, 100);
%     plot(x, novelrate_DCF_UB(1:100), 'r-','LineWidth',lineWidth);
    
     % plot DCF-UB with STRING
    load('ScoreMatrix_DCF-UB_novel_with_STRING');
    novelrate_DCF_UB_with_STRING = recall(full(test1), ScoreMatrix, 100) .* 100;
    plot(x, novelrate_DCF_UB_with_STRING(1:100), '-','LineWidth', lineWidth,'Color', [0.827450980392157 0.341176470588235 0.341176470588235],'linesmoothing','on');

    % plot DCF with STRING
    load('ScoreMatrix_DCF_novel_with_STRING');
    novelrate_DCF_with_STRING = recall(full(test1), ScoreMatrix, 100) .* 100;
    plot(x, novelrate_DCF_with_STRING(1:100), '-','LineWidth',lineWidth,'Color',[1 1 1]);
    
   
    
    % legend and save
    legend( 'Katz','Catapult', 'IMC', 'DCF-UB', 'DCF','Location','best','EdgeColor','w');
	save('novelrate_recall_all.mat','novelrate_catapult','novelrate_katz', 'novelrate_IMC','novelrate_DCF_with_STRING', 'novelrate_DCF_UB_with_STRING');
	send_mail_upon_finished('evaluation for novel diseases finished', 'finished', '18850544602@163.com');
    set(axes1,'FontSize',12,'LineWidth',1.5,'TickLength',[0 0],...
    'TitleFontSizeMultiplier',0.01,'TitleFontWeight','normal');
    xlabel('Number of genes looked at');
    ylabel('P(hidden gene among genes looked at)');

catch ME 
    rethrow(ME)
   	disp('Error');
    	% please enter create you own check conditions .
	% send_mail_upon_finished('something wroing during train or evaluation', ME.message, 'your email address');
	
end
