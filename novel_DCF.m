function ScoreMatrix = novel_DCF(k, featureRank, networkRank, lambda, alpha, fname, indicators, frac)
% INPUT PARAMS:
% split       : One of the k-fold splits to be hidden.  %% 留一验证法, 记住用于训练的集合
% k           : Rank of the parameter matrix Z. %% 作为参数的矩阵的秩
% featureRank : Number of compressed dimensions (microarray gene expression data). %% 希望得到的数据降维之后的维度
% networkRank : Number of latent features from similarity networks. %% 相似性网络的维度
% loss	    : loss function to be used (10 - for squared loss) %% 损失函数的参数
% lambda      : Regularization parameter l in the objective. %% ?
% alpha       : the parameter r for PU formulation  %% ?
% mode : specify bootstrapping features in every bit, only one bit can be
%           1, in the sequence ('microarray', 'OMIM', 'orthologous', 'genes', 
%               'phenotypes'), correspoding to 16,8,4,2,1
% frac : select fraction of origin features

%%%% OUTPUT:
%%% ScoreMatrix : Matrix of size # genes x # diseases; Each column gives the scores estimated by IMC for the corresponding OMIM disease.	 %% 评分矩阵

%%% Data files are provided in the same directory.

% load data from file
tic
genesPhenes = load ('genesPhenes.mat');  %% 基因表型
geneFeatures = load ('geneFeatures.mat');	%% 基因特征
clinicalFeatures = load ('clinicalFeaturesTFidf.mat'); %% 临床特征
load('novel_training.mat')
% specify filenames to save .mat files


POSTFIX = '.mat';

BOOTSTRP_FEA_NAME = '';
BOOTSTRP_FEA = [];

% 取第一个物种基因集作为训练集
genePheneTraining = novel_training;

%genePheneTraining(:,25) = 0;
% if (split_index > 0) %% else use all data to m ake predictions.
%     splits = load('splitsUniform.mat');
%     numPhenes = size(splits.splits{split_index},2); %%返回splits 中维度2的长度,得到基因表型的维度
%     genePheneTraining = genePheneTraining(:,1:numPhenes) - splits.splits{split_index}; %% 为什么要做减法?
% end
features = [];  %% 行特征 -- gene
col_features = []; %% 列特征 -- disease

% assign mode to dropping feature
% dropping_feature = @(feature) rand(size(feature));
dropping_feature = @(feature) zeros(size(feature));

% 分别对 Gene 和 Disease 进行降维
if (featureRank == Inf)	%%如果需要的维度很大, 那么不需要降维, 直接返回默认矩阵
    features = geneFeatures.GeneFeatures; %% SVD 的两个分解矩阵?
    col_features = clinicalFeatures.colFeatures;	%%
elseif (featureRank > 0)
    disp ('Reducing feature dimensionality..');
%     Reducing dimensionality of microarray
%     MICROARRAY_NAME = 'GeneFeaSDAE.mat';
    MICROARRAY_NAME = 'microarray_features_100.mat';
      
%     geneFeatures.GeneFeatures = normalization(geneFeatures.GeneFeatures); %% 归一化
    
%     if ~exist(microarray_filename,'file')   % 2 specified file type
%         mysdae(geneFeatures.GeneFeatures, featureRank, 0, MICROARRAY_NAME);
%     end
    load(MICROARRAY_NAME);
    features = H;
    if bitand(indicators,16) 
        disp('erasing microarray...')
%         indices = bootstrp_features(frac, size(col_features));
%         features = features(:,indices);
        features = dropping_feature(features);
%         features = rand(size(features));
    end
    clear H;
    
    
    % Reducing dimensionality of OMIM pages
    col_featuresFilename = 'pre_colFeatures.mat';
    if ~exist(col_featuresFilename,'file')
        [U,S] = svds(sparse(clinicalFeatures.F(1:numPhenes,:)), featureRank); %% 使用 spaese 创建稀疏矩阵
        col_features = U;
        save(col_featuresFilename,'colFeatures')
        if bitand(indicators,8) 
              disp('erasing OMIM features..')
%             indices = bootstrp_features(frac, size(col_features));
%             col_features = col_features(:,indices);
              col_features = dropping_feature(col_features);
        end
    else
        load(col_featuresFilename);
        col_features = colFeatures;
        if bitand(indicators,8) 
%             indices = bootstrp_features(frac, size(col_features));
%             col_features = col_features(:,indices);
              col_features = dropping_feature(col_features);
        end
    end
end

if (networkRank > 0)
    disp ('Reducing gene network dimensionality..');
    % Reducing dimensionality of gene-gene interaction network
%    network_filename = 'GeneGeneHsp1q2.mat';
%    network_features = load(network_filename);
%    features = [features network_features.features(1:genesPhenes.numGenes,:) ];
     network_features = load(fname);

     if bitand(indicators,4) 
        disp('erasing HumanNet...') 
%       indices = bootstrp_features(frac, size(network_features.features));
%       network_features.features = network_features.features(:,indices);    
        network_features.H = dropping_feature(network_features.H);
     end
     features = [features network_features.H(1:genesPhenes.numGenes,:) ];
     
     string = load('STRING_discreted.mat');
     features = [features string.features];
     % Reducing dimensionality of orthologous phenotypes, 循环处理 GenePhene
   % 中的每个 Cell 同源基因
    PHENOTYPES_NAME = 'GP_features_100.mat';
  
    if ~exist(PHENOTYPES_NAME,'file')
        GP_sp = [];
            for sp=2:numel(genesPhenes.GenePhene)
               GP = genesPhenes.GenePhene{sp};
               for i=1:size(GP,2)
                   if(sum(GP(:,i)) > 0)
                       GP(:,i) = GP(:,i)/norm(GP(:,i));
                   end
               end
               GP_sp = [GP_sp GP];
            end
         mysdae(full(GP_sp), networkRank, 1, PHENOTYPES_NAME);
    end
    sdae_data = load(PHENOTYPES_NAME);
    sdae_H = sdae_data.H;
    if bitand(indicators,2)
%        indices = bootstrp_features(frac, size(sdae_H));
%        sdae_H = sdae_H(:,indices);
         sdae_H = dropping_feature(sdae_H);
    end
    
    features = [features sdae_H];
    
    clear sdae_H
    %     
    % Reducing dimensionality of  disease similarities network
    uFilename = 'HN2vec_U.mat';
    if ~exist(uFilename,'file')
        [U,S] = svds(genesPhenes.PhenotypeSimilaritiesLog(1:numPhenes,1:numPhenes), networkRank);
        save(uFilename, 'U');
        clear U S
    end
    load(uFilename);
    if bitand(indicators,1)
%        indices = bootstrp_features(frac, size(U));
%        U = U(:,indices);
         U = dropping_feature(U);
    end
    
    col_features = [col_features U];
    clear U
    
end

for i=1:size(features,2)
    if (norm(features(:,i)) > 0)
        features(:,i) = features(:,i)/norm(features(:,i)); %% 返回矢量或者矩阵的2-范数
    end
end
for i=1:size(col_features,2)
    if (norm(col_features(:,i))>0)
        col_features(:,i) = col_features(:,i)/norm(col_features(:,i));
    end
end

% If no features are supplied (i.e. networkRank & featureRank are 0), use default matrix completion.
if (isempty(features))
    features = eye(genesPhenes.numGenes);
end
if (isempty(col_features))
    col_features = eye(numPhenes);
end
fprintf('Using %d row, %d col features.\n', size(features,2), size(col_features,2));

matrix_filename = 'HN2vec_P.mat';
% [W,H,~,~] = imfTrain(sparse(double(genePheneTraining)), sparse(double(features)), sparse(double(col_features)),  ['-n 16 -t 15 -T 5 -g 20 -p 3 -a 0 -s 11 -k ', num2str(k), ' -l ', num2str(lambda),  ' -r ', num2str(alpha)]);
 [W, H, ~] = train_mf(sparse(double(genePheneTraining)), sparse(double(features)), sparse(double(col_features)), [' -l ' num2str(lambda) ' -k ' num2str(k) ' -t 10 ' ' -s ' num2str(10)]); 
save(matrix_filename,'W','H');
save('features.mat', 'features');
save('colFeatures.mat','col_features');
% for DCF-UB
ScoreMatrix = features * W' * H * col_features';
% for DCF
% ScoreMatrix = features * W *H' * col_features';
% save 'novel_DCF.mat' 'ScoreMatrix'
save 'novel_DCF_UB.mat' 'ScoreMatrix'


toc
end
