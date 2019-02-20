function mysdae(data, rank, network,filename)
para_gpuDevice=4;
length = size(data, 2);
if  network == 0
    para_layers=[length 1000 500 200 rank];
    para_blayers=[1 1 1 1 1];
elseif network == 1
    para_layers=[length 200 rank];
    para_blayers=[1 1 1];
end
para_lv=10;
para_lu=1;
para_ln=1e3;
para_pretrain=0;
% para_save=999;
para_folder=45;
para_dropout=0.1;
para_from=1;
para_sdae_n_epoch=500;
para_save_lag=501;
fprintf('The pid is: %d\n',feature('getpid'));
sdae_worker(para_gpuDevice,para_layers,para_blayers,...
    para_lv,para_lu,para_ln,para_pretrain,...
    para_folder,para_dropout,para_from,para_sdae_n_epoch,...
    para_save_lag, data, network,filename);
%%exit;
