
%% how well does the vicon and video agree
clear
conditions = {'FaNoOcc_1','FaNoOcc_2','FaNoOcc_3','FaOcc_1','FaOcc_2','FaOcc_3'};
%
subset = @(x) x{1};

cond_name = cellfun(@(x) subset(strsplit(x,'_')),conditions,'UniformOutput',0);
cond_name = strjoin(unique(cond_name),'_');


varargin_table = table();
varargin_table.behavior ='video_openpose_landmarks_manual_cleaned_in_vicon_units';
varargin_table.analysis_type = 'nancorr_';

varargin_table.corr_fun = @(x,y) nancorr(x,y);

[results_corr_avg_cond_body,~,results_corr] = IBS_vicon_video_permute_corr(conditions,varargin_table);

%%
varargin_table = table();
varargin_table.behavior ='video_openpose_landmarks_manual_cleaned_in_vicon_units';
varargin_table.analysis_type = 'ICC_';

varargin_table.corr_fun = @(x,y) icc21([x' y']);

[results_corr_avg_cond_body,~,results_corr] = IBS_vicon_video_corr(conditions,varargin_table);


%%

varargin_table = table();
varargin_table.behavior ='video_openpose_landmarks_manual_cleaned_in_vicon_units';
varargin_table.analysis_type = 'ICC_';

varargin_table.corr_fun = @(x,y) icc21([x' y']);

[results_corr_avg_cond_body,~,results_corr] = IBS_vicon_video_permute_corr(conditions,varargin_table);

%%


%%

distance = IBS_vicon_video_get_corresp_manipulations({'no_manipulation','distance','distance_025'});




%% how much do the dyads move?
[results_mag_mov_avg_cond_body,body_parts] = IBS_vicon_mov_size(conditions);


