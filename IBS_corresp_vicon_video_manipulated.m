
%% how well does the vicon and video agree
clear
conditions = {'FaNoOcc_1','FaNoOcc_2','FaNoOcc_3','FaOcc_1','FaOcc_2','FaOcc_3'};
%
subset = @(x) x{1};

cond_name = cellfun(@(x) subset(strsplit(x,'_')),conditions,'UniformOutput',0);
cond_name = strjoin(unique(cond_name),'_');
%%
varargin_table = table();

manipulation = 'distance';
manipulation_dst_05 = 'distance_05';
manipulation_dst_025 = 'distance_025';



manipulations = {manipulation,manipulation_dst_05,manipulation_dst_025,};

for manipulation_no = 1:length(manipulations)
    cur_manipulation = manipulations{manipulation_no};
    
    varargin_table.manipulation = cur_manipulation;
    
    [results_corr_avg_cond_body] = IBS_vicon_manipulated_video_corr(conditions,varargin_table);

end