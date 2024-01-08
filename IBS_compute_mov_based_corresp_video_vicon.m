function corr_small_med_large = IBS_compute_mov_based_corresp_video_vicon(Dyad_no,Sub_no,condition)
%IBS_COMPUTE_MOV_BASED_CORRESP_VIDEO_VICON 
%
% SYNOPSIS: IBS_compute_mov_based_corresp_video_vicon
%
% INPUT Function to check the correspondences between vicon and video but based on how much the subjects move
%
% OUTPUT 
%
% REMARKS
%
% created with MATLAB ver.: 9.8.0.1359463 (R2020a) Update 1 on Microsoft Windows 10 Pro Version 10.0 (Build 19042)
%
% created by: Atesh
% DATE: 23-Dec-2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_data_dir = 'D:\\iannettilab_dropbox\\Dropbox\\Koul_Atesh\\IBS\\';
varargin_table = table();
varargin_table.behav_analysis = 'joint';
varargin_table.analysis_sub_type = '_insta_abs_detrend';


body_parts = {'Head','Torso','Left_Shoulder','Left_Elbow','Left_Wrist',...
    'Left_Knee','Left_Feet','Right_Shoulder','Right_Elbow','Right_Wrist','Right_Knee','Right_Feet'};


switch(Sub_no)
    case 0
        vicon_table = table({'S0_Head1','S0_Head2','S0_Head3','S0_Head4'},...
            {'S0_Torso1', 'S0_Torso2', 'S0_Torso3'},...
            {'S0_Left_Arm1'}, {'S0_Left_Arm2'}, {'S0_Left_Arm3'},{ 'S0_Left_Leg1', 'S0_Left_Leg2'},{ 'S0_Left_Leg3'},...
            {'S0_Right_Arm1'},{ 'S0_Right_Arm2'},{ 'S0_Right_Arm3'},{'S0_Right_Leg1', 'S0_Right_Leg2'},{ 'S0_Right_Leg3'},...
            'VariableNames',{'Head','Torso','Left_Shoulder','Left_Elbow','Left_Wrist',...
            'Left_Knee','Left_Feet','Right_Shoulder','Right_Elbow','Right_Wrist','Right_Knee','Right_Feet'});
        
        video_cam_timestamp = 'timestamps_0';
        video_cam_data = 'ALL_0';
    case 1
        vicon_table = table({'S1_Head1', 'S1_Head2', 'S1_Head3'},{ 'S1_Torso1', 'S1_Torso2', 'S1_Torso3'},...
            {'S1_Left_Arm1'}, {'S1_Left_Arm2'},{ 'S1_Left_Arm3'},...
            {'S1_Right_Arm1'},{ 'S1_Right_Arm2'},{ 'S1_Right_Arm3'},...
            {'S1_Right_Leg1', 'S1_Right_Leg2'},{ 'S1_Right_Leg3'},...
            {'S1_Left_Leg1', 'S1_Left_Leg2'},{ 'S1_Left_Leg3'},...
            'VariableNames',{'Head','Torso','Left_Shoulder','Left_Elbow','Left_Wrist',...
            'Right_Shoulder','Right_Elbow','Right_Wrist','Right_Knee','Right_Feet','Left_Knee','Left_Feet'});
        video_cam_timestamp = 'timestamps_1';
        video_cam_data = 'ALL_1';
end



video_table = table({'Nose','REye','LEye','REar','LEar'},...
    {'Neck','MidHip','RHip','LHip'},...
    {'LShoulder'}, {'LElbow'}, {'LWrist'},{ 'LKnee'},{ 'LAnkle','LBigToe','LSmallToe','LHeel'},...
    {'RShoulder'},{ 'RElbow'},{ 'RWrist'},{'RKnee'},{ 'RAnkle','RBigToe','RSmallToe', 'RHeel'},...
    'VariableNames',{'Head','Torso','Left_Shoulder','Left_Elbow','Left_Wrist',...
    'Left_Knee','Left_Feet','Right_Shoulder','Right_Elbow','Right_Wrist','Right_Knee','Right_Feet'});


vicon = IBS_get_sub_behavior_data('body_landmarks_vicon',Dyad_no,Sub_no,condition,raw_data_dir,varargin_table);
%%
normalizing_dim = 1; % for each column
average_dim = 2;% across the columns
func = @(x) mean(normalize(table2array(vicon(:,contains(vicon.Properties.VariableNames,vicon_table.(x)))),normalizing_dim,'zscore'),average_dim,'omitnan');

averaged_vicon = cellfun(@(x) func(x),vicon_table.Properties.VariableNames,'UniformOutput',0);

averaged_vicon = array2table(cat(2,averaged_vicon{:}),'VariableNames',vicon_table.Properties.VariableNames);
averaged_vicon.timepoints = vicon.timepoints;

%%
video_cam_data = cellfun(@(x) IBS_get_video_cam_all_body_parts(Dyad_no,Sub_no,...
    condition,raw_data_dir,varargin_table,video_table,x),body_parts,'UniformOutput',0);
video_cam_data = cat(2,video_cam_data{:});

video_cam_data = [video_cam_data(:,1:2:end) video_cam_data(:,2)];
video_cam_data = array2table(video_cam_data,'VariableNames',[body_parts 'timepoints']);

%%
[~,~,vicon_movement] = IBS_compute_body_movement_vicon(Dyad_no,Sub_no,condition);

diff_fun = @(x) x-min(x);
vicon_movement = varfun(diff_fun,vicon_movement);
%%
mov_sizes = @(x) {[min(x) quantile(x,1/3)],[quantile(x,1/3) quantile(x,2*(1/3))],...
    [quantile(x,2*(1/3)) max(x)]};  


savemat_fname = ['E:\\Projects\\IBS\\Results\\Vicon\\vicon_video_corresp\\global_range_vicon.mat'];
load(savemat_fname,'all_subs')
all_subs = all_subs(:,1:12);
mov_size_data = varfun(mov_sizes,all_subs);

%%
vicon_interp_fun = @(x) interp1(averaged_vicon.timepoints,x,0:0.1:120);
video_interp_fun = @(x) interp1(video_cam_data.timepoints,x,0:0.1:120);

% varfun(vicon_interp_fun,averaged_vicon)
averaged_vicon_interp = arrayfun(@(x) vicon_interp_fun(table2array(averaged_vicon(:,x))),...
    1:(size(averaged_vicon,2)-1),'UniformOutput',0);
averaged_vicon_interp = cat(1,averaged_vicon_interp{:})';

averaged_vicon_mov_interp = arrayfun(@(x) vicon_interp_fun(table2array(vicon_movement(:,x))),...
    1:(size(vicon_movement,2)-1),'UniformOutput',0);
averaged_vicon_mov_interp = cat(1,averaged_vicon_mov_interp{:})';

averaged_video_interp = arrayfun(@(x) video_interp_fun(table2array(video_cam_data(:,x))),...
    1:(size(video_cam_data,2)-1),'UniformOutput',0);
averaged_video_interp = cat(1,averaged_video_interp{:})';

averaged_vicon_interp= array2table(averaged_vicon_interp,'VariableNames',body_parts);
averaged_vicon_mov_interp= array2table(averaged_vicon_mov_interp,'VariableNames',body_parts);
averaged_video_interp= array2table(averaged_video_interp,'VariableNames',body_parts);
%%

% 
corr_small_med_large = cellfun(@(x) cellfun(@(y) IBS_compute_vicon_video_piecewise_corr(averaged_vicon_interp,...
    averaged_video_interp,averaged_vicon_mov_interp,y,x),mov_size_data.(x),'UniformOutput',0),...
    mov_size_data.Properties.VariableNames,'UniformOutput',0);


corr_small_med_large = cell2mat(cat(1,corr_small_med_large{:}));
corr_small_med_large = corr_small_med_large';


end



function video_cam_data = IBS_get_video_cam_all_body_parts(Dyad_no,Sub_no,condition,raw_data_dir,varargin_table,video_table,body_part)
varargin_table.labels = video_table.(body_part);

video_cam_data = IBS_get_sub_behavior_data('video_openpose_landmarks_manual_cleaned',Dyad_no,Sub_no,condition,raw_data_dir,varargin_table);

video_cam_data = table2array(video_cam_data);
end


function correlation_coeff = IBS_compute_vicon_video_piecewise_corr(vicon_data,video_data,mov_data,data_split,body_partname)
if ~sum(any(isnan(data_split)))
    body_part = strrep(body_partname,'Fun_','');
    cur_vicon = vicon_data.(body_part);
    cur_video = video_data.(body_part);
    mov_data = mov_data.(body_part);
    low_lim_indices = find(mov_data>=data_split(1));
    upper_lim_indices = find(mov_data<=data_split(2));
       
    cur_vicon = cur_vicon(intersect(low_lim_indices,upper_lim_indices));
    cur_video = cur_video(intersect(low_lim_indices,upper_lim_indices));
    
    correlation_coeff = nancorr(cur_vicon',cur_video);
else
    correlation_coeff = NaN;
end

end

