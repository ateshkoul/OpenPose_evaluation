function correlation_coeff = IBS_compute_corresp_vicon_video(Dyad_no,Sub_no,condition,body_part,varargin_table)
%IBS_COMPUTE_CORRESP_VICON_VIDEO
%
% SYNOPSIS: IBS_compute_corresp_vicon_video
%
% INPUT function to check the correspondence between the vicon and video based landmark detection
%
% OUTPUT
%
% REMARKS
%
% created with MATLAB ver.: 9.8.0.1359463 (R2020a) Update 1 on Microsoft Windows 10 Pro Version 10.0 (Build 19042)
%
% created by: Atesh
% DATE: 10-Dec-2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_data_dir = 'D:\\iannettilab_dropbox\\Dropbox\\Koul_Atesh\\IBS\\';

if nargin < 5
    varargin_table = table();
    varargin_table.behav_analysis = 'joint';
    varargin_table.analysis_sub_type = '_insta_abs_detrend';
    varargin_table.behavior ='video_openpose_landmarks_manual_cleaned'; 
    
    
    
end


if ~ismember('corr_fun',varargin_table.Properties.VariableNames); varargin_table.corr_fun = @nancorr;end


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


varargin_table.labels = video_table.(body_part);


vicon = IBS_get_sub_behavior_data('body_landmarks_vicon',Dyad_no,Sub_no,condition,raw_data_dir,varargin_table);


video_cam=IBS_get_sub_behavior_data(varargin_table.behavior,Dyad_no,Sub_no,condition,raw_data_dir,varargin_table);
%%

% func = @(x) mean(table2array(vicon(:,contains(vicon.Properties.VariableNames,vicon_table.(x)))),2,'omitnan');
normalizing_dim = 1; % for each column
average_dim = 2;% across the columns


switch(varargin_table.analysis_type)
    case 'ICC_'
        func = @(x) mean(table2array(vicon(:,contains(vicon.Properties.VariableNames,vicon_table.(x)))),average_dim,'omitnan');
    otherwise
        func = @(x) mean(normalize(table2array(vicon(:,contains(vicon.Properties.VariableNames,vicon_table.(x)))),normalizing_dim,'zscore'),average_dim,'omitnan');
end
averaged_vicon = cellfun(@(x) func(x),vicon_table.Properties.VariableNames,'UniformOutput',0);

averaged_vicon = array2table(cat(2,averaged_vicon{:}),'VariableNames',vicon_table.Properties.VariableNames);
averaged_vicon.timepoints = vicon.timepoints;

body_data_video_cam = interp1(video_cam.(video_cam_timestamp),video_cam.(video_cam_data),0:0.1:120);
averaged_body_data_vicon = interp1(averaged_vicon.timepoints,averaged_vicon.(body_part),0:0.1:120);

nan_values = find(isnan(body_data_video_cam));

if ~(length(nan_values)==length(body_data_video_cam))
    
    switch(varargin_table.analysis_type)
        case 'ICC_'
            averaged_body_data_vicon(nan_values) = [];
            body_data_video_cam(nan_values) = [];            
            body_data_video_cam = body_data_video_cam+(averaged_body_data_vicon(1,1)-body_data_video_cam(1,1));
    end
    correlation_coeff = varargin_table.corr_fun(averaged_body_data_vicon,body_data_video_cam);
    
else
    correlation_coeff = NaN;
    
end

end
