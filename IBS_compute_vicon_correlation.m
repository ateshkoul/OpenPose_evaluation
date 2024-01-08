function correlation_coeff = IBS_compute_vicon_correlation(Dyad_no,condition)
%IBS_compute_vicon_correlation
%
% SYNOPSIS: IBS_compute_vicon_correlation
%
% INPUT function to check the correlation for vicon movements
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
varargin_table = table();
varargin_table.behav_analysis = 'joint';
varargin_table.analysis_sub_type = '_insta_abs_detrend';


vicon_table_S0 = table({'S0_Head1','S0_Head2','S0_Head3','S0_Head4'},...
    {'S0_Torso1', 'S0_Torso2', 'S0_Torso3'},...
    {'S0_Left_Arm1'}, {'S0_Left_Arm2'}, {'S0_Left_Arm3'},{ 'S0_Left_Leg1', 'S0_Left_Leg2'},{ 'S0_Left_Leg3'},...
    {'S0_Right_Arm1'},{ 'S0_Right_Arm2'},{ 'S0_Right_Arm3'},{'S0_Right_Leg1', 'S0_Right_Leg2'},{ 'S0_Right_Leg3'},...
    'VariableNames',{'Head','Torso','Left_Shoulder','Left_Elbow','Left_Wrist',...
    'Left_Knee','Left_Feet','Right_Shoulder','Right_Elbow','Right_Wrist','Right_Knee','Right_Feet'});

vicon_table_S1 = table({'S1_Head1', 'S1_Head2', 'S1_Head3'},{ 'S1_Torso1', 'S1_Torso2', 'S1_Torso3'},...
    {'S1_Left_Arm1'}, {'S1_Left_Arm2'},{ 'S1_Left_Arm3'},...
    {'S1_Right_Arm1'},{ 'S1_Right_Arm2'},{ 'S1_Right_Arm3'},...
    {'S1_Right_Leg1', 'S1_Right_Leg2'},{ 'S1_Right_Leg3'},...
    {'S1_Left_Leg1', 'S1_Left_Leg2'},{ 'S1_Left_Leg3'},...
    'VariableNames',{'Head','Torso','Left_Shoulder','Left_Elbow','Left_Wrist',...
    'Right_Shoulder','Right_Elbow','Right_Wrist','Right_Knee','Right_Feet','Left_Knee','Left_Feet'});

marker_names = {'Head','Torso','Left_Shoulder','Left_Elbow','Left_Wrist',...
    'Right_Shoulder','Right_Elbow','Right_Wrist','Right_Knee','Right_Feet','Left_Knee','Left_Feet'};

%%
varargin_table.output_type = 'change';

vicon_S0 = IBS_get_sub_behavior_data('body_landmarks_vicon',Dyad_no,0,condition,raw_data_dir,varargin_table);
vicon_S1 = IBS_get_sub_behavior_data('body_landmarks_vicon',Dyad_no,1,condition,raw_data_dir,varargin_table);

%%
% func = @(x) mean(table2array(vicon(:,contains(vicon.Properties.VariableNames,vicon_table.(x)))),2,'omitnan');
normalizing_dim = 1; % for each column
average_dim = 2;% across the columns

%%
func_S0 = @(x) mean(normalize(table2array(vicon_S0(:,contains(vicon_S0.Properties.VariableNames,vicon_table_S0.(x)))),...
    normalizing_dim,'zscore'),average_dim,'omitnan');

averaged_vicon_S0 = cellfun(@(x) func_S0(x),vicon_table_S0.Properties.VariableNames,'UniformOutput',0);

averaged_vicon_S0 = array2table(cat(2,averaged_vicon_S0{:}),'VariableNames',vicon_table_S0.Properties.VariableNames);
averaged_vicon_S0.timepoints = vicon_S0.timepoints;

%%
func_S1 = @(x) mean(normalize(table2array(vicon_S1(:,contains(vicon_S1.Properties.VariableNames,vicon_table_S1.(x)))),...
    normalizing_dim,'zscore'),average_dim,'omitnan');

averaged_vicon_S1 = cellfun(@(x) func_S1(x),vicon_table_S1.Properties.VariableNames,'UniformOutput',0);

averaged_vicon_S1 = array2table(cat(2,averaged_vicon_S1{:}),'VariableNames',vicon_table_S1.Properties.VariableNames);
averaged_vicon_S1.timepoints = vicon_S1.timepoints;


correlation_coeff = cell2mat(cellfun(@(x) nancorr(averaged_vicon_S1.(x),averaged_vicon_S0.(x)),marker_names,'UniformOutput',0));

end
