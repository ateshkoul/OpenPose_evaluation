function [results_corr_avg_cond_body,results_corr_avg_cond_body_table,results_corr] = IBS_vicon_video_corr(conditions,varargin_table)
%IBS_VICON_VIDEO_CORR computes the correlation between the vicon and video data
%
% SYNOPSIS: IBS_vicon_video_corr
%
% INPUT
%
% OUTPUT
%
% REMARKS
%
% created with MATLAB ver.: 9.8.0.1359463 (R2020a) Update 1 on Microsoft Windows 10 Pro Version 10.0 (Build 19042)
%
% created by: Atesh
% DATE: 09-May-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if nargin <1
    conditions = {'FaNoOcc_1','FaNoOcc_2','FaNoOcc_3','FaOcc_1','FaOcc_2','FaOcc_3'};
    varargin_table = table();
end
if nargin <2
    varargin_table = table();
end
subset = @(x) x{1};


%%

if ~ismember('behav_analysis',varargin_table.Properties.VariableNames); varargin_table.behav_analysis = 'joint';end
if ~ismember('analysis_sub_type',varargin_table.Properties.VariableNames); varargin_table.analysis_sub_type = '_insta_abs_detrend';end
if ~ismember('behavior',varargin_table.Properties.VariableNames);varargin_table.behavior ='video_openpose_landmarks_manual_cleaned';end


%%
cond_name = cellfun(@(x) subset(strsplit(x,'_')),conditions,'UniformOutput',0);
cond_name = strjoin(unique(cond_name),'_');

try varargin_table.analysis_type
    cond_name = [varargin_table.analysis_type cond_name];
catch
end

body_parts = {'Head','Torso','Left_Shoulder','Left_Elbow','Left_Wrist',...
    'Left_Knee','Left_Feet','Right_Shoulder','Right_Elbow','Right_Wrist','Right_Knee','Right_Feet'};

savefname_mat = ['E:\\Projects\\IBS\\Results\\Vicon\\vicon_video_corresp\\figures\\vicon_video_corresp_' cond_name ' .mat'];



% body_parts = {'Head'};
if exist(savefname_mat,'file')
    
    load(savefname_mat,'results_corr','results_corr_avg','results_corr_avg_cond',...
        'results_corr_avg_cond_body','cond_name')
    
else
    % this is over conditions, subjects, dyads and then body parts as we
    % are not computing anything related to the dyad
    results_corr = cellfun(@(body_part) arrayfun(@(Dyads) arrayfun(@(sub) cellfun(@(cond) IBS_compute_corresp_vicon_video(Dyads,sub,cond,body_part,varargin_table),...
        conditions,'UniformOutput',0),0:1,'UniformOutput',0),1:23,'UniformOutput',0),body_parts,'UniformOutput',0);
    
    results_corr_avg = cellfun(@(corr_dyad) cellfun(@(corr_sub) cellfun(@(corr_cond) mat2cell(mean(cell2mat(corr_cond),'omitnan'),1,1),...
        corr_sub,'UniformOutput',0),corr_dyad,'UniformOutput',0),results_corr,'UniformOutput',0);
    
    results_corr_avg_cond = cellfun(@(body_part) cellfun(@(cond) cat(2,cond{:}),...
        body_part,'UniformOutput',0),results_corr_avg,'UniformOutput',0);
    
    results_corr_avg_cond_body = cellfun(@(body_part) cell2mat(cat(2,body_part{:})),...
        results_corr_avg_cond,'UniformOutput',0);
    results_corr_avg_cond_body = cat(1,results_corr_avg_cond_body{:});
    
    
    save(savefname_mat,'results_corr','results_corr_avg','results_corr_avg_cond',...
        'results_corr_avg_cond_body','cond_name')
end

savefname = ['E:\\Projects\\IBS\\Results\\Vicon\\vicon_video_corresp\\figures\\vicon_video_corresp_all_' cond_name '.tif'];
savefname = ['E:\\Projects\\IBS\\Results\\Vicon\\vicon_video_corresp\\figures\\25_hz_vicon_video_corresp_all_' cond_name '.tif'];

if exist(savefname,'file')
else
    % plot
    n_sub = 46;
    n_bodyparts = 12;
    x = repmat(1:n_bodyparts,[n_sub 1])';
    figure('units','normalized','outerposition',[0 0 1 1]);
    bar(mean(results_corr_avg_cond_body,2,'omitnan'))
    hold on
    arrayfun(@(t) scatter(x(t,:)+normrnd(0,1,[1 n_sub])/10,results_corr_avg_cond_body(t,:)),1:size(x,1),'UniformOutput',0)
    xticklabels(body_parts)
    set(gca,'TickLabelInterpreter','none')
    title('Correlation between vicon and video')
    ylabel('correlation value')
%     saveas(gcf,savefname)
end
end
