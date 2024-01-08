function [gof_cell] = IBS_vicon_video_mov_size(results_corr,results_mag_mov,conditions)
%IBS_VICON_VIDEO_MOV_SIZE evaluates the relationship of the profile of vicon-video depending of the size of movement performed
%
% SYNOPSIS: IBS_vicon_video_mov_size
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

if nargin <3
    conditions = {'FaNoOcc_1','FaNoOcc_2','FaNoOcc_3','FaOcc_1','FaOcc_2','FaOcc_3'};
    
end

subset = @(x) x{1};

cond_name = cellfun(@(x) subset(strsplit(x,'_')),conditions,'UniformOutput',0);
cond_name = strjoin(unique(cond_name),'_');


marker_names = {'Head','Torso','Left_Shoulder','Left_Elbow','Left_Wrist',...
    'Left_Knee','Left_Feet','Right_Shoulder','Right_Elbow','Right_Wrist','Right_Knee','Right_Feet'};
gof_cell = cell(1,length(marker_names));


savefname_mat = ['E:\\Projects\\IBS\\Results\\Vicon\\vicon_video_corresp\\figures\\relationship_mov_corr_' cond_name '.mat'];


if exist(savefname_mat,'file')
    load(savefname_mat)
else
    for plot_no = 1:length(marker_names)
        a = results_corr(plot_no,:);
        b = results_mag_mov(plot_no,:);
        isFinite = ~isnan(a) & ~isnan(b);
        a = a(isFinite);
        b = b(isFinite);
        %     b(b<35) = nan;
        [fit, gof] = createFits(b,a);
         gof_cell{plot_no} = gof;
         
         %% 13-06-23 added this below:
         
         
         
    end
    save(savefname_mat,'results_corr','results_mag_mov','cond_name','gof_cell')
end


end