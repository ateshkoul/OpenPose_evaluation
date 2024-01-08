function IBS_mov_based_corresp_video_vicon
%IBS_MOV_BASED_CORRESP_VIDEO_VICON
%
% SYNOPSIS: IBS_mov_based_corresp_video_vicon
%
% INPUT function to compute the correspondence between vicon and video for each movement type on a trial to trial basis.
%
% OUTPUT
%
% REMARKS
%
% created with MATLAB ver.: 9.8.0.1359463 (R2020a) Update 1 on Microsoft Windows 10 Pro Version 10.0 (Build 19042)
%
% created by: Atesh
% DATE: 17-Jan-2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data_split = 'median_split';
data_split = 'three_splits';
analysis_type = 'absolute';%'relative'
savemat_fname = ['E:\\Projects\\IBS\\Results\\Vicon\\vicon_video_corresp\\' analysis_type '_' data_split '_mov_based_corresp_video_vicon.mat'];
% savemat_fname = ['E:\\Projects\\IBS\\Results\\Vicon\\vicon_video_corresp\\relative_mov_based_corresp_video_vicon.mat'];

if exist(savemat_fname,'file')
    load(savemat_fname,'results_avg','results_corr','body_parts','conditions');
else
    body_parts = {'Head','Torso','Left_Shoulder','Left_Elbow','Left_Wrist',...
        'Left_Knee','Left_Feet','Right_Shoulder','Right_Elbow','Right_Wrist','Right_Knee','Right_Feet'};
    Dyads = 1:23;
    
    % body_parts = {'Head'};
    conditions = {'FaNoOcc_1','FaNoOcc_2','FaNoOcc_3','FaOcc_1','FaOcc_2','FaOcc_3'};
    
    results_corr = arrayfun(@(Dyads) arrayfun(@(sub) cellfun(@(cond) IBS_compute_mov_based_corresp_video_vicon(Dyads,sub,cond),...
        conditions,'UniformOutput',0),0:1,'UniformOutput',0),1:23,'UniformOutput',0);
    
    
    
    results_corr = cellfun(@(Dyads) cellfun(@(sub) mean(cat(3,sub{:}),3,'omitnan'),...
        Dyads,'UniformOutput',0),results_corr,'UniformOutput',0);
    results_corr = cat(2,results_corr{:});
    results_corr = cat(3,results_corr{:});
    
    
    % results_corr_sizewise = arrayfun(@(x) squeeze(results_corr(x,:,:)),1:3,'UniformOutput',0);
    % cellfun(@(x) bar(mean(x,2,'omitnan')),results_corr,'UniformOutput',0);
    %
    % for mov_size = 1:3
    %     subplot(1,3,mov_size)
    %     bar(mean(results_corr{mov_size},2,'omitnan'))
    %     ylim([-0.01 0.5])
    % end
    
    results_corr = arrayfun(@(x) squeeze(results_corr(:,x,:)),1:12,'UniformOutput',0);
    
    results_avg = cell2mat(cellfun(@(x) mean(x,2,'omitnan'),results_corr,'UniformOutput',0))';
    %%
    results_corr_sizewise = cat(3,results_corr{:});
    results_corr_sizewise_perm = permute(results_corr_sizewise,[3 2 1]);
    results_corr_table = array2table(reshape(results_corr_sizewise_perm,length(body_parts),length(Dyads)*2*3)',...
        'VariableNames',body_parts);
    
    mov_sizes = repmat({'low','med','high'},length(Dyads)*2,1);
    results_corr_table.Mov_size = mov_sizes(:);
    save(savemat_fname,'results_avg','results_corr','body_parts','conditions','results_corr_table')
    savecsv_fname = ['E:\\Projects\\IBS\\Results\\Vicon\\vicon_video_corresp\\' analysis_type '_' data_split '_mov_based_corresp_video_vicon.csv'];
    writetable(results_corr_table,savecsv_fname)
    
end
%%

%%
figure('units','normalized','outerposition',[0 0 1 1]);
b = bar(results_avg);
ylim([-0.01 1])
results_corr_sizewise = cat(3,results_corr{:});

n_sub = 23;
n_bodyparts = 12;
x = repmat(1:n_bodyparts,[n_sub 1]);
for plot= 1:length(b)
    hold on
    cur_b = b(:,plot);
    Xend_points = cur_b.XEndPoints;
    x = repmat(Xend_points,46,1);
    arrayfun(@(t) scatter(x(:,t),squeeze(results_corr_sizewise(plot,:,t))),1:size(Xend_points,2),'UniformOutput',0)
end

xticks(1:12)
xticklabels(body_parts)
set(gca,'TickLabelInterpreter','none')
savefname = ['E:\\Projects\\IBS\\Results\\Vicon\\vicon_video_corresp\\figures\\' analysis_type '_' data_split '_mov_based_corresp_video_vicon.tif'];

% savefname = ['E:\\Projects\\IBS\\Results\\Vicon\\vicon_video_corresp\\figures\\relative_mov_based_corresp_video_vicon.tif'];
saveas(gcf,savefname)
end