clear

body_parts = {'Head','Torso','Left_Shoulder','Left_Elbow','Left_Wrist',...
    'Left_Knee','Left_Feet','Right_Shoulder','Right_Elbow','Right_Wrist','Right_Knee','Right_Feet'};
conditions_FaNoOcc = {'FaNoOcc_1','FaNoOcc_2','FaNoOcc_3'};
correlations_vicon_FaNoOcc = arrayfun(@(x) cellfun(@(cond) IBS_compute_vicon_correlation(x,cond),...
    conditions_FaNoOcc,'UniformOutput',0),1:23,'UniformOutput',0);


correlations_avg_FaNoOcc = cellfun(@(x) num2cell(mean(cat(1,x{:}),'omitnan')),correlations_vicon_FaNoOcc,'UniformOutput',0);
correlations_avg_FaNoOcc = cell2mat(cat(1,correlations_avg_FaNoOcc{:}));


conditions_FaOcc = {'FaOcc_1','FaOcc_2','FaOcc_3'};
correlations_vicon_FaOcc = arrayfun(@(x) cellfun(@(cond) IBS_compute_vicon_correlation(x,cond),...
    conditions_FaOcc,'UniformOutput',0),1:23,'UniformOutput',0);


correlations_avg_FaOcc = cellfun(@(x) num2cell(mean(cat(1,x{:}),'omitnan')),correlations_vicon_FaOcc,'UniformOutput',0);
correlations_avg_FaOcc = cell2mat(cat(1,correlations_avg_FaOcc{:}));

savemat_fname = ['E:\\Projects\\IBS\\Results\\Vicon\\vicon_video_corresp\\results_vicon_correlation_all.mat'];

correlations_avg_FaOcc_table = array2table(correlations_avg_FaOcc,'VariableNames',body_parts);
correlations_avg_FaNoOcc_table = array2table(correlations_avg_FaNoOcc,'VariableNames',body_parts);
save(savemat_fname,'correlations_avg_FaOcc_table','correlations_avg_FaNoOcc_table','body_parts','conditions_FaOcc','conditions_FaNoOcc')

savecsv_fname = ['E:\\Projects\\IBS\\Results\\Vicon\\vicon_video_corresp\\results_corr_vicon_FaNoOcc.csv'];
writetable(correlations_avg_FaNoOcc_table,savecsv_fname)
