function [behav_data] = IBS_get_sub_behavior_data(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)

% raw_data_dir = 'D:\\iannettilab_dropbox\\Dropbox\\Koul_Atesh\\IBS\\';
% varargin_table = table();
% varargin_table.behav_analysis = 'joint';
% varargin_table.analysis_sub_type = '_insta_abs_detrend';
% % IBS_get_sub_behavior_data('Gaze_nose_dist',1,1,'FaNoOcc_1',raw_data_dir,varargin_table)
% s = IBS_get_sub_behavior_data('body_landmarks_vicon',10,1,'FaNoOcc_1',raw_data_dir,varargin_table);
% behavior_params = IBS_get_behavior_type(behavior,varargin_table.behav_analysis);
behavior_params = IBS_get_behavior_type(behavior,varargin_table.behav_analysis,varargin_table.analysis_sub_type);

norm_fun = behavior_params.norm_fun{1,1};
%
% e.g. import_fun = @IBS_import_pupil_csv_data
% e.g. import_fun = @IBS_import_gaze_positions_on_body_data
% varargin_table.norm_method = @(x) normalize(x);
varargin_table.norm_method = 'zscore';
switch(behavior)
    case 'Eye_tracker_pupil'
        behav_data = IBS_load_Eye_tracker_pupil(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
        %     case 'Eye_tracker_gaze_positions_on_body'
    case 'Eye_tracker_gaze'
        behav_data = IBS_load_Eye_tracker_gaze(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
    case 'Gaze_nose_dist'
        behav_data = IBS_load_Eye_tracker_gaze_nose_dist(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
    case 'Eye_tracker_blink'
        behav_data = IBS_load_Eye_tracker_blink(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
    case 'Video_landmarks'
        behav_data = IBS_load_Video_landmarks(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
    case 'manipulated_Video_landmarks'
        behav_data = IBS_load_manipulated_Video_landmarks(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
    case 'Video_manual_labelled'
        behav_data = IBS_load_Video_manual_labelled(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
    case 'ECG'
        behav_data = IBS_load_ECG(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
        
        
    case 'norm_ECG'
        behav_data = IBS_load_norm_ECG(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
    case 'EDA'
        behav_data =IBS_load_EDA(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
        
    case 'Resp'
        behav_data =IBS_load_resp(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
        
    case 'Smile_auto'
        behav_data = IBS_load_smile_auto(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
    case 'video_openpose_landmarks_manual_cleaned'
        behav_data = IBS_load_video_openpose_landmarks_manual_cleaned(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
        
    case 'video_openpose_landmarks_manual_cleaned_in_vicon_units'
        behav_data = IBS_load_video_openpose_landmarks_manual_cleaned_in_vicon_units(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
   
    case 'video_openpose_landmarks_manual_cleaned_add_body'
        behav_data = IBS_load_video_openpose_landmarks_manual_cleaned_add_body(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
    case 'manipulated_video_openpose_landmarks_manual_cleaned'
        behav_data = IBS_load_manipulated_video_openpose_landmarks_manual_cleaned(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
        
    case 'body_landmarks_vicon'
        behav_data = IBS_load_body_landmarks_vicon(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
    case 'vicon_video_landmark'
        behav_data = IBS_load_vicon_video_landmark(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
        
    case 'Abs_video_landmarks'
        behav_data = IBS_load_abs_video_openpose_landmarks(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table);
        
end



end




function behav_data = IBS_load_Eye_tracker_pupil(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)

import_fun = @IBS_import_pupil_csv_data;
varargin_table.cols_to_keep = {'pupil_timestamp','diameter'};
sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];
behav_data = import_fun([sub_dir '\\pupil_positions.csv']);
behav_data = behav_data(:,ismember(behav_data.Properties.VariableNames,varargin_table.cols_to_keep));
behav_data = rmmissing(behav_data);
%         behav_data.pupil_timestamp = behav_data.pupil_timestamp - behav_data.pupil_timestamp(1,1);

behav_data_starttime = IBS_import_eye_tracker_blink_info([sub_dir '\\export_info.csv']);


%         behav_data.time_stamps = behav_data.time_stamps - behav_data.time_stamps(1,1);
if Dyad_no == 7 && Sub == 1 && strcmp(condition,'FaNoOcc_3')
    % there is a delay of 28.4s as the acquisition started at that point
    behav_data.pupil_timestamp = ((behav_data.pupil_timestamp - behav_data_starttime)+ 28.4) + IBS_get_behav_delay(Dyad_no,Sub,condition);
    
elseif contains(condition,{'Task_1','Task_2','Task_3'})
    
    behav_data.pupil_timestamp = (behav_data.pupil_timestamp - behav_data_starttime) + IBS_get_behav_delay_task(Dyad_no,Sub,condition);
else
    behav_data.pupil_timestamp = (behav_data.pupil_timestamp - behav_data_starttime)+ IBS_get_behav_delay(Dyad_no,Sub,condition);
    
    
    
end


%         behav_data.pupil_timestamp = (behav_data.pupil_timestamp - behav_data.pupil_timestamp(1,1)) + IBS_get_behav_delay(Dyad_no,Sub,condition);
% behav_data.pupil_timestamp = (behav_data.pupil_timestamp - behav_data_starttime) + IBS_get_behav_delay(Dyad_no,Sub,condition);

behav_data(behav_data.diameter==0,:) = [];
behav_data(behav_data.diameter<0,:) = [];
behav_data(abs(normalize(behav_data.diameter))>2.5,:) = [];
behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);


behav_data.Properties.VariableNames = strrep(behav_data.Properties.VariableNames,'pupil_timestamp','time_stamps');

% keep as a table
%         behav_data = table2array(behav_data);
%         behav_data(:,1) = behav_data(:,1)-behav_data(1,1);
% may be it's better to normalize them just before the glm
%         behav_data(:,2:end) = norm_fun(behav_data(:,2:end),varargin_table.norm_method);
% this is because in python the script is nFrames -1 so the last frame is
% skipped
behav_data = behav_data(1:end-1,:);

%% outlier removal
if ismember('outlier_removal',varargin_table.Properties.VariableNames)    
    dep_measure = [varargin_table.cols_to_keep{1,2} '_'];
    behav_data{abs(behav_data.([dep_measure num2str(Sub)]))>(mean(behav_data.([dep_measure num2str(Sub)])) + ...
        varargin_table.outlier_removal*std(behav_data.([dep_measure num2str(Sub)]))),[dep_measure num2str(Sub)]} = NaN;
end

end




function behav_data = IBS_load_Eye_tracker_gaze(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)
% sub is 0 and 1
% this doesn't work for a table because table.import_fun(x) where x
% is an argument to function is taken as an index for a table
%         varargin_table.import_fun = @IBS_import_pupil_csv_data;
import_fun = @IBS_import_gaze_positions_on_body_data;
% eye_on_body for a lot of subjects is always 1
%         varargin_table.cols_to_keep = {'time_stamps','eye_on_body','eye_on_face','eye_on_face_ellipse'};
%         varargin_table.cols_to_keep = {'time_stamps','eye_on_face','eye_on_face_ellipse','eye_on_eye_ellipse'};
varargin_table.cols_to_keep = {'time_stamps','eye_on_face'};
% check for compatibility with iscience
% varargin_table.cols_to_keep = {'time_stamps','eye_on_face_ellipse'};

sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];
behav_data = import_fun([sub_dir '\\results_ellipse_eye_gaze_api_full.csv']);
behav_data = behav_data(:,ismember(behav_data.Properties.VariableNames,varargin_table.cols_to_keep));
behav_data = rmmissing(behav_data);

behav_data_starttime = IBS_import_eye_tracker_blink_info([sub_dir '\\export_info.csv']);


%         behav_data.time_stamps = behav_data.time_stamps - behav_data.time_stamps(1,1);
if Dyad_no == 7 && Sub == 1 && strcmp(condition,'FaNoOcc_3')
    % there is a delay of 28.4s as the acquisition started at that point
    behav_data.time_stamps = ((behav_data.time_stamps - behav_data_starttime)+ 28.4) + IBS_get_behav_delay(Dyad_no,Sub,condition);
    
elseif Dyad_no == 11 && Sub == 1 && strcmp(condition,'Task_1')
    behav_data.time_stamps = ((behav_data.time_stamps - behav_data_starttime)+ 39.668) + IBS_get_behav_delay(Dyad_no,Sub,condition);
    
else
    behav_data.time_stamps = (behav_data.time_stamps - behav_data_starttime)+ IBS_get_behav_delay(Dyad_no,Sub,condition);
end


%         behav_data.time_stamps = behav_data.time_stamps + IBS_get_behav_delay(Dyad_no,Sub,condition);

% this is useful when merging columns across behavioral tasks
behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);
% to reduce the name in lme

behav_data.Properties.VariableNames = strrep(behav_data.Properties.VariableNames,'_on_','');
behav_data.Properties.VariableNames = strrep(behav_data.Properties.VariableNames,'_ellipse','_e');

% this helps in keeping them unique
%         behav_data = table2array(behav_data);
%         behav_data(:,1) = behav_data(:,1)-behav_data(1,1);
% may be it's better to normalize them just before the glm
%         behav_data(:,2:end) = norm_fun(behav_data(:,2:end),varargin_table.norm_method);

end

function behav_data = IBS_load_Eye_tracker_gaze_nose_dist(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)

% this doesn't work for a table because table.import_fun(x) where x
% is an argument to function is taken as an index for a table
%         varargin_table.import_fun = @IBS_import_pupil_csv_data;
import_fun = @IBS_import_gaze_nose_dist;
% eye_on_body for a lot of subjects is always 1
%         varargin_table.cols_to_keep = {'time_stamps','eye_on_body','eye_on_face','eye_on_face_ellipse'};
%         varargin_table.cols_to_keep = {'time_stamps','eye_on_face','eye_on_face_ellipse','eye_on_eye_ellipse'};
varargin_table.cols_to_keep = {'time_stamps','eye_gaze_distance'};

sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];
behav_data = import_fun([sub_dir '\\results_ellipse_eye_gaze_api_full_gaze_nose_distance.csv']);
behav_data = behav_data(:,ismember(behav_data.Properties.VariableNames,varargin_table.cols_to_keep));
behav_data = rmmissing(behav_data);

behav_data_starttime = IBS_import_eye_tracker_blink_info([sub_dir '\\export_info.csv']);


%         behav_data.time_stamps = behav_data.time_stamps - behav_data.time_stamps(1,1);
if Dyad_no == 7 && Sub == 1 && strcmp(condition,'FaNoOcc_3')
    % there is a delay of 28.4s as the acquisition started at that point
    behav_data.time_stamps = ((behav_data.time_stamps - behav_data_starttime)+ 28.4) + IBS_get_behav_delay(Dyad_no,Sub,condition);
    
elseif contains(condition,{'Task_1','Task_2','Task_3'})
    
    behav_data.time_stamps = (behav_data.time_stamps - behav_data_starttime) + IBS_get_behav_delay_task(Dyad_no,Sub,condition);
else
    behav_data.time_stamps = (behav_data.time_stamps - behav_data_starttime)+ IBS_get_behav_delay(Dyad_no,Sub,condition);
    
    
    
end


%         behav_data.time_stamps = behav_data.time_stamps + IBS_get_behav_delay(Dyad_no,Sub,condition);
% inverse the data as lower values mean better social behav
% behav_data.eye_gaze_distance = 1./behav_data.eye_gaze_distance;
% this is because in python the script is nFrames -1 so the last frame is
% skipped
behav_data = behav_data(1:end-1,:);


%%

% %% Filt_2 : Remove outliers (filter too high or low data points)
% [ indices ]             = Giac_StdOverMeanVectorFilter( behav_data.eye_gaze_distance', 3 );   %3 std away from m
% behav_data.eye_gaze_distance(indices)       = nan;
% %
% %
% % %% Filt_3: Interpolate nans (using Dav_interp_nans function)
% %
% behav_data.eye_gaze_distance                     = Dav_interp_nans(behav_data.eye_gaze_distance);
%
% %%
% % behav_data.mouth_size = smoothdata(behav_data.mouth_size,'gaussian',1,'SamplePoints',behav_data.time_stamps_rel);
% standard_time_stamps =linspace(behav_data.time_stamps(1),behav_data.time_stamps(end),length(behav_data.time_stamps));
% behav_data.eye_gaze_distance = interp1(behav_data.time_stamps,behav_data.eye_gaze_distance,standard_time_stamps)';
%
%
% % behav_data.eye_gaze_distance = smoothdata(behav_data.eye_gaze_distance,'lowess',1,'SamplePoints',standard_time_stamps);
%
% % behav_data.eye_gaze_distance = smoothdata(behav_data.eye_gaze_distance,'movmean',1,'SamplePoints',standard_time_stamps);
% % behav_data.eye_gaze_distance = movmean(behav_data.eye_gaze_distance,[0.15 0],'SamplePoints',standard_time_stamps);
% % behav_data.eye_gaze_distance = movmean(behav_data.eye_gaze_distance,0.1,'SamplePoints',standard_time_stamps);
% % behav_data.eye_gaze_distance = movmean(behav_data.eye_gaze_distance,[0 0.2],'SamplePoints',standard_time_stamps);
% % behav_data.eye_gaze_distance = movmean(behav_data.eye_gaze_distance,0.05,'SamplePoints',standard_time_stamps);
% % behav_data.eye_gaze_distance = movmean(behav_data.eye_gaze_distance,[0.3 0],'SamplePoints',standard_time_stamps);
% behav_data.eye_gaze_distance = movmean(behav_data.eye_gaze_distance,[0 0.1],'SamplePoints',standard_time_stamps);

%% lowpass the signal
% standard_time_stamps =linspace(behav_data.time_stamps(1),behav_data.time_stamps(end),length(behav_data.time_stamps));
% % frame_rate = length(behav_data.time_stamps)/(behav_data.time_stamps(end)-behav_data.time_stamps(1));
%
% xt = table();
% xt.eye_gaze_distance = interp1(behav_data.time_stamps,behav_data.eye_gaze_distance,standard_time_stamps)';
% xt.time_stamps = seconds(standard_time_stamps');
% xt = table2timetable(xt,'RowTimes','time_stamps');
% lowpass_hz = 2;
% rmpath('C:\\Users\\Atesh\\Downloads\\fieldtrip-20200327\\external\\signal\\')
% xt = lowpass(xt,lowpass_hz);
%
% addpath('C:\\Users\\Atesh\\Downloads\\fieldtrip-20200327\\external\\signal\\')
% behav_data = timetable2table(xt);
% behav_data.time_stamps = seconds(behav_data.time_stamps);
% clear xt
%%
% this is useful when merging columns across behavioral tasks
behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);
% to reduce the name in lme


% this helps in keeping them unique
%         behav_data = table2array(behav_data);
%         behav_data(:,1) = behav_data(:,1)-behav_data(1,1);
% may be it's better to normalize them just before the glm
%         behav_data(:,2:end) = norm_fun(behav_data(:,2:end),varargin_table.norm_method);


end



function behav_data = IBS_load_Eye_tracker_blink(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)

import_fun = @IBS_import_blink;
varargin_table.cols_to_keep = {'start_timestamp','end_timestamp','duration'};
sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];
behav_data = import_fun([sub_dir '\\blinks.csv']);
behav_data = behav_data(:,ismember(behav_data.Properties.VariableNames,varargin_table.cols_to_keep));
behav_data = rmmissing(behav_data);
%         behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);
% keep as a table
%         behav_data = table2array(behav_data);
%         behav_data(:,1) = behav_data(:,1)-behav_data(1,1);
% may be it's better to normalize them just before the glm
%         behav_data(:,2:end) = norm_fun(behav_data(:,2:end),varargin_table.norm_method);

%         import_fun_info = @(x) jsondecode(fileread(x));
import_fun_info = @IBS_import_eye_tracker_blink_info;


behav_data_starttime = import_fun_info([sub_dir '\\export_info.csv']);
%         import_fun_pupil = @IBS_import_pupil_csv_data;
%         varargin_table.cols_to_keep = {'pupil_timestamp'};
%         sub_dir = IBS_get_sub_behavior_dir('Eye_tracker_pupil',Dyad_no,Sub,condition);
%         sub_dir = [raw_data_dir sub_dir];
%         behav_data_pupil = import_fun_pupil([sub_dir '\\pupil_positions.csv']);
%         behav_data_pupil = behav_data_pupil(:,ismember(behav_data_pupil.Properties.VariableNames,varargin_table.cols_to_keep));
%         behav_data_starttime = behav_data_pupil.pupil_timestamp(1,1);



behav_data.start_timestamp = round(behav_data.start_timestamp - behav_data_starttime,4);
behav_data.end_timestamp = round(behav_data.end_timestamp - behav_data_starttime,4);

end



function behav_data = IBS_load_Video_landmarks(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)


try
    body_part = varargin_table.body_part;
    import_fun = varargin_table.import_fun;
    
catch
    body_part = 'body';
    import_fun = @IBS_import_video_landmarks;
    
end

%         varargin_table.cols_to_keep = {'time_stamps','eye_on_face','eye_on_face_ellipse','eye_on_eye_ellipse'};

sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];
behav_data = import_fun([sub_dir '\\results_video_' num2str(Sub) '_pose_' body_part '_processed.csv']);
%         behav_data = behav_data(:,ismember(behav_data.Properties.VariableNames,varargin_table.cols_to_keep));

behav_data = behav_data(contains(behav_data.Unnamed0,condition),:);
%         behav_data = IBS_clean_video_landmark_data(behav_data);
%         behav_data.timepoints = linspace(0,120,size(behav_data,1))';

if Dyad_no == 10 && strcmp(condition,'NeOcc_3') && (Sub == 1)
    behav_data = array2table(nan(3176,size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
elseif Dyad_no == 19 && strcmp(condition,'NeNoOcc_3')
    behav_data = array2table(nan(3176,size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
    
end

end

function behav_data = IBS_load_manipulated_Video_landmarks(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)


try
    body_part = varargin_table.body_part;
    import_fun = varargin_table.import_fun;
    
catch
    body_part = 'body';
    import_fun = @IBS_import_video_landmarks;
    
end

%         varargin_table.cols_to_keep = {'time_stamps','eye_on_face','eye_on_face_ellipse','eye_on_eye_ellipse'};

sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];

% update 05-10-2022
% old way - this is not nice because of the fact that it adds col names on
% top. instead readtable reads the columns from the file itself
% behav_data = import_fun([sub_dir '\\results_distance_manipulation_video_' num2str(Sub) '_pose_' body_part '.csv']);
%         behav_data = behav_data(:,ismember(behav_data.Properties.VariableNames,varargin_table.cols_to_keep));
%% update 02-12-2022:
% this is rendered important in the case of Dyad 12 when the first frame
% doesn't have data. I tried all possible name pair combinations but the
% only thing that works is to use opts like this. the parameter that seems
% to do the trick is DataLines = [2 Inf]. however u can't give this as name
% pair.




opts = detectImportOptions([sub_dir '\\results_' varargin_table.manipulation ...
    '_manipulation_video_' num2str(Sub) '_pose_' body_part '.csv']);


opts.DataLines = [2 Inf];
% important to have delimiter for cases where there is no data at all!
% this was old way
% behav_data = readtable([sub_dir '\\results_' varargin_table.manipulation  '_manipulation_video_' num2str(Sub) '_pose_' body_part '.csv'],'Delimiter',',');
behav_data = readtable([sub_dir '\\results_' varargin_table.manipulation ...
    '_manipulation_video_' num2str(Sub) '_pose_' body_part '.csv'],opts);

switch(varargin_table.body_part)
    case 'body'
        if size(behav_data,2) ~= 55
            
            behav_data = readtable([sub_dir '\\results_' varargin_table.manipulation  '_manipulation_video_' num2str(Sub) '_pose_' body_part '.csv'],'Delimiter',',');
            
            assert(size(behav_data,2) == 55)
        end
        
end




% so this is problematic because there are cases where the images are not
% loaded correctly and hence there are missing frames. To overcome this,
% this is what I do:
% behav_data = behav_data(contains(behav_data.Condition_block,condition),:);
behav_data = behav_data(contains(behav_data.Var1,condition),:);

behav_data(isnan(behav_data.Dyd_no),{'Dyd_no','Sub_no','Condition','Condition_block'}) = repmat(behav_data(1,{'Dyd_no','Sub_no','Condition','Condition_block'}),sum(isnan(behav_data.Dyd_no)),1);


%         behav_data = IBS_clean_video_landmark_data(behav_data);
%         behav_data.timepoints = linspace(0,120,size(behav_data,1))';

if Dyad_no == 18 && strcmp(condition,'FaOcc_1')
    behav_data = array2table(nan(3101,size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
elseif Dyad_no == 16 && strcmp(condition,'FaNoOcc_3') && Sub == 1
    behav_data = array2table(nan(3112,size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
    
end


%% this is important so that I don't get errors later
sub_dir = IBS_get_sub_behavior_dir('Video_landmarks',Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];
timestamps = IBS_import_video_timestamps(Sub,condition,sub_dir);

% assert(size(behav_data,1)==size(timestamps,2))

if size(behav_data,1)~=size(timestamps,2)
    behav_data = array2table(nan(size(timestamps,2),size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
    
end




end

function behav_data = IBS_load_Video_manual_labelled(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)

import_fun = @IBS_load_video_manual_labelled;
%         varargin_table.label_types = {'Smile','Hands_Feet','Mov_head'};
try
    varargin_table.label_retain_types = varargin_table.label_retain_types;
catch
    varargin_table.label_retain_types = {'Time','Smile','All_mov'};
    varargin_table.label_types = {'Smile','Hands_Feet','Mov_head','Mov_arms','Mov_trunk','Mov_legs'};
    
end
% varargin_table.label_retain_types = {'Time','All_mov'};


mapObj = containers.Map({0,1},{1,2});

behav_data = import_fun(Dyad_no,mapObj(Sub),condition,varargin_table.label_types);

behav_data = behav_data(:,contains(behav_data.Properties.VariableNames,varargin_table.label_retain_types));


end



function behav_data = IBS_load_ECG(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)

import_fun = @IBS_import_ECG_EDA;

varargin_table.cols_to_keep = {'ECG_Rate','ECG_Phase_Completion_Ventricular','ECG_Phase_Completion_Atrial'};

sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
warning('not using input raw dir')
sub_dir = [ sub_dir];
behav_data = import_fun([sub_dir '\\results_ecg_' num2str(Sub) '_' condition '.csv']);
behav_data = behav_data(:,ismember(behav_data.Properties.VariableNames,varargin_table.cols_to_keep));
behav_data = behav_data(2:end,:); % in python there seems to be 1 row more.



behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);

behav_data = rmmissing(behav_data);
load('ecg_time')
behav_data.time_stamps = ecg_time';

%% outlier correction
% behav_data{abs(behav_data.(['ECG_Rate_' num2str(Sub)]))>(mean(behav_data.(['ECG_Rate_' num2str(Sub)])) + 2.5*std(behav_data.(['ECG_Rate_' num2str(Sub)]))),['ECG_Rate_' num2str(Sub)]} = NaN;
%% outlier removal
if ismember('outlier_removal',varargin_table.Properties.VariableNames)    
    dep_measure = [varargin_table.cols_to_keep{1,1} '_'];
    behav_data{abs(behav_data.([dep_measure num2str(Sub)]))>(mean(behav_data.([dep_measure num2str(Sub)])) + ...
        varargin_table.outlier_removal*std(behav_data.([dep_measure num2str(Sub)]))),[dep_measure num2str(Sub)]} = NaN;
end


end




function behav_data = IBS_load_norm_ECG(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)

import_fun = @IBS_import_ECG_EDA;

varargin_table.cols_to_keep = {'ECG_Rate','ECG_Phase_Completion_Ventricular','ECG_Phase_Completion_Atrial'};

sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
warning('not using input raw dir')
sub_dir = [ sub_dir];
behav_data = import_fun([sub_dir '\\results_ecg_' num2str(Sub) '_' condition '.csv']);
behav_data = behav_data(:,ismember(behav_data.Properties.VariableNames,varargin_table.cols_to_keep));
behav_data = behav_data(2:end,:); % in python there seems to be 1 row more.


behav_data_Baseline_Start_1 = import_fun([sub_dir '\\results_ecg_' num2str(Sub) '_Baseline-Start_1.csv']);
behav_data_Baseline_Start_2 = import_fun([sub_dir '\\results_ecg_' num2str(Sub) '_Baseline-Start_2.csv']);
behav_data_Baseline_Start_3 = import_fun([sub_dir '\\results_ecg_' num2str(Sub) '_Baseline-Start_3.csv']);
behav_data_Baseline_End_1 = import_fun([sub_dir '\\results_ecg_' num2str(Sub) '_Baseline-End_1.csv']);
behav_data_Baseline_End_2 = import_fun([sub_dir '\\results_ecg_' num2str(Sub) '_Baseline-End_2.csv']);
behav_data_Baseline_End_3 = import_fun([sub_dir '\\results_ecg_' num2str(Sub) '_Baseline-End_3.csv']);



behav_data_Baseline = mean([mean(behav_data_Baseline_Start_1.ECG_Rate),...
    mean(behav_data_Baseline_Start_2.ECG_Rate),...
    mean(behav_data_Baseline_Start_3.ECG_Rate),...
    mean(behav_data_Baseline_End_1.ECG_Rate),...
    mean(behav_data_Baseline_End_2.ECG_Rate),...
    mean(behav_data_Baseline_End_3.ECG_Rate)]);


behav_data.ECG_Rate = behav_data.ECG_Rate./behav_data_Baseline;

behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);

behav_data = rmmissing(behav_data);
load('ecg_time')
behav_data.time_stamps = ecg_time';

%% outlier correction
% behav_data{abs(behav_data.(['ECG_Rate_' num2str(Sub)]))>(mean(behav_data.(['ECG_Rate_' num2str(Sub)])) + 2.5*std(behav_data.(['ECG_Rate_' num2str(Sub)]))),['ECG_Rate_' num2str(Sub)]} = NaN;



end


function behav_data = IBS_load_resp(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)

import_fun = @IBS_import_resp;

varargin_table.cols_to_keep = {'RSP_Rate','RSP_Amplitude'};

sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
warning('not using input raw dir')
sub_dir = [ sub_dir];
behav_data = import_fun([sub_dir '\\results_resp_' num2str(Sub) '_' condition '.csv']);
behav_data = behav_data(:,ismember(behav_data.Properties.VariableNames,varargin_table.cols_to_keep));
behav_data = behav_data(2:end,:); % in python there seems to be 1 row more.
behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);

behav_data = rmmissing(behav_data);
load('ecg_time')
behav_data.time_stamps = ecg_time';

%% outlier removal
if ismember('outlier_removal',varargin_table.Properties.VariableNames)    
    dep_measure = [varargin_table.cols_to_keep{1,1} '_'];
    behav_data{abs(behav_data.([dep_measure num2str(Sub)]))>(mean(behav_data.([dep_measure num2str(Sub)])) + ...
        varargin_table.outlier_removal*std(behav_data.([dep_measure num2str(Sub)]))),[dep_measure num2str(Sub)]} = NaN;
end
end


function behav_data =IBS_load_EDA(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)

import_fun = @IBS_import_ECG_EDA;

varargin_table.cols_to_keep = {'EDA_Tonic','EDA_Phasic'};
% varargin_table.cols_to_keep = {'EDA_Tonic','EDA_Phasic'};

sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
warning('not using input raw dir')
sub_dir = [ sub_dir];
behav_data = import_fun([sub_dir '\\results_ecg_' num2str(Sub) '_' condition '.csv']);
behav_data = behav_data(:,ismember(behav_data.Properties.VariableNames,varargin_table.cols_to_keep));
behav_data = behav_data(2:end,:); % in python there seems to be 1 row more.
behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);

behav_data = rmmissing(behav_data);
load('ecg_time')
behav_data.time_stamps = ecg_time';

%% outlier removal
if ismember('outlier_removal',varargin_table.Properties.VariableNames)    
    dep_measure = [varargin_table.cols_to_keep{1,1} '_'];
    behav_data{abs(behav_data.([dep_measure num2str(Sub)]))>(mean(behav_data.([dep_measure num2str(Sub)])) + ...
        varargin_table.outlier_removal*std(behav_data.([dep_measure num2str(Sub)]))),[dep_measure num2str(Sub)]} = NaN;
end
end



function behav_data = IBS_load_smile_auto(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)

Subs = [0,1];
video_sub = setdiff(Subs,Sub);
% video_sub = Sub;
import_fun = @IBS_import_smile_detection;
% varargin_table.cols_to_keep = {'time_stamps_rel','smile','frame_no','smile_pca','smile_sd','mouth_size','mouth_size_inner','mouth_size_x','mouth_size_y'};
% varargin_table.cols_to_keep = {'time_stamps_rel','smile_sd'};
varargin_table.cols_to_keep = {'time_stamps_rel','mouth_size'};

sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,video_sub,condition);

% sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];

behav_data = import_fun([sub_dir '\\results_smile_detection_from_eye_tracker_pca.csv']);
%                 behav_data = import_fun([sub_dir '\\results_smile_detection_from_eye_tracker_pca_Near.csv']);
%         behav_data = IBS_smile_auto_mouth_size(behav_data);
% behav_data = IBS_smile_auto_mouth_size(behav_data,sub_dir);




% behav_data.mouth_size(behav_data.mouth_size <=0) = nan;
%
% behav_data.mouth_size(behav_data.mouth_size > nanmean(behav_data.mouth_size)+3*nanstd(behav_data.mouth_size)) = nan;
% behav_data.mouth_size(behav_data.mouth_size < nanmean(behav_data.mouth_size)-3*nanstd(behav_data.mouth_size)) = nan;

% this is necessary for find sudden change function - it doesn't handle
% missing data
% behav_data.mouth_size(isnan(behav_data.mouth_size)) = nanmean(behav_data.mouth_size);
% this seemed to preserve a bit the latency



% %%
%% Filt_2 : Remove outliers (filter too high or low data points)
[ indices ]             = Giac_StdOverMeanVectorFilter( behav_data.mouth_size', 3 );   %3 std away from m
behav_data.mouth_size(indices)       = nan;
%
%
% %% Filt_3: Interpolate nans (using Dav_interp_nans function)
%
behav_data.mouth_size                     = Dav_interp_nans(behav_data.mouth_size);

%%
% behav_data.mouth_size = smoothdata(behav_data.mouth_size,'gaussian',1,'SamplePoints',behav_data.time_stamps_rel);
standard_time_stamps =linspace(behav_data.time_stamps_rel(1),behav_data.time_stamps_rel(end),length(behav_data.time_stamps_rel));
behav_data.mouth_size = interp1(behav_data.time_stamps_rel,behav_data.mouth_size,standard_time_stamps)';

frame_rate = length(behav_data.time_stamps_rel)/(behav_data.time_stamps_rel(end)-behav_data.time_stamps_rel(1));

behav_data.mouth_size = movingmean(behav_data.mouth_size,round(frame_rate)*1,1,[]);
% behav_data.mouth_size = smoothdata(behav_data.mouth_size,'sgolay','SamplePoints',standard_time_stamps);
% behav_data.mouth_size = smoothdata(behav_data.mouth_size,'lowess',1,'SamplePoints',standard_time_stamps);
% behav_data.mouth_size = smoothdata(behav_data.mouth_size,'loess',1,'SamplePoints',standard_time_stamps);
% behav_data.mouth_size = smoothdata(behav_data.mouth_size,'sgolay',1,'SamplePoints',standard_time_stamps);
% behav_data.mouth_size = smoothdata(behav_data.mouth_size,'movmean',1,'SamplePoints',standard_time_stamps);
% behav_data.mouth_size = movmean(behav_data.mouth_size,[0.15 0],'SamplePoints',standard_time_stamps);
% behav_data.mouth_size = movmean(behav_data.mouth_size,0.1,'SamplePoints',standard_time_stamps);
% behav_data.mouth_size = movmean(behav_data.mouth_size,[0 0.2],'SamplePoints',standard_time_stamps);
% behav_data.mouth_size = movmean(behav_data.mouth_size,0.1,'SamplePoints',standard_time_stamps);
% behav_data.mouth_size = movmean(behav_data.mouth_size,[0.3 0],'SamplePoints',standard_time_stamps);
% behav_data.mouth_size = movmean(behav_data.mouth_size,[0 0.1],'SamplePoints',standard_time_stamps);


%% alternate way

% this procedure is very specifically used for being similar across all the
% movement analysis obtained from openpose. Davide and Giacomo do the same
% stuff for the x and y in movement data. even the windowsize is the same.
% Although I have reservations about using moving window, as it moves the
% time, it doesn't seem to shift the curve too much.
% window_time = 1; %s
% behav_data.mouth_size_x(behav_data.mouth_size_x <=0) = nan;
% behav_data.mouth_size_y(behav_data.mouth_size_y <=0) = nan;
% % maybe doesn't make sense to remove outliers as smile is rare event
% % behav_data.mouth_size_x(behav_data.mouth_size_x >nanmean(behav_data.mouth_size_x)+ 3*nanstd(behav_data.mouth_size_x)) = nan;
% % behav_data.mouth_size_y(behav_data.mouth_size_y >nanmean(behav_data.mouth_size_y)+ 3*nanstd(behav_data.mouth_size_y)) = nan;
%
% % this procedure ensures that the length of the data remains the same as
% % the input as well as I can do the moving window over 1 sec based on
% % computed frame rate. Here I don't assume that the frames are captured
% % equally but I input the data from the timestamps.
% standard_time_stamps =linspace(behav_data.time_stamps_rel(1),behav_data.time_stamps_rel(end),length(behav_data.time_stamps_rel));
% frame_rate = length(behav_data.time_stamps_rel)/(behav_data.time_stamps_rel(end)-behav_data.time_stamps_rel(1));
%
% behav_data.mouth_size_x = interp1(behav_data.time_stamps_rel,behav_data.mouth_size_x,standard_time_stamps)';
% behav_data.mouth_size_y = interp1(behav_data.time_stamps_rel,behav_data.mouth_size_y,standard_time_stamps)';
% % this ensures that the mouth size is the same in joint and joint XOR and
% % joint AND conditions.
% behav_data.mouth_size = sqrt(smoothdata(behav_data.mouth_size_x,'movmean',frame_rate*window_time).^2 + smoothdata(smoothdata(behav_data.mouth_size_y,'movmean',frame_rate*window_time)).^2);
%%
switch(varargin_table.behav_analysis)
    
    case 'joint_XOR'
        
        behav_data = IBS_smile_auto_mouth_size(behav_data,sub_dir);
    case 'joint_AND'
        behav_data = IBS_smile_auto_mouth_size(behav_data,sub_dir);
        
        
end

%%
% don't correct this for checking smile correspondences

behav_data_starttime = IBS_import_eye_tracker_blink_info([sub_dir '\\export_info.csv']);
%         behav_data.time_stamps = behav_data.time_stamps - behav_data.time_stamps(1,1);
if Dyad_no == 7 && video_sub == 1 && strcmp(condition,'FaNoOcc_3')
    % there is a delay of 28.4s as the acquisition started at that point
    %     behav_data.time_stamps_rel = ((behav_data.time_stamps_raw - behav_data_starttime)+ 28.4) + IBS_get_behav_delay(Dyad_no,Sub,condition);
    behav_data.time_stamps_rel = ((behav_data.time_stamps_raw - behav_data_starttime)+ 28.4) + IBS_get_behav_delay(Dyad_no,video_sub,condition);
    
elseif contains(condition,{'Task_1','Task_2','Task_3'})
    
    %     behav_data.time_stamps_rel = (behav_data.time_stamps_raw - behav_data_starttime) + IBS_get_behav_delay_task(Dyad_no,Sub,condition);
    behav_data.time_stamps_rel = (behav_data.time_stamps_raw - behav_data_starttime) + IBS_get_behav_delay_task(Dyad_no,video_sub,condition);
    
else
    %     behav_data.time_stamps_rel = (behav_data.time_stamps_raw - behav_data_starttime)+ IBS_get_behav_delay(Dyad_no,Sub,condition);
    behav_data.time_stamps_rel = (behav_data.time_stamps_raw - behav_data_starttime)+ IBS_get_behav_delay(Dyad_no,video_sub,condition);
    
    
    
end
%%

% %         behav_data.time_stamps = behav_data.time_stamps - behav_data.time_stamps(1,1);
% if Dyad_no == 7 && Sub == 1 && strcmp(condition,'FaNoOcc_3')
%     % there is a delay of 28.4s as the acquisition started at that point
%     behav_data.time_stamps_rel = (behav_data.time_stamps_rel + 28.4) + IBS_get_behav_delay(Dyad_no,Sub,condition);
%
% elseif Dyad_no == 11 && Sub == 1 && strcmp(condition,'Task_1')
%         behav_data.time_stamps_rel = (behav_data.time_stamps_rel+ 39.668) + IBS_get_behav_delay(Dyad_no,Sub,condition);
%
% else
%     behav_data.time_stamps_rel = behav_data.time_stamps_rel+ IBS_get_behav_delay(Dyad_no,Sub,condition);
% end
%
%%

behav_data = behav_data(:,ismember(behav_data.Properties.VariableNames,varargin_table.cols_to_keep));

behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);


end


function behav_data = IBS_load_video_openpose_landmarks_manual_cleaned(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)


sub_dir = IBS_get_sub_behavior_dir('Video_landmarks',Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];
timestamps = IBS_import_video_timestamps(Sub,condition,sub_dir);

% import_fun = @IBS_load_video_landmarks_manual_cleaned;
import_fun = @IBS_load_video_landmarks_manual_cleaned_abs;
try
    labels = varargin_table.labels;
catch
    labels = 'ALL';
    %     labels = {'LAnkle','LBigToe','LSmallToe','LHeel','RAnkle','RBigToe','RSmallToe', 'RHeel'};
    
    %     labels = {'Nose','REye','LEye','REar','LEar'};
    %     labels = {'RKnee','LKnee'};
    %     labels = {'LWrist','RWrist'};
    %     labels = {'RElbow','LElbow'};
    %     labels = {'RHip','LHip'};
    %     labels = {'LShoulder','RShoulder'};
end
varargin_table.cols_to_keep = {'timestamps','ALL'};
try
    raw_data_dir = varargin_table.raw_data_dir;
catch
    raw_data_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
    
end
behav_data = import_fun(Dyad_no,Sub,condition,raw_data_dir,labels,varargin_table);
% behav_data = import_fun(Dyad_no,Sub,condition,raw_data_dir);

% this is because for these conditions, there wasn't any movement data (all
% nans), so davide created fake data of all zeros with length = 25 fps *
% 120 sec = 3000 fake frames.

% Davide:
% I substituted these rejected files with nans for both S0 and S1.
% The problem of missing frames is probably due to the fact that I set the
% length of these files equal to 120 (time in sec)* 25 (frames per secs) = 3000,
% instead of their original length
if Dyad_no == 7 && strcmp(condition,'FaNoOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
elseif Dyad_no == 15 && strcmp(condition,'FaNoOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 19 && strcmp(condition,'NeNoOcc_3')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 9 && strcmp(condition,'FaOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
elseif Dyad_no == 7 && strcmp(condition,'FaOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
    
%% this was uncommented pre 15-02-2023 for manipulated_Video_landmarks analysis
% for this subject, there are frames that didn't load and so I had to
% equate it by putting nans equal to the subs. otherwise this doesn't work.
% uncomment for manipulated analysis

% elseif Dyad_no == 18 && strcmp(condition,'FaOcc_1')
%     behav_data = array2table(nan(3100,size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
    
else
    
    behav_data.timestamps = [timestamps(2:end)]';
end



switch(varargin_table.behav_analysis)
    case 'joint_XOR'
        
        behav_data.ALL = logical(1*(abs(behav_data.ALL) > nanmean(behav_data.ALL) + nanstd(behav_data.ALL)));
        
    case 'joint_AND'
        behav_data.ALL = logical(1*(abs(behav_data.ALL) > nanmean(behav_data.ALL) + nanstd(behav_data.ALL)));
        
        
        
end


behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);



end


function behav_data = IBS_load_video_openpose_landmarks_manual_cleaned_in_vicon_units(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)


sub_dir = IBS_get_sub_behavior_dir('Video_landmarks',Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];
timestamps = IBS_import_video_timestamps(Sub,condition,sub_dir);

% import_fun = @IBS_load_video_landmarks_manual_cleaned;
import_fun = @IBS_load_video_landmarks_manual_cleaned_abs;
try
    labels = varargin_table.labels;
catch
    labels = 'ALL';
    %     labels = {'LAnkle','LBigToe','LSmallToe','LHeel','RAnkle','RBigToe','RSmallToe', 'RHeel'};
    
    %     labels = {'Nose','REye','LEye','REar','LEar'};
    %     labels = {'RKnee','LKnee'};
    %     labels = {'LWrist','RWrist'};
    %     labels = {'RElbow','LElbow'};
    %     labels = {'RHip','LHip'};
    %     labels = {'LShoulder','RShoulder'};
end
varargin_table.cols_to_keep = {'timestamps','ALL'};
try
    raw_data_dir = varargin_table.raw_data_dir;
catch
    raw_data_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
    
end
varargin_table.output_type = 'absolute';
behav_data = import_fun(Dyad_no,Sub,condition,raw_data_dir,labels,varargin_table);

%% 22-09 update for ICC
%% from eucl
% vicon_conversion = 2.3109; % from IBS_compute_vicon_video_conversion.m 78 cm table old value 0.2311 (assuming vicon in cm)
% behav_data.ALL = behav_data.ALL * vicon_conversion;
%% from internet as to how much a 96 dpi image should have
%% confirm that even with diff u need to multiply the factor
% a = 25*rand(1,100);
% b = a;
% b = a/2.5+rand(1,100)/10;
% icc21([diff(a)' diff(b)']) % 0.6919
% icc21([diff(a)' (diff(b).*2.5)']) % 1


vicon_conversion = 0.264583333; % 
behav_data.ALL = behav_data.ALL * vicon_conversion;


% behav_data = import_fun(Dyad_no,Sub,condition,raw_data_dir);

% this is because for these conditions, there wasn't any movement data (all
% nans), so davide created fake data of all zeros with length = 25 fps *
% 120 sec = 3000 fake frames.

% Davide:
% I substituted these rejected files with nans for both S0 and S1.
% The problem of missing frames is probably due to the fact that I set the
% length of these files equal to 120 (time in sec)* 25 (frames per secs) = 3000,
% instead of their original length
if Dyad_no == 7 && strcmp(condition,'FaNoOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
elseif Dyad_no == 15 && strcmp(condition,'FaNoOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 19 && strcmp(condition,'NeNoOcc_3')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 9 && strcmp(condition,'FaOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
elseif Dyad_no == 7 && strcmp(condition,'FaOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
    
%% this was uncommented pre 15-02-2023 for manipulated_Video_landmarks analysis
% for this subject, there are frames that didn't load and so I had to
% equate it by putting nans equal to the subs. otherwise this doesn't work.
% uncomment for manipulated analysis

% elseif Dyad_no == 18 && strcmp(condition,'FaOcc_1')
%     behav_data = array2table(nan(3100,size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
    
else
    
    behav_data.timestamps = [timestamps(2:end)]';
end



switch(varargin_table.behav_analysis)
    case 'joint_XOR'
        
        behav_data.ALL = logical(1*(abs(behav_data.ALL) > nanmean(behav_data.ALL) + nanstd(behav_data.ALL)));
        
    case 'joint_AND'
        behav_data.ALL = logical(1*(abs(behav_data.ALL) > nanmean(behav_data.ALL) + nanstd(behav_data.ALL)));
        
        
        
end


behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);



end

function behav_data = IBS_load_manipulated_video_openpose_landmarks_manual_cleaned(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)


sub_dir = IBS_get_sub_behavior_dir('manipulated_Video_landmarks',Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];
timestamps = IBS_import_video_timestamps(Sub,condition,sub_dir);

% import_fun = @IBS_load_video_landmarks_manual_cleaned;
import_fun = @IBS_load_manipulated_video_landmarks_manual_cleaned_abs;
try
    labels = varargin_table.labels;
catch
    %     labels = 'ALL';
    %     labels = {'LAnkle','LBigToe','LSmallToe','LHeel','RAnkle','RBigToe','RSmallToe', 'RHeel'};
    
    %     labels = {'Nose','REye','LEye','REar','LEar'};
    %     labels = {'RKnee','LKnee'};
    %     labels = {'LWrist','RWrist'};
    %     labels = {'RElbow','LElbow'};
    %     labels = {'RHip','LHip'};
    labels = {'LShoulder','RShoulder'};
end
varargin_table.cols_to_keep = {'timestamps','ALL'};

% try
%     raw_data_dir = varargin_table.raw_data_dir;
% catch
%     raw_data_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
% end


behav_data = import_fun(Dyad_no,Sub,condition,raw_data_dir,labels,varargin_table);
% behav_data = import_fun(Dyad_no,Sub,condition,raw_data_dir);

% this is because for these conditions, there wasn't any movement data (all
% nans), so davide created fake data of all zeros with length = 25 fps *
% 120 sec = 3000 fake frames.

% Davide:
% I substituted these rejected files with nans for both S0 and S1.
% The problem of missing frames is probably due to the fact that I set the
% length of these files equal to 120 (time in sec)* 25 (frames per secs) = 3000,
% instead of their original length
if Dyad_no == 7 && strcmp(condition,'FaNoOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
elseif Dyad_no == 15 && strcmp(condition,'FaNoOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 19 && strcmp(condition,'NeNoOcc_3')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 9 && strcmp(condition,'FaOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
elseif Dyad_no == 7 && strcmp(condition,'FaOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
else
    
    behav_data.timestamps = [timestamps(2:end)]';
end



switch(varargin_table.behav_analysis)
    case 'joint_XOR'
        
        behav_data.ALL = logical(1*(abs(behav_data.ALL) > nanmean(behav_data.ALL) + nanstd(behav_data.ALL)));
        
    case 'joint_AND'
        behav_data.ALL = logical(1*(abs(behav_data.ALL) > nanmean(behav_data.ALL) + nanstd(behav_data.ALL)));
        
        
        
end


behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);



end
function behav_data = IBS_load_video_openpose_landmarks_manual_cleaned_add_body(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)


sub_dir = IBS_get_sub_behavior_dir('Video_landmarks',Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];
timestamps = IBS_import_video_timestamps(Sub,condition,sub_dir);

% import_fun = @IBS_load_video_landmarks_manual_cleaned;
import_fun = @IBS_load_video_landmarks_manual_cleaned_abs;
try
    labels = varargin_table.labels;
catch
    %     labels = 'ALL';
    labels = {'LAnkle','LBigToe','LSmallToe','LHeel','RAnkle','RBigToe','RSmallToe', 'RHeel'};
    
    %     labels = {'Nose','REye','LEye','REar','LEar'};
    %     labels = {'RKnee','LKnee'};
    %     labels = {'LWrist','RWrist'};
    %     labels = {'RElbow','LElbow'};
    %     labels = {'RHip','LHip'};
    %     labels = {'Nose','REye','LEye','REar','LEar'};
end
varargin_table.cols_to_keep = {'timestamps','ALL'};
try
    raw_data_dir = varargin_table.raw_data_dir;
catch
    raw_data_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
    
end
behav_data = import_fun(Dyad_no,Sub,condition,raw_data_dir,labels);
% behav_data = import_fun(Dyad_no,Sub,condition,raw_data_dir);

% this is because for these conditions, there wasn't any movement data (all
% nans), so davide created fake data of all zeros with length = 25 fps *
% 120 sec = 3000 fake frames.

% Davide:
% I substituted these rejected files with nans for both S0 and S1.
% The problem of missing frames is probably due to the fact that I set the
% length of these files equal to 120 (time in sec)* 25 (frames per secs) = 3000,
% instead of their original length
if Dyad_no == 7 && strcmp(condition,'FaNoOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
elseif Dyad_no == 15 && strcmp(condition,'FaNoOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 19 && strcmp(condition,'NeNoOcc_3')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 9 && strcmp(condition,'FaOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
elseif Dyad_no == 7 && strcmp(condition,'FaOcc_1')
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
else
    
    behav_data.timestamps = [timestamps(2:end)]';
end



switch(varargin_table.behav_analysis)
    case 'joint_XOR'
        
        behav_data.ALL = logical(1*(abs(behav_data.ALL) > nanmean(behav_data.ALL) + nanstd(behav_data.ALL)));
        
    case 'joint_AND'
        behav_data.ALL = logical(1*(abs(behav_data.ALL) > nanmean(behav_data.ALL) + nanstd(behav_data.ALL)));
        
        
        
end


behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);


behav_data.Properties.VariableNames = strrep(behav_data.Properties.VariableNames,'ALL','Other');

end


function behav_data = IBS_load_abs_video_openpose_landmarks(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)


sub_dir = IBS_get_sub_behavior_dir('Video_landmarks',Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];

timestamps = IBS_import_video_timestamps(Sub,condition,sub_dir);

import_fun = @IBS_load_abs_video_landmarks;


try
    labels = varargin_table.labels;
catch
    labels = 'ALL';
    %     labels = {'LAnkle','LBigToe','LSmallToe','LHeel','RAnkle','RBigToe','RSmallToe', 'RHeel'};
    
    %     labels = {'Nose','REye','LEye','REar','LEar'};
    %     labels = {'RKnee','LKnee'};
    %     labels = {'LWrist','RWrist'};
    %     labels = {'RElbow','LElbow'};
    %     labels = {'RHip','LHip'};
    %     labels = {'LShoulder','RShoulder'};
end
varargin_table.cols_to_keep = {'timestamps','ALL'};


try
    raw_data_dir = varargin_table.raw_data_dir;
catch
    raw_data_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
    
end
behav_data = import_fun(Dyad_no,Sub,condition,raw_data_dir,labels,varargin_table);
% behav_data = import_fun(Dyad_no,Sub,condition,raw_data_dir);

% this is because for these conditions, there wasn't any movement data (all
% nans), so davide created fake data of all zeros with length = 25 fps *
% 120 sec = 3000 fake frames.

% Davide:
% I substituted these rejected files with nans for both S0 and S1.
% The problem of missing frames is probably due to the fact that I set the
% length of these files equal to 120 (time in sec)* 25 (frames per secs) = 3000,
% instead of their original length
if Dyad_no == 7 && strcmp(condition,'FaNoOcc_1')
    % overwrite the loaded data because this is skipped trial based on
    % Davide cleaning
    behav_data = array2table(nan(3000,size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
    
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
elseif Dyad_no == 15 && strcmp(condition,'FaNoOcc_1')
    behav_data = array2table(nan(3000,size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
    
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 19 && strcmp(condition,'NeNoOcc_3')
    behav_data = array2table(nan(3000,size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
    
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 9 && strcmp(condition,'FaOcc_1')
    behav_data = array2table(nan(3000,size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
    
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 7 && strcmp(condition,'FaOcc_1')
    behav_data = array2table(nan(3000,size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
    
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 10 && strcmp(condition,'NeOcc_2')
    behav_data = array2table(nan(3000,size(behav_data,2)),'VariableNames',behav_data.Properties.VariableNames);
    
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
elseif Dyad_no == 10 && strcmp(condition,'NeOcc_3') && (Sub == 1)
    
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
elseif Dyad_no == 19 && strcmp(condition,'NeNoOcc_3') && (Sub == 1)
    behav_data.timestamps = linspace(0,120,size(behav_data,1))';
    
else
    
    %     behav_data.timestamps = [timestamps(2:end)]';
    % use all timestamps as it's not diff
    behav_data.timestamps = timestamps';
    
end



switch(varargin_table.behav_analysis)
    case 'joint_XOR'
        
        behav_data.ALL = logical(1*(abs(behav_data.ALL) > nanmean(behav_data.ALL) + nanstd(behav_data.ALL)));
        
    case 'joint_AND'
        behav_data.ALL = logical(1*(abs(behav_data.ALL) > nanmean(behav_data.ALL) + nanstd(behav_data.ALL)));
        
        
        
end


behav_data.Properties.VariableNames = cellfun(@(x) [x '_' num2str(Sub)],behav_data.Properties.VariableNames,'UniformOutput',false);



end


function behav_data = IBS_load_body_landmarks_vicon(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)
if nargin<5
    raw_data_dir = 'D:\\iannettilab_dropbox\\Dropbox\\Koul_Atesh\\IBS\\';
end
import_fun = @IBS_import_vicon_data;
% varargin_table.cols_to_keep = {'time_stamps_rel','smile','frame_no','smile_pca','smile_sd','mouth_size','mouth_size_inner','mouth_size_x','mouth_size_y'};
% varargin_table.cols_to_keep = {'time_stamps_rel','smile_sd'};
switch(Sub)
    case 0
        %         % not including S0_Head4 as there is no equivalent in S1 - was used
        %         % to distinguish the two models.
        %         varargin_table.cols_to_keep = {'timepoints','S0_Head1', 'S0_Head2', 'S0_Head3', 'S0_Torso1', 'S0_Torso2', 'S0_Torso3',...
        %             'S0_Left_Arm1', 'S0_Left_Arm2', 'S0_Left_Arm3', 'S0_Left_Leg1', 'S0_Left_Leg2', 'S0_Left_Leg3',...
        %             'S0_Right_Arm1', 'S0_Right_Arm2', 'S0_Right_Arm3','S0_Right_Leg1', 'S0_Right_Leg2', 'S0_Right_Leg3'};
        %
        varargin_table.cols_to_keep = {'timepoints','S0_Head1', 'S0_Head2', 'S0_Head3', 'S0_Head4','S0_Torso1', 'S0_Torso2', 'S0_Torso3',...
            'S0_Left_Arm1', 'S0_Left_Arm2', 'S0_Left_Arm3', 'S0_Left_Leg1', 'S0_Left_Leg2', 'S0_Left_Leg3',...
            'S0_Right_Arm1', 'S0_Right_Arm2', 'S0_Right_Arm3','S0_Right_Leg1', 'S0_Right_Leg2', 'S0_Right_Leg3'};
        %                 varargin_table.cols_to_keep = {'timepoints','S0_Head1', 'S0_Head2', 'S0_Head3', 'S0_Head4', 'S0_Torso1', 'S0_Torso2', 'S0_Torso3',...
        %             'S0_Left_Arm1', 'S0_Left_Arm2', 'S0_Left_Arm3', 'S0_Left_Leg1', 'S0_Left_Leg2', 'S0_Left_Leg3',...
        %             'S0_Right_Arm1', 'S0_Right_Arm2', 'S0_Right_Arm3','S0_Right_Leg1', 'S0_Right_Leg2', 'S0_Right_Leg3'};
    case 1
        varargin_table.cols_to_keep = {'timepoints','S1_Head1', 'S1_Head2', 'S1_Head3', 'S1_Torso1', 'S1_Torso2', 'S1_Torso3',...
            'S1_Left_Arm1', 'S1_Left_Arm2', 'S1_Left_Arm3',...
            'S1_Right_Arm1', 'S1_Right_Arm2', 'S1_Right_Arm3',...
            'S1_Right_Leg1', 'S1_Right_Leg2', 'S1_Right_Leg3',...
            'S1_Left_Leg1', 'S1_Left_Leg2', 'S1_Left_Leg3'};
end
sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);

% sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];
if (Dyad_no == 17 && strcmp(condition,'FaOcc_2')) || (Dyad_no == 18 && strcmp(condition,'FaNoOcc_1')) || (Dyad_no == 7 && strcmp(condition,'FaOcc_2'))
    behav_data = array2table(nan(100,length(varargin_table.cols_to_keep)),'VariableNames',varargin_table.cols_to_keep);
    behav_data.timepoints = linspace(0,120,100)';
    % elseif contains(condition,{'Task_1','Task_2','Task_3'})
    %
    %     %     behav_data.time_stamps_rel = (behav_data.time_stamps_raw - behav_data_starttime) + IBS_get_behav_delay_task(Dyad_no,Sub,condition);
    %     behav_data.time_stamps_rel = (behav_data.time_stamps_raw - behav_data_starttime) + IBS_get_behav_delay_task(Dyad_no,video_sub,condition);
    %
else
    %     if (Dyad_no ==14) || (Dyad_no ==16)  || (Dyad_no ==19)|| (Dyad_no ==21) || (Dyad_no ==23)|| (Dyad_no ==8)
    % the way some vicon
    if ismember(Dyad_no,[3,4,5,6,8,14,16,17,19,21])
        %         behav_data = import_fun([sub_dir sprintf('Dyd_%0.3d',Dyad_no) '_' condition '.csv'],@IBS_import_vicon_data_S1_S0,varargin_table);
        varargin_table.import_fun = @IBS_import_vicon_data_S1_S0;
        behav_data = import_fun([sub_dir sprintf('Dyd_%0.3d',Dyad_no) '_' condition '.csv'],varargin_table);
        
    else
        %         behav_data = import_fun([sub_dir sprintf('Dyd_%0.3d',Dyad_no) '_' condition '.csv'],@IBS_import_vicon_data_S0_S1);
        varargin_table.import_fun = @IBS_import_vicon_data_S0_S1;
        behav_data = import_fun([sub_dir sprintf('Dyd_%0.3d',Dyad_no) '_' condition '.csv'],varargin_table);
    end
end
behav_data = behav_data(:,ismember(behav_data.Properties.VariableNames,varargin_table.cols_to_keep));
% behav_data.Dyad_no = repmat(Dyad_no,[1 size(behav_data,1)])
end

function vicon_video_body_pose = IBS_load_vicon_video_landmark(behavior,Dyad_no,Sub,condition,raw_data_dir,varargin_table)
if nargin<5
    raw_data_dir = 'D:\\iannettilab_dropbox\\Dropbox\\Koul_Atesh\\IBS\\';
end
import_fun = @IBS_import_vicon_video_landmark;

sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);

% sub_dir = IBS_get_sub_behavior_dir(behavior,Dyad_no,Sub,condition);
sub_dir = [raw_data_dir sub_dir];


vicon_video_body_pose = import_fun([sub_dir 'results_' condition '_vicon_video_pose_body.csv']);

%% need to remove the approx 1 sec before recording
trigger_fname = [sub_dir sprintf('Dyd_%0.3d',Dyad_no) '_' condition '.txt'];
[triggers,last_frame] = import_vicon_trigger(trigger_fname);

vicon_video_fs = 25;

frames_to_skip = round(triggers(1)/vicon_video_fs);
frames_to_clip = round((last_frame-triggers(2))/vicon_video_fs);
% remove the initial
vicon_video_body_pose = vicon_video_body_pose(frames_to_skip:(end-frames_to_clip),:);


vicon_video_body_pose.timepoints = linspace(0,size(vicon_video_body_pose,1)/vicon_video_fs,size(vicon_video_body_pose,1))';
vicon_video_body_pose = movevars(vicon_video_body_pose,'timepoints','Before','S1_Nose_x');

%%
% varargin_table.cols_to_keep = {'time_stamps_rel','smile','frame_no','smile_pca','smile_sd','mouth_size','mouth_size_inner','mouth_size_x','mouth_size_y'};
% varargin_table.cols_to_keep = {'time_stamps_rel','smile_sd'};

markers =["Nose","Neck","RShoulder","RElbow","RWrist","LShoulder","LElbow",...
    "LWrist","MidHip","RHip","RKnee","RAnkle","LHip","LKnee","LAnkle",...
    "REye","LEye","REar","LEar","LBigToe","LSmallToe","LHeel","RBigToe",...
    "RSmallToe", "RHeel"];
switch(Sub)
    case 0
        varargin_table.cols_to_keep = cat(2,{'timepoints'},...
            strcat(repmat("S0_",1,length(markers)),markers),...
            {'Dyd_no' ,'Condition','Condition_block'});
    case 1
        varargin_table.cols_to_keep = cat(2,{'timepoints'},...
            strcat(repmat("S1_",1,length(markers)),markers),...
            {'Dyd_no' ,'Condition','Condition_block'});
end

% if (Dyad_no == 17 && strcmp(condition,'FaOcc_2')) || (Dyad_no == 18 && strcmp(condition,'FaNoOcc_1'))
%     behav_data = array2table(nan(100,length(varargin_table.cols_to_keep)),'VariableNames',varargin_table.cols_to_keep);
%     behav_data.timepoints = linspace(0,120,100)';
%     % elseif contains(condition,{'Task_1','Task_2','Task_3'})
%     %
%     %     %     behav_data.time_stamps_rel = (behav_data.time_stamps_raw - behav_data_starttime) + IBS_get_behav_delay_task(Dyad_no,Sub,condition);
%     %     behav_data.time_stamps_rel = (behav_data.time_stamps_raw - behav_data_starttime) + IBS_get_behav_delay_task(Dyad_no,video_sub,condition);
%     %
% else
%     behav_data = import_fun([sub_dir sprintf('Dyd_%0.3d',Dyad_no) '_' condition '.csv']);
% end
vicon_video_body_pose = vicon_video_body_pose(:,ismember(vicon_video_body_pose.Properties.VariableNames,varargin_table.cols_to_keep));



end