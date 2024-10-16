% Based on Bom MS, Brak AMA, Raemaekers M, Ramsey NF, Vansteensel MJ,
% Branco MP (2024): Large-scale fMRI dataset for the design of motor-based
% Brain-Computer Interfaces https://doi.org/10.1101/2024.07.29.24311044

% doc MrImage

%% Add Paths 
addpath(genpath("tapas-master"));
addpath(genpath("spm12"));
addpath(genpath("Patient 23 - 3T"));


%% Load subject data for 3T Data
func_3T = MrImage('Patient 23 - 3T/sub-023_ses-3T_task-Motor2Class_run-1_bold.nii');
anat_3T = MrImage('Patient 23 - 3T/sub-023_ses-3T_task-mansfield_run-1_T1w.nii');
task_3T = tdfread('Patient 23 - 3T/sub-023_ses-3T_task-Motor2Class_run-1_events.tsv');

%% Big Brain GLM Parameter Pre-processing
% Segment out data of duration and onset that only corresponds with 
% "move"condition

% Data variables 
trial_type = task_3T.trial_type;  % 9x4 cell array
duration = task_3T.duration;      % 9x1 numeric array
onset = task_3T.onset;            % 9x1 numeric array

% New arrays for used onset and duration variables
glm_trial_type = {};  
glm_onset = [];       
glm_duration = [];    

% Loop through each row in trial_type
for i = 1:size(trial_type, 1)
    condition = strtrim(trial_type(i, :));  % Remove any extra spaces

    if strcmp(condition, 'move')
        % Store corresponding trial type, duration, and onset for "move"
        glm_trial_type = [glm_trial_type; {condition}]; 
        glm_duration = [glm_duration; duration(i)];      
        glm_onset = [glm_onset; onset(i)];               
    end
end

% convert glm_trial_type to character array
glm_trial_type = char(glm_trial_type);

%% Load subject data for 7T Data
%func_7T = MrImage('sub-023_ses-7T_task-Mapping5Fingers_run-1_bold.nii');
%anat_7T = MrImage('sub-023_ses-7T_task-t1_run-1_T1w.nii');
%task_7T = tdfread('sub-023_ses-7T_task-Motor2Class_run-1_events.tsv');


%% Visualisation for 3T
% plot of the first five volumes of the functional image time series
func_3T.plot('t', 1:5, 'rotate90', 1);

%% plot of the anatomical image
anat_3T.plot('rotate90', 1);
title('Anatomical Image of Patient');

% interactive plot of anatomical and mean of functional image time series
anat_3T.plot('plotType', 'spmi', 'overlayImages', func_3T.mean);
title('Interactive plot of functional time series patient image');


%% Visualisation for 7T
% plot of the first five volumes of the functional image time series
%func_7T.plot('t', 1:5, 'rotate90', 1);

% plot of the anatomical image
%anat_7T.plot('rotate90', 1);

% interactive plot of anatomical and mean of functional image time series
%anat_7T.plot('plotType', 'spmi', 'overlayImages', func_7T.mean);


%% Quality check
% plot the mean of the functional image time series (1 point)
func_3T.mean(4).plot('rotate90', 1);    %Is this right?
%func_3T.mean(3).plot('rotate90', 1);

% plot the snr of the functional image time series (1 point)
func_3T.snr(4).plot('rotate90', 1);     %Is this right?
%func_3T.snr(3).plot('rotate90', 1);


%% Realignment
% perform relignment (1 point)
% NOTE: For faster relignment, increase separation and reduce quality
% Optimising values:
% 'Separation', 8
% 'quality', 0.6
[r_func, realignment_parameters] = func_3T.realign('quality', 0.6, 'separation', 6);
r_func.plot('rotate90', 1); 

%% Quality check 
% plot the mean of the realigned functional image time series (1 point)
r_func.mean(4).plot('rotate90', 1);

% plot the snr of the realigned functional image time series (1 point)
r_func.snr(4).plot('rotate90', 1);

% plot realignment parameters with title (1 point translation, 1 point rotation)
% Increase maximum Y ranges
figure;
hold on;
plot(realignment_parameters(:, 1:3));
title('Translational realignment parameters of image');
xlabel('Image number (n)')
ylabel('Translational Displacement (mm)')
legend('x', 'y', 'z');

figure;
hold on;
plot(realignment_parameters(:, 4:6));
title('Rotational realignment parameters of image');
xlabel('Image number (n)')
ylabel('Rotational Displacement (Deg)')
legend('pitch', 'roll', 'yaw');


%% Coregistration 
% coregister anatomy to mean of *realigned* functional (1 point)
[coreg_anat, affineCoregistrationGeometry] = anat_3T.coregister_to(r_func.mean, 'applyTransformation', 'geometry');
% save coregistered anatomy
coreg_anat.parameters.save.fileName = 'coreg_T1.nii';
coreg_anat.parameters.save.path = 'Patient 23 - 3T'; % change to your subject folder here
coreg_anat.save();

%% Quality check
% visualise aligment between realigned functional mean and coregistered anatomical
% image (1 point)
r_func.mean(4).plot('rotate90', 1); 
coreg_anat.mean(4).plot('rotate90', 1);

%% Smoothing
% apply smoothing (1 point) at twice the voxel size (1 point)
Smoothed_func = r_func.smooth('fwhm', 2*[4 4 4]);
Smoothed_func.plot('rotate90', 1);


%% Quality check
% plot the mean of the smoothed functional image time series (1 point)
Smoothed_func.mean(4).plot('rotate90', 1); 

% plot the snr of the smoothed functional image time series (1 point)
Smoothed_func.snr(4).plot('rotate90', 1);


%% GLM 
% timing (1 point)
% names (1 point)
% durations (1 point)
% onsets (1 point)

S = MrSeries();

S.parameters.save.path = './Patient 23 - 3T'; % change to your subject directory
S.data = Smoothed_func; % add your smoothed and realigned data
S.glm.regressors.realign = realignment_parameters; % add realignment parameters

% Segmented conditions with only "move"
S.glm.conditions.names = {glm_trial_type}; % add condition name
S.glm.conditions.durations = {glm_duration}; % add task durations
S.glm.conditions.onsets = {glm_onset}; % add task onsets

S.glm.timingUnits = 'secs'; % don't touch
S.glm.repetitionTime = Smoothed_func.geometry.TR_s; % add repetition time (TR)
S.glm.hrfDerivatives = [0 0]; % don't touch
% estimate GLM
S.specify_and_estimate_1st_level(); % don't touch

%%
% Saskia Notes:
%
% "Don't follow tutorial"
% Find good cluster level threshold
% 

%% Display Results
% interactive during in-person presentation
% statistical threshold (1 point)
% 0.001

% contrast (1 point)
% For Move: [1 0]
% For Resting: [0 1]

% design matrix (1 point)

% General Linear Model

% glass brain (1 point)
% coregistered anatomy underlay (1 point)

%% Grading
% total: 22 points

