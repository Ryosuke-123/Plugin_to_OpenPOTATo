%% Parameter
SamplingTime = 0.1;
SamplingFreq = 10;
AllTrialData = 6000;

%% Read Sample fNIRS Data
BrainData = readmatrix('sample.csv');

%% Pre-Processing
% 1. Butter Worth Low Pass Filter[2Hz]
CutOff = 2;
order  = 3;
[b, a] = butter(order, CutOff/(SamplingFreq/2), 'low');

Lowpass_BrainData = filtfilt(b,a,BrainData);

% 2. Motion Artifact Detection
% 2-1. Motion Artifact detection from the gradient signal obtained by the Sobel Filter.
Flag_Gradient_MotionArtifactData = zeros(AllTrialData,1);
SobelKernel = [-1, 0, 1];

Gradient_Lowpass_BrainData = conv(Lowpass_BrainData,SobelKernel,'same');

Sorted_GradientData = sort(Gradient_Lowpass_BrainData);
Q1Index = floor(1 * length(Sorted_GradientData) / 4);
Q2Index = floor(2 * length(Sorted_GradientData) / 4);
Q3Index = floor(3 * length(Sorted_GradientData) / 4);
Q1  = Sorted_GradientData(Q1Index);
Q2  = Sorted_GradientData(Q2Index);
Q3  = Sorted_GradientData(Q3Index);
IQR = Q3 - Q1;
Outlier1 = Q1 - (1.5 * IQR);
Outlier2 = Q3 + (1.5 * IQR);

for allTrialData = 1:AllTrialData
    if Gradient_Lowpass_BrainData(allTrialData) > Outlier2
        Flag_Gradient_MotionArtifactData(allTrialData) = 1;
    elseif Gradient_Lowpass_BrainData(allTrialData) < Outlier1
        Flag_Gradient_MotionArtifactData(allTrialData) = 1;
    else
        Flag_Gradient_MotionArtifactData(allTrialData) = 0;
    end
end

% 2-2. Motion Artifact is detected from the standard deviation of the movement of the heartbeat (Window Size 1 s).
SD_Lowpass_BrainData = zeros(AllTrialData,1);
Flag_SD_MotionArtifactData = zeros(AllTrialData,1);

window_size = 10;

for allTrialData = 1:AllTrialData
    if allTrialData < AllTrialData - window_size + 1
        % Extract data within Window Size.
        window_data = Lowpass_BrainData(allTrialData:allTrialData + window_size - 1);

        % Calculation of standard deviation
        SD_Lowpass_BrainData(allTrialData,1) = std(window_data);
    end
end

Sorted_SDData(1:AllTrialData - window_size + 1) = sort(SD_Lowpass_BrainData(1:AllTrialData - window_size + 1,1));

Q1Index = floor(1 * length(Sorted_SDData(1:AllTrialData - window_size + 1)) / 4);
Q2Index = floor(2 * length(Sorted_SDData(1:AllTrialData - window_size + 1)) / 4);
Q3Index = floor(3 * length(Sorted_SDData(1:AllTrialData - window_size + 1)) / 4);
Q1 = Sorted_SDData(Q1Index);
Q2 = Sorted_SDData(Q2Index);
Q3 = Sorted_SDData(Q3Index);
IQR = Q3 - Q1;
Outlier1 = Q1 - (1.5 * IQR);
Outlier2 = Q3 - (1.5 * IQR);

for allTrialData = 1:AllTrialData
    if allTrialData < AllTrialData - window_size + 1
        if SD_Lowpass_BrainData(allTrialData) > Outlier2
            Flag_SD_MotionArtifactData(allTrialData:allTrialData + window_size - 1) = 1;
        elseif SD_Lowpass_BrainData(allTrialData) < Outlier1
            Flag_SD_MotionArtifactData(allTrialData:allTrialData + window_size - 1) = 1;
        end
    end
end

% 4. Baseline Shift Detection
% To detect the maximum amplitude change due to a heartbeat, the data sequence is shifted 0.5[s] (approximately half of the heartbeat cycle) forward
Shift_Number = 5;

Shift_Lowpass_BrainData = circshift(Lowpass_BrainData,-1 * Shift_Number);

% Amplitude change (absolute value) due to heartbeat calculated from original data and Shift data
HeartBeat_Change = abs(Shift_Lowpass_BrainData - Lowpass_BrainData);

% Calculate the maximum value in the "amplitude change (absolute value) due to heartbeat" during the Motion Artifact-free period, and use it as the threshold for detecting Baseline Shift. This is used as the threshold for detecting Baseline Shift.
MAx_HeartBeat_Change = zeros(1,1);

for allTrialData = 1:AllTrialData 
    if allTrialData < AllTrialData - BackCutData - window_size + 1
        % No Motion Artifact
        if Flag_Gradient_MotionArtifactData(allTrialData) == 0 && Flag_Gradient_MotionArtifactData(allTrialData + 1) == 0 && Flag_Gradient_MotionArtifactData(allTrialData + 2) == 0 && Flag_Gradient_MotionArtifactData(allTrialData + 3) == 0 && Flag_Gradient_MotionArtifactData(allTrialData + 4) == 0
            if HeartBeat_Change(allTrialData) > Max_HeartBeat_Change
                Max_HeartBeat_Change = HeartBeat_Change(allTrialData);
            end
        end
    end
end

% Determines whether the Motion Artifact is a Baseline Shift or not based on the calculated threshold value.
Flag_BaselineShiftData = zeros(AllTrialData,1);

for allTrialData = 1:AllTrialData
    if allTrialData < AllTrialData - window_size + 1
        % Motion Artifact
        if  Flag_Gradient_MotionArtifactData(allTrialData) == 1
            if HeartBeat_Change(allTrialData) > Max_HeartBeat_Change
                Flag_BaselineShiftData(allTrialData:allTrialData + Shift_Number - 1) = 1;
            end
        end
    end
end

% 5. Baseline Shift Detection
% Modelling with Spline interpolation for periods of Baseline Shift
Spline_BrainData = zeros(AllTrialData,1);

for allTrialData = 1:AllTrialData
    if allTrialData < AllTrialData - window_size + 1
        % Baseline Shift
        if Flag_BaselineShiftData(allTrialData) == 1
            Spline_BrainData(allTrialData - 1) = Lowpass_BrainData(allTrialData - 1);
            Spline_BrainData(allTrialData)     = Lowpass_BrainData(allTrialData);
        end
    end
end
        
SplineCorr_BrainData = zeros(AllTrialData,1);
smoothness = 0.99;

count = 0;
dataPoints = [];

for allTrialData = 1:AllTrialData
    if allTrialData < AllTrialData - Shift_Number
        if Spline_BrainData(allTrialData) ~= 0
            count = count + 1;
            dataPoints(count) = Spline_BrainData(allTrialData);
        elseif Spline_BrainData(allTrialData) == 0 && count > 0
            if count >= 2
                SplineCorr_BrainData(allTrialData - count:allTrialData - 1) = csaps(1:count,dataPoints,smoothness,1:count);
            end
            count = 0;
            dataPoints = [];
        end
    end
end

% Subtracts the modelled spline interpolation signal from the original signal, removing the Baseline Shift component
BaselineShift_Corr_BrainData = zeros(AllTrialData,1);

for allTrialData = 1:AllTrialData
    if allTrialData < AllTrialData - Shift_Number
        if Flag_BaselineShiftData(allTrialData) == 1
            BaselineShift_Corr_BrainData(allTrialData - 1) = Lowpass_BrainData(allTrialData - 1) - SplineCorr_BrainData(allTrialData - 1);
            BaselineShift_Corr_BrainData(allTrialData)     = Lowpass_BrainData(allTrialData) - SplineCorr_BrainData(allTrialData);
        end
    end
end

% Replaces data with Baseline Shift components with data with Baseline Shift components removed and corrects vertical shift
Replace_BaselineShift_BrainData = zeros(AllTrialData,1);
Corr_BaselineShift_BrainData    = zeros(AllTrialData,1);

alpha = round(SamplingFreq / 3);
beta  = 2 * SamplingFreq;
mean_a = 0;
mean_b = 0;
theta  = 0; 

flag_first_segment = 1;

pre_count  = 1;
post_count = 1;
pre_segment       = [];
post_segment      = [];
corr_post_segment = [];

% 5-1. Replaces data with Baseline Shift components with data with Baseline Shift components removed.
for allTrialData = 1:AllTrialData
    if allTrialData < AllTrialData - Shift_Number
        % Baseline Shift
        if Flag_BaselineShiftData(allTrialData) == 0
            Replace_BaselineShift_BrainData(allTrialData) = Lowpass_BrainData(allTrialData);
        elseif Flag_BaselineShiftData(allTrialData) == 1
            Replace_BaselineShift_BrainData(allTrialData) = BaselineShift_Corr_BrainData(allTrialData);
        end
    end
end

% 5-2. Correction of vertical shift
for allTrialData = 1:AllTrialData
    if allTrialData < AllTrialData - Shift_Number
        % No BS Segment
        if Flag_BaselineShiftData(allTrialData) == 0
            % Yes 1st Segment
            if flag_first_segment == 1
                Corr_BaselineShift_BrainData(allTrialData) = Replace_BaselineShift_BrainData(allTrialData);

                pre_segment(pre_count) = Replace_BaselineShift_BrainData(allTrialData);
                pre_count = pre_count + 1;
            % No 1st Segment
            elseif flag_first_segment == 0
                post_segment(post_count) = Replace_BaselineShift_BrainData(allTrialData);
                % post_count = post_count + 1;

                % Next Sampling data, Yes BS Segment
                if Flag_BaselineShiftData(allTrialData + 1) == 1

                    % Correction of vertical shift
                    % Condition 1
                    if length(pre_segment) <= alpha && length(post_segment) <= alpha
                        sum_a = 0;
                        sum_b = 0;

                        for i = 1:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:length(post_segment)
                            sum_b = sum_b + pre_segment(i);
                        end

                        mean_a = sum_a / length(pre_segment);
                        mean_b = sum_b / length(post_segment);
                    % Condition 2
                    elseif length(pre_segment) <= alpha && length(post_segment) > alpha && length(post_segment) < beta
                        sum_a = 0;
                        sum_b = 0;

                        for i = 1:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:alpha
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / length(pre_segment);
                        mean_b = sum_b / alpha;

                        theta = mean_a - mean_b;
                    % Condition 3
                    elseif length(pre_segment) <= alpha && length(post_segment) >= beta
                        sum_a = 0;
                        sum_b = 0;
                        theta_post = ceil(length(post_segment) / 10);

                        for i = 1:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:theta_post
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / length(pre_segment);
                        mean_b = sum_b / theta_post;

                        theta = mean_a - mean_b;
                   % Condition 4
                    elseif length(pre_segment) > alpha && length(pre_segment) < beta && length(post_segment) <= alpha
                        sum_a = 0;
                        sum_b = 0;

                        for i = length(pre_segment) - alpha:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:length(post_segment)
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / (alpha + 1);
                        mean_b = sum_b / length(post_segment);

                        theta = mean_a - mean_b;
                    % Condition 5
                    elseif length(pre_segment) > alpha && length(pre_segment) < beta && length(post_segment) > alpha && length(post_segment) < beta
                        sum_a = 0;
                        sum_b = 0;

                        for i = length(pre_segment) - alpha:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:alpha
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / (alpha + 1);
                        mean_b = sum_b / alpha;

                        theta = mean_a - mean_b;
                    % Condition 6
                    elseif length(pre_segment) > alpha && length(pre_segment) < beta && length(post_segment) >= beta
                        sum_a = 0;
                        sum_b = 0;
                        theta_post = ceil(length(post_segment) / 10);

                        for i = length(pre_segment) - alpha:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:theta_post
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / (alpha + 1);
                        mean_b = sum_b / theta_post;

                        theta = mean_a - mean_b;
                    % Condition 7
                    elseif length(pre_segment) >= beta && length(post_segment) <= alpha
                        sum_a = 0;
                        sum_b = 0;
                        theta_pre = ceil(length(pre_segment) / 10);

                        for i = length(pre_segment) - theta_pre:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:length(post_segment)
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / (theta_pre + 1);
                        mean_b = sum_b / length(post_segment);

                        theta = mean_a - mean_b;
                    % Condition 8
                    elseif length(pre_segment) >= beta && length(post_segment) > alpha && length(post_segment) < beta
                        sum_a = 0;
                        sum_b = 0;
                        theta_pre = ceil(length(pre_segment) / 10);

                        for i = length(pre_segment) - theta_pre:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:alpha
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / (theta_pre + 1);
                        mean_b = sum_b / alpha;

                        theta = mean_a - mean_b;
                    % Condition 9
                    elseif length(pre_segment) >= beta && length(post_segment) >= beta
                        sum_a = 0;
                        sum_b = 0;
                        theta_pre  = ceil(length(pre_segment) / 10);
                        theta_post = ceil(length(post_segment) / 10);

                        for i = length(pre_segment) - theta_pre:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:theta_post
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / (theta_pre + 1);
                        mean_b = sum_b / theta_post;

                        theta = mean_a - mean_b;
                    end

                    % Correction of vertical shift
                    Corr_BaselineShift_BrainData(allTrialData - post_count + 1:allTrialData) = Replace_BaselineShift_BrainData(allTrialData - post_count + 1:allTrialData) + theta;
                    corr_post_segment(1:post_count) = Corr_BaselineShift_BrainData(allTrialData - post_count + 1:allTrialData);

                    % Initialise PreSegment and move the contents of Corr_PostSegment
                    pre_segment  = [];
                    post_segment = [];
                    post_count   = 1;
                    pre_segment = corr_post_segment;
                    corr_post_segment = [];

                elseif allTrialData + 1 == AllTrialData - Shift_Number
                    % Correct vertical shift
                    % Condition 1
                    if length(pre_segment) <= alpha && length(post_segment) <= alpha
                        sum_a = 0;
                        sum_b = 0;

                        for i = 1:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:length(post_segment)
                            sum_b = sum_b + pre_segment(i);
                        end

                        mean_a = sum_a / length(pre_segment);
                        mean_b = sum_b / length(post_segment);

                        theta = mean_a - mean_b;
                    % Condition 2
                    elseif length(pre_segment) <= alpha && length(post_segment) > alpha && length(post_segment) < beta
                        sum_a = 0;
                        sum_b = 0;

                        for i = 1:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:alpha
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / length(pre_segment);
                        mean_b = sum_b / alpha;

                        theta = mean_a - mean_b;
                    % Condition 3
                    elseif length(pre_segment) <= alpha && length(post_segment) >= beta
                        sum_a = 0;
                        sum_b = 0;
                        theta_post = ceil(length(post_segment) / 10);

                        for i = 1:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:theta_post
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / length(pre_segment);
                        mean_b = sum_b / theta_post;

                        theta = mean_a - mean_b;
                    % Condition 4
                    elseif length(pre_segment) > alpha && length(pre_segment) < beta && length(post_segment) <= alpha
                        sum_a = 0;
                        sum_b = 0;

                        for i = length(pre_segment) - alpha:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:length(post_segment)
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / (alpha + 1);
                        mean_b = sum_b / length(post_segment);

                        theta = mean_a - mean_b;
                    % Condition 5
                    elseif length(pre_segment) > alpha && length(pre_segment) < beta && length(post_segment) > alpha && length(post_segment) < beta
                        sum_a = 0;
                        sum_b = 0;

                        for i = length(pre_segment) - alpha:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:alpha
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / (alpha + 1);
                        mean_b = sum_b / alpha;

                        theta = mean_a - mean_b;
                    % Condition 6
                    elseif length(pre_segment) > alpha && length(pre_segment) < beta && length(post_segment) >= beta
                        sum_a = 0;
                        sum_b = 0;
                        theta_post = ceil(length(post_segment) / 10);

                        for i = length(pre_segment) - alpha:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:theta_post
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / (alpha + 1);
                        mean_b = sum_b / theta_post;

                        theta = mean_a - mean_b;
                    % Condition 7
                    elseif length(pre_segment) >= beta && length(post_segment) <= alpha
                        sum_a = 0;
                        sum_b = 0;
                        theta_pre = ceil(length(pre_segment) / 10);

                        for i = length(pre_segment) - theta_pre:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:length(post_segment)
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / (theta_pre + 1);
                        mean_b = sum_b / length(post_segment);

                        theta = mean_a - mean_b;                                        
                    % Condition 8
                    elseif length(pre_segment) >= beta && length(post_segment) > alpha && length(post_segment) < beta
                        sum_a = 0;
                        sum_b = 0;
                        theta_pre = ceil(length(pre_segment) / 10);

                        for i = length(pre_segment) - theta_pre:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:alpha
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / (theta_pre + 1);
                        mean_b = sum_b / alpha;

                        theta = mean_a - mean_b;                                        
                    % Condition 9
                    elseif length(pre_segment) >= beta && length(post_segment) >= beta
                        sum_a = 0;
                        sum_b = 0;
                        theta_pre  = ceil(length(pre_segment) / 10);
                        theta_post = ceil(length(post_segment) / 10);

                        for i = length(pre_segment) - theta_pre:length(pre_segment)
                            sum_a = sum_a + pre_segment(i);
                        end
                        for i = 1:theta_post
                            sum_b = sum_b + post_segment(i);
                        end

                        mean_a = sum_a / (theta_pre + 1);
                        mean_b = sum_b / theta_post;

                        theta = mean_a - mean_b;
                    end

                    % Correction of vertical shift
                    Corr_BaselineShift_BrainData(allTrialData - post_count + 1:allTrialData) = Replace_BaselineShift_BrainData(allTrialData - post_count + 1:allTrialData) + theta;
                    corr_post_segment(1:post_count) = Corr_BaselineShift_BrainData(allTrialData - post_count + 1:allTrialData);

                    % Initialise PreSegment and move the contents of Corr_PostSegment
                    pre_segment  = [];
                    post_segment = [];
                    post_count   = 1;
                    pre_segment  = corr_post_segment;
                    corr_post_segment = [];
                else
                    post_count = post_count + 1;
                end
            end

        % Yes BS Segment
        elseif Flag_BaselineShiftData(allTrialData) == 1
            % Lower the flag to determine the first Segment
            flag_first_segment = 0;
            pre_count = 1;

            post_segment(post_count) = Replace_BaselineShift_BrainData(allTrialData);
            % post_count = post_count + 1;

            % Next Sampling data, No BS Segment
            if Flag_BaselineShiftData(allTrialData + 1) == 0
                % Correct vertical shift
                % Condition 1
                if length(pre_segment) <= alpha && length(post_segment) <= alpha
                    sum_a = 0;
                    sum_b = 0;

                    for i = 1:length(pre_segment)
                        sum_a = sum_a + pre_segment(i);
                    end
                    for i = 1:length(post_segment)
                        sum_b = sum_b + pre_segment(i);
                    end

                    mean_a = sum_a / length(pre_segment);
                    mean_b = sum_b / length(post_segment);

                    theta = mean_a - mean_b;
                % Condition 2
                elseif length(pre_segment) <= alpha && length(post_segment) > alpha && length(post_segment) < beta
                    sum_a = 0;
                    sum_b = 0;

                    for i = 1:length(pre_segment)
                        sum_a = sum_a + pre_segment(i);
                    end
                    for i = 1:alpha
                        sum_b = sum_b + post_segment(i);
                    end

                    mean_a = sum_a / length(pre_segment);
                    mean_b = sum_b / alpha;

                    theta = mean_a - mean_b;
                % Condition 3
                elseif length(pre_segment) <= alpha && length(post_segment) >= beta
                    sum_a = 0;
                    sum_b = 0;
                    theta_post = ceil(length(post_segment) / 10);

                    for i = 1:length(pre_segment)
                        sum_a = sum_a + pre_segment(i);
                    end
                    for i = 1:theta_post
                        sum_b = sum_b + post_segment(i);
                    end

                    mean_a = sum_a / length(pre_segment);
                    mean_b = sum_b / theta_post;

                    theta = mean_a - mean_b;
                % Condition 4
                elseif length(pre_segment) > alpha && length(pre_segment) < beta && length(post_segment) <= alpha
                    sum_a = 0;
                    sum_b = 0;

                    for i = length(pre_segment) - alpha:length(pre_segment)
                        sum_a = sum_a + pre_segment(i);
                    end
                    for i = 1:length(post_segment)
                        sum_b = sum_b + post_segment(i);
                    end

                    mean_a = sum_a / (alpha + 1);
                    mean_b = sum_b / length(post_segment);

                    theta = mean_a - mean_b;
                % Condition 5
                elseif length(pre_segment) > alpha && length(pre_segment) < beta && length(post_segment) > alpha && length(post_segment) < beta
                    sum_a = 0;
                    sum_b = 0;

                    for i = length(pre_segment) - alpha:length(pre_segment)
                        sum_a = sum_a + pre_segment(i);
                    end
                    for i = 1:alpha
                        sum_b = sum_b + post_segment(i);
                    end

                    mean_a = sum_a / (alpha + 1);
                    mean_b = sum_b / alpha;

                    theta = mean_a - mean_b;
                % Condition 6
                elseif length(pre_segment) > alpha && length(pre_segment) < beta && length(post_segment) >= beta
                    sum_a = 0;
                    sum_b = 0;
                    theta_post = ceil(length(post_segment) / 10);

                    for i = length(pre_segment) - alpha:length(pre_segment)
                        sum_a = sum_a + pre_segment(i);
                    end
                    for i = 1:theta_post
                        sum_b = sum_b + post_segment(i);
                    end

                    mean_a = sum_a / (alpha + 1);
                    mean_b = sum_b / theta_post;

                    theta = mean_a - mean_b;
                % Condition 7
                elseif length(pre_segment) >= beta && length(post_segment) <= alpha
                    sum_a = 0;
                    sum_b = 0;
                    theta_pre = ceil(length(pre_segment) / 10);

                    for i = length(pre_segment) - theta_pre:length(pre_segment)
                        sum_a = sum_a + pre_segment(i);
                    end
                    for i = 1:length(post_segment)
                        sum_b = sum_b + post_segment(i);
                    end

                    mean_a = sum_a / (theta_pre + 1);
                    mean_b = sum_b / length(post_segment);

                    theta = mean_a - mean_b;
                % Condition 8
                elseif length(pre_segment) >= beta && length(post_segment) > alpha && length(post_segment) < beta
                    sum_a = 0;
                    sum_b = 0;
                    theta_pre = ceil(length(pre_segment) / 10);

                    for i = length(pre_segment) - theta_pre:length(pre_segment)
                        sum_a = sum_a + pre_segment(i);
                    end
                    for i = 1:alpha
                        sum_b = sum_b + post_segment(i);
                    end

                    mean_a = sum_a / (theta_pre + 1);
                    mean_b = sum_b / alpha;

                    theta = mean_a - mean_b;
                % Condition 9
                elseif length(pre_segment) >= beta && length(post_segment) >= beta
                    sum_a = 0;
                    sum_b = 0;
                    theta_pre  = ceil(length(pre_segment) / 10);
                    theta_post = ceil(length(post_segment) / 10);

                    for i = length(pre_segment) - theta_pre:length(pre_segment)
                        sum_a = sum_a + pre_segment(i);
                    end
                    for i = 1:theta_post
                        sum_b = sum_b + post_segment(i);
                    end

                    mean_a = sum_a / (theta_pre + 1);
                    mean_b = sum_b / theta_post;

                    theta = mean_a - mean_b;
                end

                % Correction of vertical shift
                Corr_BaselineShift_BrainData(allTrialData - post_count + 1:allTrialData) = Replace_BaselineShift_BrainData(allTrialData - post_count + 1:allTrialData) + theta;
                corr_post_segment(1:post_count) = Corr_BaselineShift_BrainData(allTrialData - post_count + 1:allTrialData);

                % Initialise PreSegment and move the contents of Corr_PostSegment
                pre_segment  = [];
                post_segment = [];
                post_count   = 1;
                pre_segment  = corr_post_segment;
                corr_post_segment = [];
            else
                post_count = post_count + 1;
            end
        end
    end
end

figure(1);
hold on;
plot(Lowpass_BrainData,'black','LineWidth',1.5);
plot(Flag_BaselineShiftData,'blue','LineWidth',1.5);
plot(Corr_BaselineShift_BrainData,'green','LineWidth',1.5);
legend('Raw Brain Data','Baseline Shift Flag','Corr Brain Data');
xlabel('Time');
ylabel('ΔHb');
