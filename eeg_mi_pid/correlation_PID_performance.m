clear
clc
cd 'D:\Entrainment__\final_functions_behavior\'

%% PARAMETERS
flows = [0.5]
fhighs = [10]
channels = "all" % "max", "cluster", "$nameOfTheChannel"
start_time = 0.5
pc = 1
N_SUBJ = 22
bad_chs=[22,28,32,41,46];
num_PID_difficulty = 2

% BEHAVIOR
meanRT_subjects_easy = load("results\meanRT_subjects_easy").meanRT_subjects_easy;
meanRT_subjects_hard = load("results\meanRT_subjects_hard").meanRT_subjects_hard;
accuracy_subjects_easy = load("results\accuracy_subjects_easy").accuracy_subjects_easy;
accuracy_subjects_hard = load("results\accuracy_subjects_hard").accuracy_subjects_hard;
lengths_easyTrials_subjects = load("results\lengths_easyTrials_subjects.mat").lengths_easyTrials_subjects;
lengths_hardTrials_subjects = load("results\lengths_hardTrials_subjects.mat").lengths_hardTrials_subjects;
%% CORRELATION
for freq=1:length(flows)
    flow = flows(freq)
    fhigh = fhighs(freq)
    %% FILE LOADING
    %clusters_red_subjects_difficulty = (load("results\PID\statisticsCB\PC1\stat_red_difficulty_offset500_" + string(flow) + "to" + string(fhigh) + "hz.mat").stat_red_difficulty.posclusterslabelmat ~= 0);
    %clusters_unq1_subjects_difficulty = (load("results\PID\statisticsCB\PC1\stat_unq1_difficulty_offset500_" + string(flow) + "to" + string(fhigh) + "hz.mat").stat_unq1_difficulty.posclusterslabelmat ~= 0);
    clusters_unq2_subjects_difficulty = (load("results\PID\statisticsCB\PC1\stat_unq2_difficulty_offset500_" + string(flow) + "to" + string(fhigh) + "hz.mat").stat_unq2_difficulty.posclusterslabelmat ~= 0);
    %clusters_syn_subjects_difficulty = (load("results\PID\statisticsCB\PC1\stat_syn_difficulty_offset500_" + string(flow) + "to" + string(fhigh) + "hz.mat").stat_syn_difficulty.posclusterslabelmat ~= 0);


    % PID
    fname_folder_PID = "results\PID\";
    fname_redundancy = "REDUNDANCY_";
    fname_unq1 = "UNIQUE1_";
    fname_unq2 = "UNIQUE2_";
    fname_synergy = "SYNERGY_";
    fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_";
    fname_difficulty = "difficultySplit_";
    fname_bandp = "trials_concat_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs";
    
    fname_red_difficulty = fname_folder_PID + fname_redundancy + fname_offset + fname_difficulty + string(num_PID_difficulty) + "_" + fname_bandp;
    fname_unq1_difficulty = fname_folder_PID + fname_unq1 + fname_offset + fname_difficulty + string(num_PID_difficulty) + "_" + fname_bandp;
    fname_unq2_difficulty = fname_folder_PID + fname_unq2 + fname_offset + fname_difficulty + string(num_PID_difficulty) + "_" + fname_bandp;
    fname_syn_difficulty = fname_folder_PID + fname_synergy + fname_offset + fname_difficulty + string(num_PID_difficulty) + "_" + fname_bandp;

    lagged_red_subjects_difficulty = load(fname_red_difficulty).f;
    lagged_unq1_subjects_difficulty = load(fname_unq1_difficulty).f;
    lagged_unq2_subjects_difficulty = load(fname_unq2_difficulty).f;
    lagged_syn_subjects_difficulty = load(fname_syn_difficulty).f;

    %% DELETE BAD CHS
    for sub=1:N_SUBJ
        lagged_red_subjects_difficulty{sub}{pc, 1}(bad_chs, :) = [];
        lagged_unq1_subjects_difficulty{sub}{pc, 1}(bad_chs, :) = [];
        lagged_unq2_subjects_difficulty{sub}{pc, 1}(bad_chs, :) = [];
        lagged_syn_subjects_difficulty{sub}{pc, 1}(bad_chs, :) = [];
    end
    
    %% COMPUTATION OF THE DIFFERENCE
    difference_accuracy_subjects = accuracy_subjects_easy - accuracy_subjects_hard;
    difference_meanRT_subjects = meanRT_subjects_easy - meanRT_subjects_hard;
    difference_combinedRTAccuracy = meanRT_subjects_easy./accuracy_subjects_easy - meanRT_subjects_hard./accuracy_subjects_hard;
    combinationRTAccuracy = [meanRT_subjects_easy./(accuracy_subjects_easy) meanRT_subjects_hard./(accuracy_subjects_hard)];
    difference_lengths = lengths_hardTrials_subjects - lengths_easyTrials_subjects;

    difference_red_subjects = zeros(1, N_SUBJ);
    difference_unq1_subjects = zeros(1, N_SUBJ);
    difference_unq2_subjects = zeros(1, N_SUBJ);
    difference_syn_subjects = zeros(1, N_SUBJ);
    unq2_subjects = zeros(1, 2*N_SUBJ);
    if(channels == "all")
        for sub=1:N_SUBJ
            difference_red_subjects(sub) = mean(lagged_red_subjects_difficulty{sub}{pc, 1}(:, 1), 1) - mean(lagged_red_subjects_difficulty{sub}{pc, 1}(:, 2), 1);
            difference_unq1_subjects(sub) = mean(lagged_unq1_subjects_difficulty{sub}{pc, 1}(:, 1), 1) - mean(lagged_unq1_subjects_difficulty{sub}{pc, 1}(:, 2), 1);
            difference_unq2_subjects(sub) = mean(lagged_unq2_subjects_difficulty{sub}{pc, 1}(:, 1), 1) - mean(lagged_unq2_subjects_difficulty{sub}{pc, 1}(:, 2), 1);
            unq2_subjects(sub) = mean(lagged_unq2_subjects_difficulty{sub}{pc, 1}(:, 1), 1);
            unq2_subjects(N_SUBJ + sub) = mean(lagged_unq2_subjects_difficulty{sub}{pc, 1}(:, 2), 1);
            difference_syn_subjects(sub) = mean(lagged_syn_subjects_difficulty{sub}{pc, 1}(:, 1), 1) - mean(lagged_syn_subjects_difficulty{sub}{pc, 1}(:, 2), 1);
        end
    elseif(channels == "cluster")
        for sub=1:N_SUBJ
            difference_unq2_subjects(sub) = mean(lagged_unq2_subjects_difficulty{sub}{pc, 1}(clusters_unq2_subjects_difficulty, 1), 1) - mean(lagged_unq2_subjects_difficulty{sub}{pc, 1}(clusters_unq2_subjects_difficulty, 2), 1);
        end
    end
     unq2_subjects([13]) = [];
     accuracies = [accuracy_subjects_easy,accuracy_subjects_hard];
     accuracies([13]) = [];
%     difference_unq2_subjects(13) = [];
%     difference_combinedRTAccuracy(13) = [];

    %% VARIABLES PREPARATION FOR SCATTER PLOTS
    Cov_combinedRTAccuracyUnq2 = cov([difference_unq2_subjects'*10000, difference_combinedRTAccuracy']);
    Cov_combinedRTAccuracyUnq2_split = cov([unq2_subjects'*10000, combinationRTAccuracy']);
    Cov_combinedRTAccuracyUnq2_split_easy = cov([unq2_subjects(1:N_SUBJ)'*10000, combinationRTAccuracy(1:N_SUBJ)']);
    Cov_combinedRTAccuracyUnq2_split_hard = cov([unq2_subjects(N_SUBJ:2*N_SUBJ)'*10000, combinationRTAccuracy(N_SUBJ:2*N_SUBJ)']);
    Cov_accuracyUnq2 = cov([difference_unq2_subjects'*10000, difference_accuracy_subjects']);
    Cov_meanRTUnq2 = cov([difference_unq2_subjects'*10000, difference_meanRT_subjects']);
    [V_meanRTUnq2, D_meanRTUnq2] = eig(Cov_meanRTUnq2);
    [V_accuracyUnq2, D_accuracyUnq2] = eig(Cov_accuracyUnq2);
    [V_combinedRTAccuracyUnq2_split, D_combinedRTAccuracyUnq2_split] = eig(Cov_combinedRTAccuracyUnq2_split);
    [V_combinedRTAccuracyUnq2_split_easy, D_combinedRTAccuracyUnq2_split_easy] = eig(Cov_combinedRTAccuracyUnq2_split_easy);
    [V_combinedRTAccuracyUnq2_split_hard, D_combinedRTAccuracyUnq2_split_hard] = eig(Cov_combinedRTAccuracyUnq2_split_hard);
    [V_combinedRTAccuracyUnq2, D_combinedRTAccuracyUnq2] = eig(Cov_combinedRTAccuracyUnq2);
    avg_meanRTUnq2 = [mean(difference_unq2_subjects*10000), mean(difference_meanRT_subjects)];
    avg_accuracyUnq2 = [mean(difference_unq2_subjects*10000), mean(difference_accuracy_subjects)];
    avg_combinedRTAccuracyUnq2 = [mean(difference_unq2_subjects*10000), mean(difference_combinedRTAccuracy)];
    avg_combinedRTAccuracyUnq2_split = [mean(unq2_subjects*10000), mean(combinationRTAccuracy)];
    avg_combinedRTAccuracyUnq2_split_easy = [mean(unq2_subjects(1:N_SUBJ)*10000), mean(combinationRTAccuracy(1:N_SUBJ))];
    avg_combinedRTAccuracyUnq2_split_hard = [mean(unq2_subjects(N_SUBJ:N_SUBJ*2)*10000), mean(combinationRTAccuracy(N_SUBJ:N_SUBJ*2))];
    figure
    scatter(unq2_subjects(1:N_SUBJ)*10000, combinationRTAccuracy(1:N_SUBJ), 'red');
    hold on
    scatter(unq2_subjects(N_SUBJ:N_SUBJ*2)*10000, combinationRTAccuracy(N_SUBJ:N_SUBJ*2), 'blue');
    line([avg_combinedRTAccuracyUnq2_split_easy(1) - V_combinedRTAccuracyUnq2_split_easy(1, 2)*2; avg_combinedRTAccuracyUnq2_split_easy(1) + V_combinedRTAccuracyUnq2_split_easy(1, 2)*2], [avg_combinedRTAccuracyUnq2_split_easy(2) - V_combinedRTAccuracyUnq2_split_easy(2, 2)*2; avg_combinedRTAccuracyUnq2_split_easy(2) + V_combinedRTAccuracyUnq2_split_easy(2, 2)*2], "color", "red");
    line([avg_combinedRTAccuracyUnq2_split_hard(1) - V_combinedRTAccuracyUnq2_split_hard(1, 2)*2; avg_combinedRTAccuracyUnq2_split_hard(1) + V_combinedRTAccuracyUnq2_split_hard(1, 2)*2], [avg_combinedRTAccuracyUnq2_split_hard(2) - V_combinedRTAccuracyUnq2_split_hard(2, 2)*2; avg_combinedRTAccuracyUnq2_split_hard(2) + V_combinedRTAccuracyUnq2_split_hard(2, 2)*2], "color", "blue");
    hold off
    axis equal
    figure
    scatter(difference_unq2_subjects * 10000, difference_combinedRTAccuracy);
    hold on
    line([avg_combinedRTAccuracyUnq2(1) - V_combinedRTAccuracyUnq2(1, 2)*2; avg_combinedRTAccuracyUnq2(1) + V_combinedRTAccuracyUnq2(1, 2)*2], [avg_combinedRTAccuracyUnq2(2) - V_combinedRTAccuracyUnq2(2, 2)*2; avg_combinedRTAccuracyUnq2(2) + V_combinedRTAccuracyUnq2(2, 2)*2]);
    hold off
    axis equal
    title("Scatter " + channels + "Unique PC1 - RT-Accuracy " + flow + " to " + fhigh)
    xlabel("Difference Unique PC1")
    ylabel("Difference RT-Accuracy")
    savefig("figures\PID\scatter_RTAccuracyCombined_Unq2_" + channels + "chs_" + string(flow) + "to" + string(fhigh) + "hz.fig");
    figure
    scatter(difference_unq2_subjects * 10000, difference_accuracy_subjects);
    hold on
    line([avg_accuracyUnq2(1) - V_accuracyUnq2(1, 2); avg_accuracyUnq2(1) + V_accuracyUnq2(1, 2)], [avg_accuracyUnq2(2) - V_accuracyUnq2(2, 2); avg_accuracyUnq2(2) + V_accuracyUnq2(2, 2)]);
    hold off
    axis equal
    title("Scatter " + channels + "Unique PC1 - Accuracy " + flow + " to " + fhigh)
    xlabel("Difference Unique PC1")
    ylabel("Difference Accuracy")
    figure
    scatter(difference_unq2_subjects * 10000, difference_meanRT_subjects);
    hold on
    line([avg_meanRTUnq2(1) - V_meanRTUnq2(1, 2); avg_meanRTUnq2(1) + V_meanRTUnq2(1, 2)], [avg_meanRTUnq2(2) - V_meanRTUnq2(2, 2); avg_meanRTUnq2(2) + V_meanRTUnq2(2, 2)]);
    hold off
    axis equal
    title("Scatter " + channels + "Unique PC1 - RT " + flow + " to " + fhigh)
    xlabel("Difference Unique PC1")
    ylabel("Difference RT")
end