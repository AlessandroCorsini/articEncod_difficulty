clear
clc
% cd 'D:\Entrainment__\final_functions_behavior\'

%% ADDING PATHS AND SETTING PARAMS
addpath '../fieldtrip-20210311'
addpath './results'
addpath '../analisi'
addpath(pwd)

%%% PARAMETERS %%%
PCs = [1];
subs = 1:22;
NUM_RANDOMIZATION = 1000;
flow=0.5;
fhigh=10;
TRIAL_CONCAT = true;
PLOT_DIFFICULTY = true;
PLOT_PERFORMANCE = false;
N_SUBJ = 22;
start_time = 0.5;
bad_chs=[22,28,32,41,46];
load 'chs_eeg_ok'
chs_eeg=chs_eeg_ok
load('freq_template')

for pc=PCs
    %% FILE LOADING
    fname_redundancy = "REDUNDANCY_";
    fname_unq1 = "UNIQUE1_";
    fname_unq2 = "UNIQUE2_";
    fname_synergy = "SYNERGY_";
    fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_";
    fname_difficulty = "difficultySplit_";
    fname_performance = "performanceSplit_";
    fname_bandp = "trials_concat_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";
    fname_matchLength ="matchedLength_";
    if(PLOT_DIFFICULTY)
        rdn_subjects_difficulty = load("results/PID/" + fname_redundancy + fname_offset + fname_difficulty + fname_bandp).Rdn;
        unq1_subjects_difficulty = load("results/PID/" + fname_unq1 + fname_offset + fname_difficulty + fname_bandp).Unq1;
        unq2_subjects_difficulty = load("results/PID/" + fname_unq2 + fname_offset + fname_difficulty + fname_bandp).Unq2;
        syn_subjects_difficulty = load("results/PID/" + fname_synergy + fname_offset + fname_difficulty + fname_bandp).Syn;
        rdn_subjects_difficulty_matchedLength = load("results/PID/" + fname_redundancy + fname_offset + fname_difficulty + fname_matchLength + fname_bandp).Rdn;
        unq1_subjects_difficulty_matchedLength = load("results/PID/" + fname_unq1 + fname_offset + fname_difficulty + fname_matchLength + fname_bandp).Unq1;
        unq2_subjects_difficulty_matchedLength = load("results/PID/" + fname_unq2 + fname_offset + fname_difficulty + fname_matchLength + fname_bandp).Unq2;
        syn_subjects_difficulty_matchedLength = load("results/PID/" + fname_synergy + fname_offset + fname_difficulty + fname_matchLength + fname_bandp).Syn;
    end

    for sub=subs
        % Mean original PID
        rdn_subjects_difficulty{sub}{pc,1}(bad_chs, :) = [];
        unq1_subjects_difficulty{sub}{pc,1}(bad_chs, :) = [];
        unq2_subjects_difficulty{sub}{pc,1}(bad_chs, :) = [];
        syn_subjects_difficulty{sub}{pc,1}(bad_chs, :) = [];

        Rdn_sub_mean_allchs = mean(rdn_subjects_difficulty{sub}{pc,1});
        Unq1_sub_mean_allchs = mean(unq1_subjects_difficulty{sub}{pc,1});
        Unq2_sub_mean_allchs = mean(unq2_subjects_difficulty{sub}{pc,1});
        Syn_sub_mean_allchs = mean(syn_subjects_difficulty{sub}{pc,1});

        % Mean matchedLength PID
        Rdn_sub_mean_allchs_matchedLength_allrandomizations = zeros(1000, 2);
        Unq1_sub_mean_allchs_matchedLength_allrandomizations = zeros(1000, 2);
        Unq2_sub_mean_allchs_matchedLength_allrandomizations = zeros(1000, 2);
        Syn_sub_mean_allchs_matchedLength_allrandomizations = zeros(1000, 2);
        for i=1:NUM_RANDOMIZATION
            rdn_subjects_difficulty_matchedLength{sub}{1, i}{pc,1}(bad_chs, :) = [];
            unq1_subjects_difficulty_matchedLength{sub}{1, i}{pc,1}(bad_chs, :) = [];
            unq2_subjects_difficulty_matchedLength{sub}{1, i}{pc,1}(bad_chs, :) = [];
            syn_subjects_difficulty_matchedLength{sub}{1, i}{pc,1}(bad_chs, :) = [];
            Rdn_sub_mean_allchs_matchedLength_allrandomizations(i, :) = mean(rdn_subjects_difficulty_matchedLength{sub}{1,i}{pc,1});
            Unq1_sub_mean_allchs_matchedLength_allrandomizations(i, :) = mean(unq1_subjects_difficulty_matchedLength{sub}{1,i}{pc,1});
            Unq2_sub_mean_allchs_matchedLength_allrandomizations(i, :) = mean(unq2_subjects_difficulty_matchedLength{sub}{1,i}{pc,1});
            Syn_sub_mean_allchs_matchedLength_allrandomizations(i, :) = mean(syn_subjects_difficulty_matchedLength{sub}{1,i}{pc,1});
        end
        Rdn_sub_mean_allchs_matchedLength = mean(Rdn_sub_mean_allchs_matchedLength_allrandomizations);
        Rdn_sub_std_allchs_matchedLength = std(Rdn_sub_mean_allchs_matchedLength_allrandomizations);
        Unq1_sub_mean_allchs_matchedLength = mean(Unq1_sub_mean_allchs_matchedLength_allrandomizations);
        Unq1_sub_std_allchs_matchedLength = std(Unq1_sub_mean_allchs_matchedLength_allrandomizations);
        Unq2_sub_mean_allchs_matchedLength = mean(Unq2_sub_mean_allchs_matchedLength_allrandomizations);
        Unq2_sub_std_allchs_matchedLength = std(Unq2_sub_mean_allchs_matchedLength_allrandomizations);
        Syn_sub_mean_allchs_matchedLength = mean(Syn_sub_mean_allchs_matchedLength_allrandomizations);
        Syn_sub_std_allchs_matchedLength = std(Syn_sub_mean_allchs_matchedLength_allrandomizations);
        
        %% CHECK PERCENTILES
%         p_5 = prctile(Unq2_sub_mean_allchs_matchedLength_allrandomizations(:,2), 10);
%         p_95 = prctile(Unq2_sub_mean_allchs_matchedLength_allrandomizations(:,2), 90);
%         if(Unq2_sub_mean_allchs(2) > p_5 && Unq2_sub_mean_allchs(2) < p_95)
%             disp("Unq2 sub: " + sub + " ok");
%         end
%         p_5 = prctile(Unq1_sub_mean_allchs_matchedLength_allrandomizations(:,2), 10);
%         p_95 = prctile(Unq1_sub_mean_allchs_matchedLength_allrandomizations(:,2), 90);
%         if(Unq1_sub_mean_allchs(2) > p_5 && Unq1_sub_mean_allchs(2) < p_95)
%             disp("Unq1 sub: " + sub + " ok");
%         end
%         p_5 = prctile(Rdn_sub_mean_allchs_matchedLength_allrandomizations(:,2), 10);
%         p_95 = prctile(Rdn_sub_mean_allchs_matchedLength_allrandomizations(:,2), 90);
%         if(Rdn_sub_mean_allchs(2) > p_5 && Rdn_sub_mean_allchs(2) < p_95)
%             disp("Rdn sub: " + sub + " ok");
%         end
%         p_5 = prctile(Syn_sub_mean_allchs_matchedLength_allrandomizations(:,2), 10);
%         p_95 = prctile(Syn_sub_mean_allchs_matchedLength_allrandomizations(:,2), 90);
%         if(Syn_sub_mean_allchs(2) > p_5 && Syn_sub_mean_allchs(2) < p_95)
%             disp("Syn sub: " + sub + " ok");
%         end

%         PLOT sub by sub
%         figure
%         histfit(Rdn_sub_mean_allchs_matchedLength_allrandomizations(:,end));
% 
%         figure
%         histfit(Unq1_sub_mean_allchs_matchedLength_allrandomizations(:,end));

        figure
        H = histfit(Unq2_sub_mean_allchs_matchedLength_allrandomizations(:,end));
        H(1).FaceColor = [0.65,0.65,0.65];
        H(2).Color = [0,0,0];
        hold on
        xline(Unq2_sub_mean_allchs(2), 'Color', 'red', 'LineWidth', 1);
%         saveas(gcf, "Unq2_1000_PID_matchedLength_sub_" + sub + ".pdf"); 

%         figure
%         histfit(Syn_sub_mean_allchs_matchedLength_allrandomizations(:,end));

        %% TTEST
%         [h, p] = ttest(Rdn_sub_mean_allchs_matchedLength_allrandomizations(:,2), Rdn_sub_mean_allchs(2), 'Alpha', 0.05);
%         if h==0 
%             disp("Rdn sub: " + sub + ", p= " + p);
%         end
%         [h, p] = ttest(Unq1_sub_mean_allchs_matchedLength_allrandomizations(:,2), Unq1_sub_mean_allchs(2), 'Alpha', 0.05);
%         if h==0
%             disp("Unq1 sub: " + sub + ", p= " + p);
%         end
%         [h, p] = ttest(Unq2_sub_mean_allchs_matchedLength_allrandomizations(:,2), Unq2_sub_mean_allchs(2), 'Alpha', 0.05);
%         if h==0 
%             disp("Unq2 sub: " + sub + ", p= " + p);
%         end
%         [h, p] = ttest(Syn_sub_mean_allchs_matchedLength_allrandomizations(:,2), Syn_sub_mean_allchs(2), 'Alpha', 0.05);
%         if h==0 
%             disp("Syn sub: " + sub + ", p= " + p);
%         end

    end
end
