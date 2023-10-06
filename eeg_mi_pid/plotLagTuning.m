clear
clc
% close all

addpath 'results/PID/'
load 'chs_eeg_ok'
load('freq_template')
pc= 1;
flow = 5;
fhigh = 7;
start_time = 0.5;
bad_chs=[22,28,32,41,46];
freq_fourier_all.label = cellstr(chs_eeg_ok);
lags = [-0.2:0.05:0.4];
NUM_SUBJ = 22;
NUM_CHS = 59;
Rdn_subs_lags_easy = zeros(length(lags),NUM_SUBJ,NUM_CHS);
Rdn_subs_lags_hard = zeros(length(lags),NUM_SUBJ,NUM_CHS);
Syn_subs_lags_easy = zeros(length(lags),NUM_SUBJ,NUM_CHS);
Syn_subs_lags_hard = zeros(length(lags),NUM_SUBJ,NUM_CHS);
Unq1_subs_lags_easy = zeros(length(lags),NUM_SUBJ,NUM_CHS);
Unq1_subs_lags_hard = zeros(length(lags),NUM_SUBJ,NUM_CHS);
Unq2_subs_lags_easy = zeros(length(lags),NUM_SUBJ,NUM_CHS);
Unq2_subs_lags_hard = zeros(length(lags),NUM_SUBJ,NUM_CHS);

for n=1:length(lags)
    lag = lags(n);
    disp("**flow: " + flow + ", fhigh: " + fhigh + "**");
    fname_redundancy = "REDUNDANCY_";
    fname_unq1 = "UNIQUE1_";
    fname_unq2 = "UNIQUE2_";
    fname_synergy = "SYNERGY_";
    fname_offset = "offset" + string(start_time * 1000) + "ms_lag" + string(lag*1000) + "ms_allFreqs_";
    fname_difficulty = "difficultySplit_";
    fname_bandp = "trials_concat_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";
    Rdn_subs = load("results/PID/" + fname_redundancy + fname_offset + fname_difficulty + fname_bandp).Rdn;
    Unq1_subs = load("results/PID/" + fname_unq1 + fname_offset + fname_difficulty + fname_bandp).Unq1;
    Unq2_subs = load("results/PID/" + fname_unq2 + fname_offset + fname_difficulty + fname_bandp).Unq2;
    Syn_subs = load("results/PID/" + fname_synergy + fname_offset + fname_difficulty + fname_bandp).Syn;
    
    for sub=1:NUM_SUBJ
        % Rdn
        Rdn_subs{sub}{pc,1}(bad_chs,:) = [];
        Rdn_subs_lags_easy(n,sub,:) = Rdn_subs{sub}{pc,1}(:,1);
        Rdn_subs_lags_hard(n,sub,:) = Rdn_subs{sub}{pc,1}(:,2);
        % Unq1
        Unq1_subs{sub}{pc,1}(bad_chs,:) = [];
        Unq1_subs_lags_easy(n,sub,:) = Unq1_subs{sub}{pc,1}(:,1);
        Unq1_subs_lags_hard(n,sub,:) = Unq1_subs{sub}{pc,1}(:,2);
        % Unq2
        Unq2_subs{sub}{pc,1}(bad_chs,:) = [];
        Unq2_subs_lags_easy(n,sub,:) = Unq2_subs{sub}{pc,1}(:,1);
        Unq2_subs_lags_hard(n,sub,:) = Unq2_subs{sub}{pc,1}(:,2);
        % Syn
        Syn_subs{sub}{pc,1}(bad_chs,:) = [];
        Syn_subs_lags_easy(n,sub,:) = Syn_subs{sub}{pc,1}(:,1);
        Syn_subs_lags_hard(n,sub,:) = Syn_subs{sub}{pc,1}(:,2);
    end
end

%% ALL SUBJ
figure
subplot(2,2,1)
shadedErrorBar(lags, mean(mean(Rdn_subs_lags_easy,3),2), std(mean(Rdn_subs_lags_easy,3),0,2)/sqrt(NUM_SUBJ),'lineProps', '-r')
hold on
shadedErrorBar(lags, mean(mean(Rdn_subs_lags_hard,3),2), std(mean(Rdn_subs_lags_hard,3),0,2)/sqrt(NUM_SUBJ),'lineProps', '-b')
legend("EASY", "HARD");
title("Rdn All")
subplot(2,2,2)
shadedErrorBar(lags, mean(mean(Unq1_subs_lags_easy,3),2), std(mean(Unq1_subs_lags_easy,3),0,2)/sqrt(NUM_SUBJ),'lineProps', '-r')
hold on
shadedErrorBar(lags, mean(mean(Unq1_subs_lags_hard,3),2), std(mean(Unq1_subs_lags_hard,3),0,2)/sqrt(NUM_SUBJ),'lineProps', '-b')
legend("EASY", "HARD");
title("SE All")
subplot(2,2,3)
shadedErrorBar(lags, mean(mean(Unq2_subs_lags_easy,3),2), std(mean(Unq2_subs_lags_easy,3),0,2)/sqrt(NUM_SUBJ),'lineProps', '-r')
hold on
shadedErrorBar(lags, mean(mean(Unq2_subs_lags_hard,3),2), std(mean(Unq2_subs_lags_hard,3),0,2)/sqrt(NUM_SUBJ),'lineProps', '-b')
legend("EASY", "HARD");
title("PC1 All")
subplot(2,2,4)
shadedErrorBar(lags, mean(mean(Syn_subs_lags_easy,3),2), std(mean(Syn_subs_lags_easy,3),0,2)/sqrt(NUM_SUBJ),'lineProps', '-r')
hold on
shadedErrorBar(lags, mean(mean(Syn_subs_lags_hard,3),2), std(mean(Syn_subs_lags_hard,3),0,2)/sqrt(NUM_SUBJ),'lineProps', '-b')
legend("EASY", "HARD");
title("Syn All")

%% SINGLE SUBJ
for sub=[1:10]
    figure
    subplot(2,2,1)
    shadedErrorBar(lags, mean(squeeze(Rdn_subs_lags_easy(:,sub,:)),2), std(squeeze(Rdn_subs_lags_easy(:,sub,:)),0,2)/sqrt(NUM_CHS),'lineProps', '-r')
    hold on
    shadedErrorBar(lags, mean(squeeze(Rdn_subs_lags_hard(:,sub,:)),2), std(squeeze(Rdn_subs_lags_hard(:,sub,:)),0,2)/sqrt(NUM_CHS),'lineProps', '-b')
    legend("EASY", "HARD");
    title("Rdn")

    subplot(2,2,2)
    shadedErrorBar(lags, mean(squeeze(Unq1_subs_lags_easy(:,sub,:)),2), std(squeeze(Unq1_subs_lags_easy(:,sub,:)),0,2)/sqrt(NUM_CHS),'lineProps', '-r')
    hold on
    shadedErrorBar(lags, mean(squeeze(Unq1_subs_lags_hard(:,sub,:)),2), std(squeeze(Unq1_subs_lags_hard(:,sub,:)),0,2)/sqrt(NUM_CHS),'lineProps', '-b')
    title("UnqSE")

    subplot(2,2,3)
    shadedErrorBar(lags, mean(squeeze(Unq2_subs_lags_easy(:,sub,:)),2), std(squeeze(Unq2_subs_lags_easy(:,sub,:)),0,2)/sqrt(NUM_CHS),'lineProps', '-r')
    hold on
    shadedErrorBar(lags, mean(squeeze(Unq2_subs_lags_hard(:,sub,:)),2), std(squeeze(Unq2_subs_lags_hard(:,sub,:)),0,2)/sqrt(NUM_CHS),'lineProps', '-b')
    title("UnqPC1")

    subplot(2,2,4)
    shadedErrorBar(lags, mean(squeeze(Syn_subs_lags_easy(:,sub,:)),2), std(squeeze(Syn_subs_lags_easy(:,sub,:)),0,2)/sqrt(NUM_CHS),'lineProps', '-r')
    hold on
    shadedErrorBar(lags, mean(squeeze(Syn_subs_lags_hard(:,sub,:)),2), std(squeeze(Syn_subs_lags_hard(:,sub,:)),0,2)/sqrt(NUM_CHS),'lineProps', '-b')
    title("Syn")
end