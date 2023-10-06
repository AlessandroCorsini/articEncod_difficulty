clear
clc
% cd 'D:/Entrainment__/final_functions_behavior/'

%% ADDING PATHS
addpath '../partial-info-decomp-master'
addpath '../fieldtrip-20210311'
addpath '../gcmi-master/matlab'
addpath '../analisi'

%% PARAMETERS
% random number genarator for PID computation
rng(2)
REALIGN_INPUTS=true;
REALIGN_EEG=true;
% Compute PID on concatenated trials splitted by difficulty
COMPUTE_PID_DIFFICULTY = true;
% Compute PID on concatenated trials splitted by performance
COMPUTE_PID_PERFORMANCE = false;
% Compute PID on concatenated trials splitted by difficulty or performance
% matching the length of the two concatenations
MATCH_LENGHT = true;
% The subjects used in the computation of PID
SELECTED_SUB = [1:22];
% If true the concatenation of the trials is done on matched length trials.
% So some samples are cutted from each trial before being concatenated.
CONCATENATE_SINGLE_TRIAL_MATCHED_LENGTH = false;
% if CONCATENATE_SINGLE_TRIAL_MATCHED_LENGTH, it sets the length of the
% single trial to be matched. 'max' mean the match is done on the shortest
% single trial. So it is the max to take all the trials. Otherwise indicate
% how many samplesGoing up in the length implies the exclusion of some trials.
% Right now only implemented for lengths <= 'max'.
SINGLE_TRIAL_LENGTH = 'mean'; %CHANGE ALSO IN THE PARFOR LOOP
% Number of iteration of PID. In combination with match_length = true at
% every computation a different segment is cutted for matching
PID_ITERATIONS = 500;
% flows = [0.5:1:8.5]
% fhighs = [2.5:1:10.5]
flows = [5];
fhighs = [7];
% The lag is between the Stimuli (speech
% envelope, PCs) and the Response (EEG channel) which is supposed to be
% delayed
Fs = 400;
lags = [0.2]; % in seconds
start_time = 0.5;
PCs = [1]; % PCs to use
pca_dir = "../pca_new_nocut";
env_dir = "../envelope";
load('cut_times')
cut_times=round(cut_times*Fs*0.001); %from ms to samples

%% LOADING SUBS NAMES
trigger_name = "startSound";
d=dir('../analisi');
subs_files ={};
subs={};
subs_n={};
for i=3:length(d)
    if endsWith(d(i).name,"_filtered_epoched" + trigger_name + "_ica_interp-epo.fif")
        subs_files=[subs_files; d(i).name];
        tmp=split(d(i).name,'_f');
        tmp2=split(tmp(1),'_');
        subs=[subs; tmp(1)];
        subs_n=[subs_n; tmp2(2)];
    end
end
subs= subs(SELECTED_SUB);
subs_files=subs_files(SELECTED_SUB);
subs_n=subs_n(SELECTED_SUB);
N_SUBJ = length(subs);

%% LOADING PCs
d_pca_em=dir(pca_dir);
pca_names={};
for i=4:length(d_pca_em)
    pca_names=[pca_names; d_pca_em(i).name];
end
pca_names=natsort(pca_names);
pcas={};
for i=1:length(pca_names);
    pca_tmp=load(pca_dir+"/"+pca_names{i});
    pcas=[pcas pca_tmp.pca_tmp];
end

%% LOADING ENVELOPE
d_env=dir(env_dir);
envelope_names={};
for i=3:length(d_env)
    envelope_names=[envelope_names; d_env(i).name];
end
envelope_names=natsort(envelope_names);
envelopes={};
for i=1:length(envelope_names)
    envelopes=[envelopes readNPY(env_dir+"/"+envelope_names{i})];
end

%% CUTTING FEW SAMPLES (20to50) AT THE END OF ENVELOPES - they are little bit longer than PCs
for i=1:length(envelopes)
    %first compute the few samples in excess in speech respect to ema(or
    %pca) due to subsampling, and eliminate them from the end of envelopes
    diff= length(envelopes{i})-length(pcas{i});
    disp("Cutting " + diff + "samples from envelope: " + envelope_names(i));
    envelopes{i}=envelopes{i}(1:end-diff);
    if REALIGN_INPUTS
        % read the time from where to cut files to eliminate silences present
        % at the beginning of sentences
        start_cutTime=cut_times(i,1);
        end_cutTime=cut_times(i,2);
        envelopes{1,i}=envelopes{1,i}(start_cutTime:end-end_cutTime);
        pcas{i}=pcas{i}(:,start_cutTime:end-end_cutTime);
    end
end

%% ANALYSIS
NUM_PERFORMANCE = length(load("../behavioral_/results/subs/" + subs{1} + "_dataSplit").dataSplit.goodPerformanceTrials);
% Performance split to use
performances = false(1, NUM_PERFORMANCE);
performances([1, 2, 3]) = true;
% CHARACTERIZATION DIFFICULTY-SPLIT
meanRT_subjects_easy = zeros(1, N_SUBJ);
meanRT_subjects_hard = zeros(1, N_SUBJ);
accuracy_subjects_easy = zeros(1, N_SUBJ);
accuracy_subjects_hard = zeros(1, N_SUBJ);
lengths_subjects_easy = zeros(1, N_SUBJ);
lengths_subjects_hard = zeros(1, N_SUBJ);
% CHARACTERIZATION PERFORMANCE-SPLITS
meanRT_subjects_good = {};
meanRT_subjects_bad = {};
accuracy_subjects_good = {};
accuracy_subjects_bad = {};
lengths_subjects_good = {};
lengths_subjects_bad = {};
disp("***START OF THE CORE PART OF THE ANALYSIS***");


%% LAG TUNING
for freq=1:length(flows)
    for lag=lags
        disp("LAG: " + lag*1000 + "ms");
        %% NAME DEFINITION
        flow = flows(freq);
        fhigh = fhighs(freq);
        disp("**flow: " + flow + ", fhigh: " + fhigh + "**");
        fname_redundancy = "REDUNDANCY_";
        fname_unq1 = "UNIQUE1_";
        fname_unq2 = "UNIQUE2_";
        fname_synergy = "SYNERGY_";
        fname_offset = "offset" + string(start_time * 1000) + "ms_lag" + string(lag*1000) + "ms_allFreqs_";
        fname_difficulty = "difficultySplit_";
        fname_performance = "performanceSplit_";
        fname_bandp = "trials_concat_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";
        if(CONCATENATE_SINGLE_TRIAL_MATCHED_LENGTH)
            fname_bandp = "trials_concat_singleTrial_matchedLength_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";
        end
        if(MATCH_LENGHT)
            fname_difficulty = fname_difficulty + "matchedLength_501_1000_";
            fname_performance = fname_performance + "matchedLength_";
        end
        Rdn_subs = cell(1,length(SELECTED_SUB));
        Unq1_subs =cell(1,length(SELECTED_SUB));
        Unq2_subs =cell(1,length(SELECTED_SUB));
        Syn_subs =cell(1,length(SELECTED_SUB));
        disp("catch");


        %% SINGLE SUBJECT ANALYSIS
        parfor n=SELECTED_SUB
            %% NAME DEFINITION
            disp("*Processing subject " + subs(n) + "*");
            % PID DIFFICULTY
            if(COMPUTE_PID_DIFFICULTY)
                Rdn={};
                Unq1={};
                Unq2={};
                Syn={};
                disp("catch");
            end

            %% METAFILES LOADING
            cfg=[];
            file=subs_files{n};
            cfg.dataset=file;
            data_epochs=ft_preprocessing(cfg)
            %load("../behavioral_/results/phrases/behaviorPhrases.mat");
            % Loading the behavioral data
            behavior = load("../behavioral_/results/subs/" + subs{n} + "_behavior").behavior;
            dataSplit = load("../behavioral_/results/subs/" + subs{n} + "_dataSplit").dataSplit;
            goodTrials_original = zeros(1, length(behavior.goodTrials));
            goodTrials_original(load("../analisi/" + subs{n} + "_goodTrials.mat").good_trials + 1) = 1;
            disp("The original good trials were: " + sum(goodTrials_original));
            disp("Now the good trials are: " + sum(behavior.goodTrials));
            validTrials_original = goodTrials_original == 1 & (behavior.answers == 0 | behavior.answers == 1);
            validTrials = behavior.goodTrials == 1 & (behavior.answers == 0 | behavior.answers == 1);

            %% ADJUSTING VECTOR OF TIME AND TRIALS USING THE SAME TIMES USED FOR ENVELOPES AND PCAs AND FIND INVALID TRIALS AND DELETING INVALID EPOCHS
            if REALIGN_EEG
                epochs_invalid = [];
                for trial=1:length(validTrials)
                    if(validTrials(trial))
                        start_cutTime = cut_times(behavior.phrases(1, trial),1);
                        end_cutTime = cut_times(behavior.phrases(1, trial),2);
                        epoch = sum(find(goodTrials_original) <= trial);
                        disp("Trial " + trial + " is valid and corresponds to epoch " + epoch);
                        t0=data_epochs.time{epoch}(1);
                        data_epochs.time{epoch}=round(data_epochs.time{epoch}(1,start_cutTime:end-end_cutTime)+ t0 - data_epochs.time{epoch}(1,start_cutTime), 5);
                        data_epochs.trial{epoch}=data_epochs.trial{epoch}(:,start_cutTime:end-end_cutTime);
                    elseif(goodTrials_original(trial) == 1)
                        epoch = sum(find(goodTrials_original) <= trial);
                        disp("Trial " + trial + "is not valid but originally good and corresponds to epoch " + epoch);
                        epochs_invalid = [epochs_invalid, epoch];
                    end
                end
            end
            data_epochs.time(epochs_invalid) = [];
            data_epochs.trial(epochs_invalid) = [];
            if(sum(validTrials) == length(data_epochs.trial))
                disp("OK. The epochs are " + length(data_epochs.trial) + ", " + length(epochs_invalid) + " were deleted");
            else
                disp("KO. There is a problem with the number of epochs (" + length(data_epochs.trial) + ")");
            end

            %% ORDERING PCAs AND ENVELOPES TRIALS (AS PROVIDED IN THE EXPERIMENT) TO MATCH THEM WITH EEG TRIALS
            all_pcas_ordered={};
            all_envelopes_ordered={};
            validPhrases = behavior.phrases(validTrials);
            for i=1:length(validPhrases)
                all_pcas_ordered=[all_pcas_ordered pcas(validPhrases(i))];
                all_envelopes_ordered=[all_envelopes_ordered envelopes(validPhrases(i))];
            end

            %% SETTING THE TIME TO CUT SEGMENTS IN ORDER TO EXCLUDE AUDITORY POTENTIALS AT THE TRIAL'S ONSET AND EXCLUSION OF CORRUPTED EMA
            corrupted_pca_trials=[];
            for i=1:length(all_pcas_ordered)
                if length(all_pcas_ordered{i}) < 2000
                    % The index in corrupted_pca_trials is referred to the index of
                    % the trial in the validTrials only
                    corrupted_pca_trials = [corrupted_pca_trials i];
                end
            end
            disp("Corrupted PC trials: " + corrupted_pca_trials);
            data_epochs.trial(corrupted_pca_trials)=[];
            data_epochs.time(corrupted_pca_trials)=[];
            all_pcas_ordered(corrupted_pca_trials)=[];
            all_envelopes_ordered(corrupted_pca_trials)=[];
            % Update valid trials. The i-th 1 in validTrials must be set to 0.
            % In this way validTrials and data_epochs are still coherent
            idxs_valid = find(validTrials);
            validTrials(idxs_valid(corrupted_pca_trials)) = 0;
            validPhrases = behavior.phrases(validTrials);
            if(COMPUTE_PID_DIFFICULTY)
                if(length(corrupted_pca_trials) == (sum(dataSplit.easyTrials(idxs_valid(corrupted_pca_trials))) + sum(dataSplit.hardTrials(idxs_valid(corrupted_pca_trials)))))
                    disp("OK. All Corrupted trials were removed from difficulty");
                else
                    disp("KO. Not all the corrupted PC were removed from difficulty");
                end
                dataSplit.easyTrials(idxs_valid(corrupted_pca_trials)) = 0;
                dataSplit.hardTrials(idxs_valid(corrupted_pca_trials)) = 0;
            end
            if(COMPUTE_PID_PERFORMANCE)
                for n_performance=1:NUM_PERFORMANCE
                    if(performances(n_performance))
                        if(length(corrupted_pca_trials) == (sum(dataSplit.goodPerformanceTrials{n_performance}(idxs_valid(corrupted_pca_trials))) + sum(dataSplit.badPerformanceTrials{n_performance}(idxs_valid(corrupted_pca_trials)))))
                            disp("OK. All Corrupted trials were removed from performance " + n_performance);
                        else
                            disp("KO. Not all the corrupted PC were removed from performance " + n_performance);
                        end
                        dataSplit.goodPerformanceTrials{n_performance}(idxs_valid(corrupted_pca_trials)) = 0;
                        dataSplit.badPerformanceTrials{n_performance}(idxs_valid(corrupted_pca_trials)) = 0;
                    end
                end
            end
            if(sum(validTrials) == length(data_epochs.trial))
                disp("OK. The epochs are " + length(data_epochs.trial) + " and are still coincident with the number of validTrials");
            else
                disp("KO. There is a problem with the number of epochs (" + length(data_epochs.trial) + ")");
            end

            %% ACCURACY,RT AND LENGTH FOR THE EASY AND FOR THE HARD AND FOR GOOD AND BAD
            if(COMPUTE_PID_DIFFICULTY)
                if(sum(dataSplit.easyTrials) == sum(dataSplit.easyTrials & validTrials') && sum(dataSplit.hardTrials) == sum(dataSplit.hardTrials & validTrials'))
                    disp("All correct");
                end
                meanRT_subjects_easy(n) = mean(behavior.RTs(dataSplit.easyTrials & validTrials'));
                meanRT_subjects_hard(n) = mean(behavior.RTs(dataSplit.hardTrials & validTrials'));
                accuracy_subjects_easy(n) = mean(behavior.answers(dataSplit.easyTrials & validTrials'));
                accuracy_subjects_hard(n) = mean(behavior.answers(dataSplit.hardTrials & validTrials'));
                %lengths_subjects_easy(n) = sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.easyTrials & validTrials')));
                %lengths_subjects_hard(n) = sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.hardTrials & validTrials')));
            end

            %% REMOVING 500ms TO EXCLUDE AUDITORY POTENTIALS AT THE TRIAL'S ONSET
            % Exclusion of stimulus-locked evoked potentials
            % (1 x length(validTrials == 1)) cell array. Cell i contains the
            % (samples x 21) double matrix referred to the PCs of the phrase number
            % validPhrases(i) after the 1.5s shift
            pcas_timeShifted={};
            % (1 x length(validTrials == 1)) cell array. Cell i contains the
            % (samples x 1) double matrix referred to the values of the envelope of
            % phrase number validPhrases(i) after the 1.5s shift
            envelopes_timeShifted={};
            lengths_timeShifted=[];
            start_samp= Fs*start_time;
            for i=1:length(all_pcas_ordered)
                pcas_timeShifted{i}=all_pcas_ordered{i}(:,start_samp:end)';
                envelopes_timeShifted{i}=all_envelopes_ordered{i}(start_samp:end);
                lengths_timeShifted=[lengths_timeShifted length(envelopes_timeShifted{i})];
            end
            T=find(data_epochs.time{1}==start_time+lag);
            cfg=[];
            cfg.endsample=[];
            cfg.begsample=[];
            for i=1:length(data_epochs.trial)
                cfg.begsample=[cfg.begsample T];
                cfg.endsample=[cfg.endsample T+lengths_timeShifted(i)-1];
            end
            % Reallignment of the EEG signals after the shift of 0.5s and the delay
            % of 0.2s given by natural time processing of the speech stimulus
            data_epochs = ft_redefinetrial(cfg, data_epochs);

            %% FILTERING INPUTS FROM flow TO fhigh
            % Input means Speech evnvelope + PCs
            data_epochs_input=data_epochs;
            data_epochs_input.elec=[];
            data_epochs_input.cfg=[];
            data_epochs_input.label=[];
            % (22 x 1) cell array with the labels of the enevelope and the PCs
            data_epochs_input.label{1,1}='Envelope';
            % Iteration over all the 21 PCs
            for i=1:size(pcas_timeShifted{1},2)
                data_epochs_input.label{i+1,1}=char('PC'+string(i));
            end
            % (1 x length(validTrials == 1)) cell array. Cell i contains a (22 x
            % samples) double matrix in which the first row is the envelope of
            % trial i, the last 21 are the PCs of trial i
            data_epochs_input.trial=[];
            for i=1:length(data_epochs.trial)
                data_epochs_input.trial{1,i}=[envelopes_timeShifted{i} pcas_timeShifted{i}]';
                data_epochs_input.time{1,i}=data_epochs_input.time{1,i}-lag;
            end
            data_epochs_input.hdr.label=data_epochs_input.label;
            data_epochs_input.hdr.nChans=length(data_epochs_input.label);
            data_epochs_input.hdr.elec=[];
            % Filtering the Inputs: envelope + PCs
            if fhigh ~=0
                cgf=[];
                cfg.bpfilter = 'yes';
                cfg.bpfreq = [flow fhigh];
                cfg.bpfiltord = 2;
                cfg.padtype ='mirror';
                cfg.bpfilttype    = "but";
                data_epochs_input_filt = ft_preprocessing(cfg, data_epochs_input);
            end
            % Assigning the filtered signals back to the variables
            for i=1:length(data_epochs_input_filt.trial)
                envelopes_timeShifted{i}=data_epochs_input_filt.trial{i}(1,:)';
                pcas_timeShifted{i}=data_epochs_input_filt.trial{i}(2:end,:)';
            end

            %% FILTERING EEG FROM flow TO fhigh
            if fhigh ~=0
                cgf=[];
                cfg.bpfilter = 'yes';
                cfg.bpfreq = [flow fhigh];
                cfg.bpfiltord = 2;
                cfg.padtype ='mirror';
                cfg.bpfilttype    = "but";
                data_epochs = ft_preprocessing(cfg, data_epochs);
            end

            %% COMPUTING PID
            data_eeg=data_epochs.trial;
            % PID DIFFICULTY
            if(COMPUTE_PID_DIFFICULTY)
                class1 = {find(dataSplit.easyTrials(validTrials == 1))};
                class2 = {find(dataSplit.hardTrials(validTrials == 1))};
            end
            % PID PERFORMANCE
            if(COMPUTE_PID_PERFORMANCE)
                class1 = {};
                class2 = {};
                for n_performance=1:NUM_PERFORMANCE
                    class1{n_performance} = find(dataSplit.goodPerformanceTrials{n_performance}(validTrials == 1));
                    class2{n_performance} = find(dataSplit.badPerformanceTrials{n_performance}(validTrials == 1));
                end
            end
            %% PID COMPUTATION
            for i=1:length(class1)
                %% PID COMPUTATION FOR ONE ENTIRE SPLIT
                for pid_iter=1:PID_ITERATIONS
                    % necessary for parfor loop
                    SINGLE_TRIAL_LENGTH = 'mean';
                    %% SINGLE PID ITERATION FOR ONE SPLIT
                    for k=PCs
                        % (length(validTrials == 1) x 1) cell array. Cell i contains the
                        % (samples x 1) double matrix representing the k-th PC.
                        pcas_timeShifted_k={};
                        min_length_epoch = length(pcas_timeShifted{1}(:, k));
                        for ep=1:length(pcas_timeShifted)
                            pcas_timeShifted_k=[pcas_timeShifted_k; pcas_timeShifted{ep}(:,k)];
                            if(length(pcas_timeShifted{ep}(:,k)) < min_length_epoch)
                                min_length_epoch = length(pcas_timeShifted{ep}(:,k));
                            end
                        end

                        %% INITIALIZATION OF TRIAL CONCATENATIONS FOR CLASS1 AND CLASS 2
                        if(pid_iter==1)
                            if(CONCATENATE_SINGLE_TRIAL_MATCHED_LENGTH && strcmp(SINGLE_TRIAL_LENGTH, 'max'))
                                SINGLE_TRIAL_LENGTH = min_length_epoch;
                            end
                            for trial=1:length(class1{i})
                                if(not(CONCATENATE_SINGLE_TRIAL_MATCHED_LENGTH))
                                    SINGLE_TRIAL_LENGTH = length(data_eeg{class1{i}(trial)});
                                end
                                if trial==1
                                    eeg_class1=data_eeg{class1{i}(trial)}(:, 1:SINGLE_TRIAL_LENGTH);
                                    env_class1=envelopes_timeShifted{class1{i}(trial)}(1:SINGLE_TRIAL_LENGTH);
                                    pc_class1=pcas_timeShifted_k{class1{i}(trial)}(1:SINGLE_TRIAL_LENGTH);
                                else
                                    eeg_class1=[eeg_class1 data_eeg{class1{i}(trial)}(:, 1:SINGLE_TRIAL_LENGTH)];
                                    env_class1=[env_class1; envelopes_timeShifted{class1{i}(trial)}(1:SINGLE_TRIAL_LENGTH)];
                                    pc_class1=[pc_class1; pcas_timeShifted_k{class1{i}(trial)}(1:SINGLE_TRIAL_LENGTH)];
                                end
                            end

                            for trial=1:length(class2{i})
                                if(not(CONCATENATE_SINGLE_TRIAL_MATCHED_LENGTH))
                                    SINGLE_TRIAL_LENGTH = length(data_eeg{class2{i}(trial)});
                                end
                                if trial==1
                                    eeg_class2=data_eeg{class2{i}(trial)}(:, 1:SINGLE_TRIAL_LENGTH);
                                    env_class2=envelopes_timeShifted{class2{i}(trial)}(1:SINGLE_TRIAL_LENGTH);
                                    pc_class2=pcas_timeShifted_k{class2{i}(trial)}(1:SINGLE_TRIAL_LENGTH);
                                else
                                    eeg_class2=[eeg_class2 data_eeg{class2{i}(trial)}(:, 1:SINGLE_TRIAL_LENGTH)];
                                    env_class2=[env_class2; envelopes_timeShifted{class2{i}(trial)}(1:SINGLE_TRIAL_LENGTH)];
                                    pc_class2=[pc_class2; pcas_timeShifted_k{class2{i}(trial)}(1:SINGLE_TRIAL_LENGTH)];
                                end
                            end
                        end
                        if(MATCH_LENGHT)
                            %% ONE SINGLE MATCHED_LENGTH COMPUTATION
                            eeg_class1_matched = eeg_class1;
                            eeg_class2_matched=eeg_class2;
                            env_class1_matched=env_class1;
                            env_class2_matched=env_class2;
                            pc_class1_matched=pc_class1;
                            pc_class2_matched=pc_class2;
                            rng(pid_iter)
                            diff = length(eeg_class2) - length(eeg_class1);
                            if(diff > 0)
                                start_sample = floor(rand(1)*length(eeg_class1));
                                eeg_class2_matched(:, start_sample+1:start_sample+diff) = [];
                                env_class2_matched(start_sample+1:start_sample+diff) = [];
                                pc_class2_matched(start_sample+1:start_sample+diff) = [];
                            end
                            if(diff < 0)
                                start_sample = rand(1)*length(eeg_class2);
                                eeg_class1_matched(:, start_sample+1:start_sample+diff) = [];
                                env_class1_matched(start_sample+1:start_sample+diff) = [];
                                pc_class1_matched(start_sample+1:start_sample+diff) = [];
                            end
                            disp("ITERATION: " + pid_iter);
                            disp("SAMPLES DELETED: (" + start_sample + ", " + (start_sample+diff) + ")");
                            disp("eeg_class1: " + length(eeg_class1_matched) + " samples");
                            disp("env_class1: " + length(env_class1_matched) + " samples");
                            disp("pc_class1: " + length(pc_class1_matched) + " samples");
                            disp("eeg_class2: " + length(eeg_class2_matched) + " samples");
                            disp("env_class2: " + length(env_class2_matched) + " samples");
                            disp("pc_class2: " + length(pc_class2_matched) + " samples");
                            [red_all_ch_class1,unq1_all_ch_class1,unq2_all_ch_class1,syn_all_ch_class1] = PID_copnorm(eeg_class1_matched,env_class1_matched,pc_class1_matched);
                            [red_all_ch_class2,unq1_all_ch_class2,unq2_all_ch_class2,syn_all_ch_class2] = PID_copnorm(eeg_class2_matched,env_class2_matched,pc_class2_matched);
                        else
                            %% ONE SINGLE NON MATCHED_LENGTH COMPUTATION
                            [red_all_ch_class1,unq1_all_ch_class1,unq2_all_ch_class1,syn_all_ch_class1] = PID_copnorm(eeg_class1,env_class1,pc_class1);
                            [red_all_ch_class2,unq1_all_ch_class2,unq2_all_ch_class2,syn_all_ch_class2] = PID_copnorm(eeg_class2,env_class2,pc_class2);
                        end
                        Rdn{k, 1} = [red_all_ch_class1, red_all_ch_class2];
                        Unq1{k, 1} = [unq1_all_ch_class1, unq1_all_ch_class2];
                        Unq2{k, 1} = [unq2_all_ch_class1, unq2_all_ch_class2];
                        Syn{k, 1} = [syn_all_ch_class1, syn_all_ch_class2];
                    end
                    Rdn_subs{n}{pid_iter} = Rdn;
                    Unq1_subs{n}{pid_iter} = Unq1;
                    Unq2_subs{n}{pid_iter} = Unq2;
                    Syn_subs{n}{pid_iter} = Syn;
                end
            end
        end
        cd 'results/PID/'
        %% SAVING RESULTS
        % LENGTH OF TRIALS CONCAT
        %     if(COMPUTE_PID_DIFFICULTY)
        %         % MEAN RT AND ACCURACY FOR EASY AND HARD
        %         save('..\meanRT_subjects_easy.mat', "meanRT_subjects_easy");
        %         save('..\meanRT_subjects_hard.mat', "meanRT_subjects_hard");
        %         save('..\accuracy_subjects_easy', "accuracy_subjects_easy");
        %         save('..\accuracy_subjects_hard', "accuracy_subjects_hard");
        %         save('..\lengths_subjects_easy', "lengths_subjects_easy");
        %         save('..\lengths_subjects_hard', "lengths_subjects_hard");
        %     end
        %     if(COMPUTE_PID_PERFORMANCE)
        %         % MEAN RT AND ACCURACY FOR GOOD AND BAD
        %         save('..\meanRT_subjects_good.mat', "meanRT_subjects_good");
        %         save('..\meanRT_subjects_bad.mat', "meanRT_subjects_bad");
        %         save('..\accuracy_subjects_good', "accuracy_subjects_good");
        %         save('..\accuracy_subjects_bad', "accuracy_subjects_bad");
        %         save('..\lengths_subjects_good', "lengths_subjects_good");
        %         save('..\lengths_subjects_bad', "lengths_subjects_bad");
        %     end

        % DIFFICULTY
        if(COMPUTE_PID_DIFFICULTY)
            % {1} Is related to the fact that there is just one split for
            % difficulty
            Rdn = Rdn_subs;
            fname = fname_redundancy + fname_offset + fname_difficulty + fname_bandp;
            save(fname, "Rdn");
            Unq1 = Unq1_subs;
            fname = fname_unq1 + fname_offset + fname_difficulty + fname_bandp;
            save(fname, "Unq1");
            Unq2 = Unq2_subs;
            fname = fname_unq2 + fname_offset + fname_difficulty + fname_bandp;
            save(fname, "Unq2");
            Syn = Syn_subs;
            fname = fname_synergy + fname_offset + fname_difficulty + fname_bandp;
            save(fname, "Syn");
        end

        % PERFORMANCE
        if(COMPUTE_PID_PERFORMANCE)
            for n_performance=1:NUM_PERFORMANCE
                if(performances(n_performance))
                    Rdn = Rdn_subs{n_performance};
                    fname = fname_redundancy + fname_offset + fname_performance + string(n_performance) + "_" + fname_bandp;
                    save(fname, "Rdn");
                    Unq1 = Unq1_subs{n_performance};
                    fname = fname_unq1 + fname_offset + fname_performance + string(n_performance) + "_" + fname_bandp;
                    save(fname, "Unq1");
                    Unq2 = Unq2_subs{n_performance};
                    fname = fname_unq2 + fname_offset + fname_performance + string(n_performance) + "_" + fname_bandp;
                    save(fname, "Unq2");
                    Syn = Syn_subs{n_performance};
                    fname = fname_synergy + fname_offset + fname_performance + string(n_performance) + "_" + fname_bandp;
                    save(fname, "Syn");
                end
            end
        end
        cd '../../'
    end
end

