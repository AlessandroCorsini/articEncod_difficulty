clear
clc
% cd 'D:\Entrainment__\final_functions_behavior\'

%% ADDING PATHS AND SETTING PARAMS
rng(1000)

addpath '../partial-info-decomp-master'
addpath '../fieldtrip-20210311'
addpath '../gcmi-master/matlab'
addpath '../analisi'

realign_inputs=true
realign_eeg=true
computePID_difficulty = true
computePID_performance = false
matchLength = true
% flows = [linspace(0.5, 9.5, 19)]; %linspace(0.5, 10, 20)];
% fhighs = [linspace(1.5, 10.5, 19)]; %linspace(1, 10.5, 20)];
flows = [5]
fhighs = [7]
% The lag is between the Stimuli (speech
% envelope, PCs) and the Response (EEG channel) which is supposed to be
% delayed
Fs=400;
delay=0.2;%s
lag=80; %delay in sample -> i.e. 200ms at 400hz
load('cut_times')
cut_times=round(cut_times*Fs*0.001); %from ms to samples

PCs = [1]; % PCs to use

%% FOLDERS

pca_dir="../pca_new_nocut";
env_dir="../envelope";
% cd results/PID/
% load('REDUNDANCY_offset500ms_allFreqs_difficultySplit_singleTrial_matchedLength_eegBandp_5to7hz_filteredInputs.mat');
% load('UNIQUE1_offset500ms_allFreqs_difficultySplit_singleTrial_matchedLength_eegBandp_5to7hz_filteredInputs.mat');
% 
% load('UNIQUE2_offset500ms_allFreqs_difficultySplit_singleTrial_matchedLength_eegBandp_5to7hz_filteredInputs.mat');
% 
% load('SYNERGY_offset500ms_allFreqs_difficultySplit_singleTrial_matchedLength_eegBandp_5to7hz_filteredInputs.mat');
% 
% 
% cd '../../'

%% LOADING DATA NAMES EEG

trigger_name = "startSound"

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


%% LOADING PCA NAMES AND DATA

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

%% LOADING ENVELOPE NAMES AND DATA

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

%% CUTTING FEW SAMPLES (20to50) AT THE END OF ENVELOPES - they are little bit longer than emas


for i=1:length(envelopes)

    %first compute the few samples in excess in speech respect to ema(or
    %pca) due to subsampling, and eliminate them from the end of envelopes
    diff= length(envelopes{i})-length(pcas{i});
    disp("Cutting " + diff + "samples from envelope: " + envelope_names(i));
    envelopes{i}=envelopes{i}(1:end-diff);

    if realign_inputs

        % read the time from where to cut files to eliminate silences present
        % at the beginning of sentences
        start_cutTime=cut_times(i,1);
        end_cutTime=cut_times(i,2);

        envelopes{1,i}=envelopes{1,i}(start_cutTime:end-end_cutTime);
        pcas{i}=pcas{i}(:,start_cutTime:end-end_cutTime);
    end

end


%% LOADING META DATA AND EEG FOR EACH SUBJECT n - CORRECT TRIALS
NUM_PERFORMANCE = length(load("../behavioral_/results/subs/" + subs{1} + "_dataSplit.mat").dataSplit.goodPerformanceTrials);
% Performance split to use
performances = false(1, NUM_PERFORMANCE);
performances([2, 3]) = true;

disp("***START OF THE CORE PART OF THE ANALYSIS***");
for freq=1:length(flows)

    % PID DIFFICULTY
    lagged_redundancy_subjects_difficulty={};
    lagged_unq1_subjects_difficulty={};
    lagged_unq2_subjects_difficulty={};
    lagged_synergy_subjects_difficulty={};

    % PID PERFORMANCE
    lagged_redundancy_subjects_performance={};
    lagged_unq1_subjects_performance={};
    lagged_unq2_subjects_performance={};
    lagged_synergy_subjects_performance={};

    flow = flows(freq);
    fhigh = fhighs(freq);
    disp("**flow: " + flow + ", fhigh: " + fhigh + "**");
    for n=1:length(subs_files)
        disp("*Processing subject " + subs(n) + "*");

        % PID DIFFICULTY
        lagged_redundancy_difficulty={};
        lagged_unq1_difficulty={};
        lagged_unq2_difficulty={};
        lagged_synergy_difficulty={};

        % PID PERFORMANCE
        lagged_redundancy_performance={};
        lagged_unq1_performance={};
        lagged_unq2_performance={};
        lagged_synergy_performance={};

        cfg=[];
        file=subs_files{n};
        cfg.dataset=file;
        data_epochs=ft_preprocessing(cfg)
        load("../behavioral_/results/phrases/behaviorPhrases.mat");
        % Loading the behavioral data
        behavior = load("../behavioral_/results/subs/" + subs{n} + "_behavior.mat").behavior;
        % Remember to write the code in order to manage more than one split:
        % each split in one column of the matrix
        dataSplit = load("../behavioral_/results/subs/" + subs{n} + "_dataSplit.mat").dataSplit;
        goodTrials_original = zeros(1, length(behavior.goodTrials));
        goodTrials_original(load("../analisi/" + subs{n} + "_goodTrials.mat").good_trials + 1) = 1;
        disp("The original good trials were: " + sum(goodTrials_original));
        disp("Now the good trials are: " + sum(behavior.goodTrials));

        validTrials_original = goodTrials_original == 1 & (behavior.answers == 0 | behavior.answers == 1);
        validTrials = behavior.goodTrials == 1 & (behavior.answers == 0 | behavior.answers == 1);

        % ADJUSTING VECTOR OF TIME AND TRIALS USING THE SAME TIMES USED FOR ENVELOPES AND PCAs
        if realign_eeg
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

        %REMOVING TRIALS WHERE THE ANSWER TO THE ATTENTIVE QUESTION WAS
        %INCORRECT

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

        %% SETTING THE TIME TO CUT SEGMENTS IN ORDER TO EXCLUDE AUDITORY POTENTIALS AT THE TRIAL'S ONSET

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

        % TO BE CUTTED
% for pc=1
%         Rdn_concat_easy = {};
%         Rdn_concat_hard = {};
%         Unq1_concat_easy = {};
%         Unq1_concat_hard = {};
%         Unq2_concat_easy = {};
%         Unq2_concat_hard = {};
%         Syn_concat_easy = {};
%         Syn_concat_hard = {};
% for i=1:length(Rdn{1}{1,n})
% Rdn_concat_easy(i) = Rdn{1}{1,n}{1,i}(pc,1);
% Unq1_concat_easy(i) = Unq1{1}{1,n}{1,i}(pc,1);
% Unq2_concat_easy(i) = Unq2{1}{1,n}{1,i}(pc,1);
% Syn_concat_easy(i) = Syn{1}{1,n}{1,i}(pc,1);
% end
% for i=1:length(Rdn{2}{1,n})
% Rdn_concat_hard(i) = Rdn{2}{1,n}{1,i}(pc,1);
% Unq1_concat_hard(i) = Unq1{2}{1,n}{1,i}(pc,1);
% Unq2_concat_hard(i) = Unq2{2}{1,n}{1,i}(pc,1);
% Syn_concat_hard(i) = Syn{2}{1,n}{1,i}(pc,1);
% end
% 
% Rdn_new{n}{pc, 1}(dataSplit.easyTrials & validTrials') = Rdn_concat_easy;
% Rdn_new{n}{pc, 1}(dataSplit.hardTrials & validTrials') = Rdn_concat_hard;
% Unq1_new{n}{pc, 1}(dataSplit.easyTrials & validTrials') = Unq1_concat_easy;
% Unq1_new{n}{pc, 1}(dataSplit.hardTrials & validTrials') = Unq1_concat_hard;
% Unq2_new{n}{pc, 1}(dataSplit.easyTrials & validTrials') = Unq2_concat_easy;
% Unq2_new{n}{pc, 1}(dataSplit.hardTrials & validTrials') = Unq2_concat_hard;
% Syn_new{n}{pc, 1}(dataSplit.easyTrials & validTrials') = Syn_concat_easy;
% Syn_new{n}{pc, 1}(dataSplit.hardTrials & validTrials') = Syn_concat_hard;
% end
%     end
% end
% 
% %         %%
% 
        validPhrases = behavior.phrases(validTrials);
        if(length(corrupted_pca_trials) == (sum(dataSplit.easyTrials(idxs_valid(corrupted_pca_trials))) + sum(dataSplit.hardTrials(idxs_valid(corrupted_pca_trials)))))
            disp("OK. All Corrupted trials were removed from difficulty");
        else
            disp("KO. Not all the corrupted PC were removed from difficulty");
        end
        dataSplit.easyTrials(idxs_valid(corrupted_pca_trials)) = 0;
        dataSplit.hardTrials(idxs_valid(corrupted_pca_trials)) = 0;
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
        if(sum(validTrials) == length(data_epochs.trial))
            disp("OK. The epochs are " + length(data_epochs.trial) + " and are still coincident with the number of validTrials");
        else
            disp("KO. There is a problem with the number of epochs (" + length(data_epochs.trial) + ")");
        end

        if(n == 1)
            lengths_easyTrials_subjects = sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.easyTrials == 1)));
            lengths_hardTrials_subjects = sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.hardTrials == 1)));
        else
            lengths_easyTrials_subjects = [lengths_easyTrials_subjects, sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.easyTrials == 1)))];
            lengths_hardTrials_subjects = [lengths_hardTrials_subjects, sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.hardTrials == 1)))];
        end
        for n_performance=1:NUM_PERFORMANCE
            if(performances(n_performance))
                if(n == 1)
                    lengths_goodPerformanceTrials_subjects{n_performance} = sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.goodPerformanceTrials{n_performance} == 1)));
                    lengths_badPerformanceTrials_subjects{n_performance} = sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.badPerformanceTrials{n_performance} == 1)));
                else
                    lengths_goodPerformanceTrials_subjects{n_performance} = [lengths_goodPerformanceTrials_subjects{n_performance},  sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.goodPerformanceTrials{n_performance} == 1)))];
                    lengths_badPerformanceTrials_subjects{n_performance} = [lengths_badPerformanceTrials_subjects{n_performance},  sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.badPerformanceTrials{n_performance} == 1)))];
                end
            end
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
        start_time=0.5;

        start_sample= Fs*start_time;

        for i=1:length(all_pcas_ordered)
            pcas_timeShifted{i}=all_pcas_ordered{i}(:,start_sample:end)';
            envelopes_timeShifted{i}=all_envelopes_ordered{i}(start_sample:end);
            lengths_timeShifted=[lengths_timeShifted length(envelopes_timeShifted{i})];
        end

        T=find(data_epochs.time{1}==start_time+delay);

        cfg=[];
        cfg.endsample=[];
        cfg.begsample=[];

        for i=1:length(data_epochs.trial)
            cfg.begsample=[cfg.begsample T];
            cfg.endsample=[cfg.endsample T+lengths_timeShifted(i)-1];
        end

        % Reallignment of the EEG signals after the shift of 1.5s and the delay
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
            data_epochs_input.time{1,i}=data_epochs_input.time{1,i}-delay;
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
        if(computePID_difficulty)
            easyTrials = find(dataSplit.easyTrials(validTrials == 1));
            hardTrials = find(dataSplit.hardTrials(validTrials == 1));
            disp("Easy Trials: " + length(easyTrials) + ", hard trials: " + length(hardTrials));
            for k=PCs
                % (length(validTrials == 1) x 1) cell array. Cell i contains the
                % (samples x 1) double matrix representing the k-th PC.
                pcas_timeShifted_k={};

                if(matchLength)
                    min_length_epoch = length(pcas_timeShifted{1}(:, k));
                    for epoch=1:length(pcas_timeShifted)
                        pcas_timeShifted_k=[pcas_timeShifted_k; pcas_timeShifted{epoch}(:,k)];
                        if(length(pcas_timeShifted{epoch}(:,k)) < min_length_epoch)
                            min_length_epoch = length(pcas_timeShifted{epoch}(:,k));
                        end
                    end
                    min_length_epoch
                    % verificare che easyTrials e hardTrials contengano gli
                    % indici e non siano i valori logici delle posizioni
                    for epoch=1:length(easyTrials)
                        disp("Computing PID sub: " + n + " easy trial: " + epoch + "/" + length(easyTrials));
                        eeg_easy = data_eeg{easyTrials(epoch)}(:, 1:min_length_epoch);
                        env_easy = envelopes_timeShifted{easyTrials(epoch)}(1:min_length_epoch);
                        pc_easy = pcas_timeShifted_k{easyTrials(epoch)}(1:min_length_epoch);
                        [red_all_ch_easy{epoch},unq1_all_ch_easy{epoch},unq2_all_ch_easy{epoch},syn_all_ch_easy{epoch}] = PID_copnorm(eeg_easy,env_easy,pc_easy);
                    end
                    for epoch=1:length(hardTrials)
                        disp("Computing PID sub: " + n + " hard trial: " + epoch + "/" + length(hardTrials));
                        eeg_hard = data_eeg{hardTrials(epoch)}(:, 1:min_length_epoch);
                        env_hard = envelopes_timeShifted{hardTrials(epoch)}(1:min_length_epoch);
                        pc_hard = pcas_timeShifted_k{hardTrials(epoch)}(1:min_length_epoch);
                        [red_all_ch_hard{epoch},unq1_all_ch_hard{epoch},unq2_all_ch_hard{epoch},syn_all_ch_hard{epoch}] = PID_copnorm(eeg_hard,env_hard,pc_hard);
                    end
                end
                % {1} = EASY, {2} = HARD
                for epoch=1:length(easyTrials)
                    lagged_redundancy_difficulty{1}{epoch}{k, 1}= red_all_ch_easy{epoch};
                    lagged_unq1_difficulty{1}{epoch}{k, 1} = unq1_all_ch_easy{epoch};
                    lagged_unq2_difficulty{1}{epoch}{k, 1} = unq2_all_ch_easy{epoch};
                    lagged_synergy_difficulty{1}{epoch}{k, 1} = syn_all_ch_easy{epoch};
                end
                for epoch=1:length(hardTrials)
                    lagged_redundancy_difficulty{2}{epoch}{k, 1} = red_all_ch_hard{epoch};
                    lagged_unq1_difficulty{2}{epoch}{k, 1} = unq1_all_ch_hard{epoch};
                    lagged_unq2_difficulty{2}{epoch}{k, 1} = unq2_all_ch_hard{epoch};
                    lagged_synergy_difficulty{2}{epoch}{k, 1} = syn_all_ch_hard{epoch};
                end
            end
        end
        clear red_all_ch_easy unq1_all_ch_easy unq2_all_ch_easy syn_all_ch_easy
        clear red_all_ch_hard unq1_all_ch_hard unq2_all_ch_hard syn_all_ch_hard
        if(computePID_performance)
            for n_performance=1:NUM_PERFORMANCE
                if(performances(n_performance))
                    goodPerformanceTrials = find(dataSplit.goodPerformanceTrials{n_performance}(validTrials == 1));
                    badPerformanceTrials = find(dataSplit.badPerformanceTrials{n_performance}(validTrials == 1));
                    %                     disp("Good performance trials: " + length(goodPerformanceTrials{n_performance}) + ", bad performance trials: " + length(badPerformanceTrials{n_performance}));
                    for k=PCs
                        % (length(validTrials == 1) x 1) cell array. Cell i contains the
                        % (samples x 1) double matrix representing the k-th PC.
                        pcas_timeShifted_k={};

                        for epoch=1:length(pcas_timeShifted)
                            pcas_timeShifted_k=[pcas_timeShifted_k; pcas_timeShifted{epoch}(:,k)];
                        end

                        for epoch=1:length(goodPerformanceTrials)
                            [red_all_ch_good{epoch},unq1_all_ch_good{epoch},unq2_all_ch_good{epoch},syn_all_ch_good{epoch}] = PID_copnorm(data_eeg(goodPerformanceTrials(epoch)),envelopes_timeShifted(goodPerformanceTrials(epoch)),pcas_timeShifted_k(goodPerformanceTrials(epoch)));
                        end

                        for epoch=1:length(badPerformanceTrials)
                            [red_all_ch_bad{epoch},unq1_all_ch_bad{epoch},unq2_all_ch_bad{epoch},syn_all_ch_bad{epoch}] = PID_copnorm(data_eeg(badPerformanceTrials(epoch)),envelopes_timeShifted(badPerformanceTrials(epoch)),pcas_timeShifted_k(badPerformanceTrials(epoch)));
                        end

                        for epopch=1:length(goodPerformanceTrials)
                            lagged_redundancy_performance{n_performance}{1}{epoch}{k, 1} = red_all_ch_good{epoch};
                            lagged_unq1_performance{n_performance}{epoch}{1}{k, 1} = unq1_all_ch_good{epoch};
                            lagged_unq2_performance{n_performance}{epoch}{1}{k, 1} = unq2_all_ch_good{epoch};
                            lagged_synergy_performance{n_performance}{epoch}{1}{k, 1} = syn_all_ch_good{epoch};
                        end

                        for epoch=1:length(badPerformanceTrials)
                            lagged_redundancy_performance{n_performance}{2}{epoch}{k, 1} = red_all_ch_bad{epoch};
                            lagged_unq1_performance{n_performance}{epoch}{2}{k, 1} = unq1_all_ch_bad{epoch};
                            lagged_unq2_performance{n_performance}{epoch}{2}{k, 1} = unq2_all_ch_bad{epoch};
                            lagged_synergy_performance{n_performance}{epoch}{2}{k, 1} = syn_all_ch_bad{epoch};
                        end
                    end
                    lagged_redundancy_subjects_performance{n_performance}{1}{n} = lagged_redundancy_performance{n_performance}{1};
                    lagged_unq1_subjects_performance{n_performance}{1}{n} = lagged_unq1_performance{n_performance}{1};
                    lagged_unq2_subjects_performance{n_performance}{1}{n} = lagged_unq2_performance{n_performance}{1};
                    lagged_synergy_subjects_performance{n_performance}{1}{n} = lagged_synergy_performance{n_performance}{1};
                end
            end
        end
        clear red_all_ch_good unq1_all_ch_good unq2_all_ch_good syn_all_ch_good
        clear red_all_ch_bad unq1_all_ch_bad unq2_all_ch_bad syn_all_ch_bad
        clear lagged_redundancy_performance lagged_unq1_performance lagged_unq2_performance lagged_synergy_performance
        if(computePID_difficulty)
            lagged_redundancy_subjects_difficulty{1}{n}=lagged_redundancy_difficulty{1};
            lagged_unq1_subjects_difficulty{1}{n}=lagged_unq1_difficulty{1};
            lagged_unq2_subjects_difficulty{1}{n}=lagged_unq2_difficulty{1};
            lagged_synergy_subjects_difficulty{1}{n}=lagged_synergy_difficulty{1};

            lagged_redundancy_subjects_difficulty{2}{n}=lagged_redundancy_difficulty{2};
            lagged_unq1_subjects_difficulty{2}{n}=lagged_unq1_difficulty{2};
            lagged_unq2_subjects_difficulty{2}{n}=lagged_unq2_difficulty{2};
            lagged_synergy_subjects_difficulty{2}{n}=lagged_synergy_difficulty{2};
        end

        clear lagged_redundancy_difficulty lagged_unq1_difficulty lagged_unq2_difficulty lagged_synergy_difficulty
    end
    cd 'results/PID/'

    %% SAVING RESULTS
    fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_noSplit_";
    fname_difficulty = "difficultySplit_";
    fname_performance = "performanceSplit_";
    fname_bandp = "singleTrial_matchedLength_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs";

    if(computePID_difficulty)
        fname = "REDUNDANCY_" + fname_offset + fname_difficulty + fname_bandp;
        Rdn = lagged_redundancy_subjects_difficulty;
        save(fname + ".mat", "Rdn");
        fname = "UNIQUE1_" + fname_offset + fname_difficulty + fname_bandp;
        Unq1 = lagged_unq1_subjects_difficulty;
        save(fname + ".mat", "Unq1");
        fname = "UNIQUE2_" + fname_offset + fname_difficulty + fname_bandp;
        Unq2 = lagged_unq2_subjects_difficulty;
        save(fname + ".mat", "Unq2");
        fname = "SYNERGY_" + fname_offset + fname_difficulty + fname_bandp;
        Syn = lagged_synergy_subjects_difficulty;
        save(fname + ".mat", "Syn");
    end

    if(computePID_performance)
        fname = "REDUNDANCY_" + fname_offset + fname_performance + fname_bandp;
        save(fname + ".mat", "lagged_redundancy_subjects_performance");
        fname = "UNIQUE1_" + fname_offset + fname_performance + fname_bandp;
        save(fname + ".mat", "lagged_unq1_subjects_performance");
        fname = "UNIQUE2_" + fname_offset + fname_performance + fname_bandp;
        save(fname + ".mat", "lagged_unq2_subjects_performance");
        fname = "SYNERGY_" + fname_offset + fname_performance + fname_bandp;
        save(fname + ".mat", "lagged_synergy_subjects_performance");
    end

    chs=data_epochs.label;
    save('channels_list','chs');
    cd 'D:/Entrainment__/final_functions_behavior/'
end

