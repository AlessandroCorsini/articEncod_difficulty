clear
clc
close all

%% ADDING PATHS
addpath '../fieldtrip-20210311'
addpath '../gcmi-master/matlab'
addpath '../analisi'

%% PARAMETERS
rng(4)
% realign_inputs true means the cut of the silent parts at the end and at 
% the beginning of the envelopes and PCs is performed.
REALIGN_INPUTS=true;
% relign_eeg true means the cut of the silent part is also applied to the 
% eeg signals 
REALIGN_EEG=true;
% Compute the MI of the split based on the diffculty of the phrases.
% Otherwise no split, but one single PID on all concatenated trials.
COMPUTE_MI_DIFFICULTY = true;
% Compute PID on concatenated trials splitted by difficulty or performance
% matching the length of the two concatenations
MATCH_LENGTH = true;
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
SINGLE_TRIAL_LENGTH = 'max'; % CHANGE ALSO IN THE PARFOPR LOOP
% Number of iteration of PID. In combination with match_length = true at
% every computation a different segment is cutted for matching
MI_ITERATIONS = 200;
flows = [0.5];
fhighs = [4];
% The lag is between the Stimuli (speech
% envelope, PCs) and the Response (EEG channel) which is supposed to be
% delayed
Fs=400;
delay=0.2;%s
lag=80; %delay in sample -> i.e. 200ms at 400hz
start_time = 0.5;
PCs = [1]; % number of PCs to use in the analysis 
pca_dir="../pca_new_nocut"; 
env_dir="../envelope";
% The cut times are in ms.
load('cut_times')
cut_times=round(cut_times*Fs*0.001); %from ms to samples

%% LOADING DATA NAMES EEG
% The usual loading of the names of the files with completed analysis and
% of the names of the subjects
trigger_name = "startSound";
d=dir('../analisi');
subs_files ={};
subs={};
subs_n={};
for i=3:length(d)
    if endsWith(d(i).name,"_filtered_epoched" + trigger_name + "_ica_interp-epo.fif")
        subs_files=[subs_files; d(i).name];
        tmp = split(d(i).name,'_f');
        tmp2 = split(tmp(1),'_');
        subs = [subs; tmp(1)];
        subs_n = [subs_n; tmp2(2)];
    end
end
subs= subs(SELECTED_SUB);
subs_files=subs_files(SELECTED_SUB);
subs_n=subs_n(SELECTED_SUB);
N_SUBJ = length(subs);

%% LOADING PCA NAMES AND DATA
d_pca_em=dir(pca_dir);
% (NUM_PHRASES x 1) cell array. The cell i contains the name of the
% file containing the values of all the 21 PCs of the phrase i.
pca_names={};
for i=4:length(d_pca_em)
     pca_names=[pca_names; d_pca_em(i).name];
end
pca_names=natsort(pca_names);
% (1 x length(pca_names)) cell array. The cell (1, i) contains a
% (21 x samples) double matrix with the values of all the 21 PCs of the
% phrase i.
pcas={};
for i=1:length(pca_names)
        pca_tmp=load(pca_dir+"/"+pca_names{i});
        pcas=[pcas pca_tmp.pca_tmp];
end

%% LOADING ENVELOPE NAMES AND DATA
d_env=dir(env_dir);
% (NUM_PHRASES x 1) cell array. The cell i contains the names of all
% the files of the envelopes of the phrases
envelope_names={};
for i=3:length(d_env)
     envelope_names=[envelope_names; d_env(i).name];
end
envelope_names=natsort(envelope_names);
% (1 x NUM_PHRASES) cell array. The cell (1, i) contains a
% (1 x samples) double matrix with the values of the envelope at each
% samples for phrase i.
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
% lengths in samples of the easy trials concatenated for each subject 
lengths_easyTrials_subjects = [];
% lengths in samples of the hard trials concatenated for each subject 
lengths_hardTrials_subjects = [];
% lengths in samples of the good trials concatenated for each subject 
lengths_goodPerformanceTrials_subjects = [];
% lengths in samples of the bad trials concatenated for each subject 
lengths_badPerformanceTrials_subjects = [];
disp("***START OF THE CORE PART OF THE ANALYSIS***");
for freq=1:length(flows)
    flow = flows(freq);
    fhigh = fhighs(freq);
    disp("**flow: " + flow + ", fhigh: " + fhigh + "**");
    fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_";
    fname_difficulty = "difficultySplit_";
    fname_bandp = "trials_concat_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";
    if(CONCATENATE_SINGLE_TRIAL_MATCHED_LENGTH)
        fname_bandp = "trials_concat_singleTrial_matchedLength_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";
    end
    if(MATCH_LENGTH)
        fname_difficulty = fname_difficulty + "matchedLength_801_1000";
    end
    % MI Difficulty or Performance
    MI_env_subs = {};
    MI_pc_subs = {};
    epochs_invalid = [];
    parfor n=SELECTED_SUB

        % repeat for 'parfor' loop
        SINGLE_TRIAL_LENGTH = 'max';
        
        disp("*Processing subject " + subs(n) + "*");
        % MI Difficulty or Performance
        MI_env = {};
        MI_pc = {};
        
        %% METAFILES LOADING
        file = subs_files{n};
        cfg = [];
        cfg.dataset=file;
        data_epochs=ft_preprocessing(cfg);
        % Loading the behavioral data
        behavior = load("../behavioral_/results/subs/" + subs{n} + "_behavior").behavior;
        dataSplit = load("../behavioral_/results/subs/" + subs{n} + "_dataSplit").dataSplit;
        goodTrials_original = zeros(1, length(behavior.goodTrials));
        goodTrials_original(load("../analisi/" + subs{n} + "_goodTrials.mat").good_trials + 1) = 1;
        disp("The original good trials were: " + sum(goodTrials_original));
        disp("Now the good trials are: " + sum(behavior.goodTrials));      
        validTrials_original = goodTrials_original == 1 & (behavior.answers == 0 | behavior.answers == 1);
        validTrials = behavior.goodTrials == 1 & (behavior.answers == 0 | behavior.answers == 1);
        
        %% ADJUSTING VECTOR OF TIME AND TRIALS USING THE SAME TIMES USED FOR ENVELOPES AND PCAs    
        if REALIGN_EEG
            epochs_invalid = [];
            for trial=1:length(validTrials) % all trials
                % Using validTrials and not validTrials_original I exclude
                % phrase 34.
                if(validTrials(trial))
                    start_cutTime = cut_times(behavior.phrases(1, trial),1);
                    end_cutTime = cut_times(behavior.phrases(1, trial),2);
                    % No problem even after the exclusion of the
                    % phrase 34 because the original goodTrials are stored in
                    % the "analisi" folder.
                    % data_epochs.time{t} is a (1 x samples_of_the_trial_t)
                    % double matrix with the time axis of the trial t.
                    % t0 contains the first value of the time axis for the
                    % trial t
                    % data_epochs is made of just the goodTrials epochs.
                    % This line of code converts the number of the trial
                    % into the corresponding number of epoch
                    epoch = sum(find(goodTrials_original) <= trial);
                    disp("Trial " + trial + " is valid and corresponds to epoch " + epoch);
                    t0 = data_epochs.time{epoch}(1);  
                    data_epochs.time{epoch}=round(data_epochs.time{epoch}(1,start_cutTime:end-end_cutTime)+ t0 - data_epochs.time{epoch}(1,start_cutTime), 5);
                    data_epochs.trial{epoch}=data_epochs.trial{epoch}(:,start_cutTime:end-end_cutTime);
                elseif(goodTrials_original(trial) == 1)
                    % Remove the invalid trials (phrase 34 and the NaN answers)
                    % that were instead considered as goodTrials in the
                    % original version
                    epoch = sum(find(goodTrials_original) <= trial);
                    disp("Trial " + trial + "is not valid but originally good and corresponds to epoch " + epoch);
                    epochs_invalid = [epochs_invalid, epoch];
                end
            end
        end
        data_epochs.time(epochs_invalid) = [];
        data_epochs.trial(epochs_invalid) = [];
        % NOW THE INDEXING IN data_epochs COINCIDES WITH THE INDEXING IN THE
        % SELECTION OF THE validTrials ONLY.
    
        %% ORDERING PCAs AND ENVELOPES TRIALS (AS PROVIDED IN THE EXPERIMENT) TO MATCH THEM WITH EEG TRIALS      
        % (1 x length(validTrials == 1)) cell array. Cell i contains the (21 x
        % samples) double matrix referred to the PCs of the phrase number
        % validPhrases(i)
        all_pcas_ordered={};
        % (1 x length(validTrials == 1)) cell array. Cell i contains the (1 x
        % samples) double matrix referred to the values of the envelope of
        % phrase number validPhrases(i)
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
        validPhrases = behavior.phrases(validTrials);
        if(COMPUTE_MI_DIFFICULTY)
            if(length(corrupted_pca_trials) == (sum(dataSplit.easyTrials(idxs_valid(corrupted_pca_trials))) + sum(dataSplit.hardTrials(idxs_valid(corrupted_pca_trials)))))
                disp("OK. All Corrupted trials were removed from difficulty");
            else
                disp("KO. Not all the corrupted PC were removed from difficulty");
            end
            dataSplit.easyTrials(idxs_valid(corrupted_pca_trials)) = 0;
            dataSplit.hardTrials(idxs_valid(corrupted_pca_trials)) = 0;
        end
        if(sum(validTrials) == length(data_epochs.trial))
            disp("OK. The epochs are " + length(data_epochs.trial) + " and are still coincident with the number of validTrials");
        else
            disp("KO. There is a problem with the number of epochs (" + length(data_epochs.trial) + ")");
        end
        % %% COMPUTING META DATA
        % if(COMPUTE_MI_DIFFICULTY)
        %     lengths_easyTrials_subjects = [lengths_easyTrials_subjects, sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.easyTrials == 1)))];
        %     lengths_hardTrials_subjects = [lengths_hardTrials_subjects, sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.hardTrials == 1)))];
        % end
        
        %% REMOVING start_time seconds TO EXCLUDE AUDITORY POTENTIALS AT THE TRIAL'S ONSET
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
                       
        %% FILTERING INPUTS (SPEECH ENVELOPE AND PCS) FROM flow TO fhigh         
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
               
        %% COMPUTING MI
        data_eeg = data_epochs.trial;
        % MI DIFFICULTY
        if(COMPUTE_MI_DIFFICULTY)            
            class1 = find(dataSplit.easyTrials(validTrials == 1));
            class2 = find(dataSplit.hardTrials(validTrials == 1));
        end
        for mi_iter=1:MI_ITERATIONS
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
                if(mi_iter==1)
                    if(CONCATENATE_SINGLE_TRIAL_MATCHED_LENGTH && strcmp(SINGLE_TRIAL_LENGTH, 'max'))
                        SINGLE_TRIAL_LENGTH = min_length_epoch;
                    end
                    for trial=1:length(class1)
                        if(not(CONCATENATE_SINGLE_TRIAL_MATCHED_LENGTH))
                            SINGLE_TRIAL_LENGTH = length(data_eeg{class1(trial)});
                        end
                        if trial==1
                            eeg_class1=data_eeg{class1(trial)}(:, 1:SINGLE_TRIAL_LENGTH);
                            env_class1=envelopes_timeShifted{class1(trial)}(1:SINGLE_TRIAL_LENGTH);
                            pc_class1=pcas_timeShifted_k{class1(trial)}(1:SINGLE_TRIAL_LENGTH);
                        else
                            eeg_class1=[eeg_class1 data_eeg{class1(trial)}(:, 1:SINGLE_TRIAL_LENGTH)];
                            env_class1=[env_class1; envelopes_timeShifted{class1(trial)}(1:SINGLE_TRIAL_LENGTH)];
                            pc_class1=[pc_class1; pcas_timeShifted_k{class1(trial)}(1:SINGLE_TRIAL_LENGTH)];
                        end
                    end
                    for trial=1:length(class2)
                        if(not(CONCATENATE_SINGLE_TRIAL_MATCHED_LENGTH))
                            SINGLE_TRIAL_LENGTH = length(data_eeg{class2(trial)});
                        end
                        if trial==1
                            eeg_class2=data_eeg{class2(trial)}(:, 1:SINGLE_TRIAL_LENGTH);
                            env_class2=envelopes_timeShifted{class2(trial)}(1:SINGLE_TRIAL_LENGTH);
                            pc_class2=pcas_timeShifted_k{class2(trial)}(1:SINGLE_TRIAL_LENGTH);
                        else
                            eeg_class2=[eeg_class2 data_eeg{class2(trial)}(:, 1:SINGLE_TRIAL_LENGTH)];
                            env_class2=[env_class2; envelopes_timeShifted{class2(trial)}(1:SINGLE_TRIAL_LENGTH)];
                            pc_class2=[pc_class2; pcas_timeShifted_k{class2(trial)}(1:SINGLE_TRIAL_LENGTH)];
                        end
                    end
                end                
                if(MATCH_LENGTH)
                    %% ONE SINGLE MATCHED_LENGTH COMPUTATION
                    eeg_class1_matched = eeg_class1;
                    eeg_class2_matched=eeg_class2;
                    env_class1_matched=env_class1;
                    env_class2_matched=env_class2;
                    pc_class1_matched=pc_class1;
                    pc_class2_matched=pc_class2;
                    rng(mi_iter)
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
                    disp("ITERATION: " + mi_iter);
                    disp("SAMPLES DELETED: (" + start_sample + ", " + (start_sample+diff) + ")");
                    disp("eeg_class1: " + length(eeg_class1_matched) + " samples");
                    disp("env_class1: " + length(env_class1_matched) + " samples");
                    disp("pc_class1: " + length(pc_class1_matched) + " samples");
                    disp("eeg_class2: " + length(eeg_class2_matched) + " samples");
                    disp("env_class2: " + length(env_class2_matched) + " samples");
                    disp("pc_class2: " + length(pc_class2_matched) + " samples");
                    if(k == PCs(1))
                        [~, MI_pc_class1, ~] = MI_copnorm(eeg_class1_matched,[],pc_class1_matched, true, false, false);
                        [MI_env_class1,~, ~] = MI_copnorm(eeg_class1_matched,env_class1_matched,[], false, true, false);
                        [~, MI_pc_class2, ~] = MI_copnorm(eeg_class2_matched,[],pc_class2_matched, true, false, false);
                        [MI_env_class2,~, ~] = MI_copnorm(eeg_class2_matched,env_class2_matched,[], false, true, false);
                        MI_env = [MI_env_class1, MI_env_class2];
                        MI_pc{k, 1} = [MI_pc_class1, MI_pc_class2];
                    else
                        [~, MI_pc_class1, ~] = MI_copnorm(eeg_class1_matched,[],pc_class1_matched, true, false, false);
                        [~, MI_pc_class2, ~] = MI_copnorm(eeg_class2_matched,[],pc_class2_matched, true, false, false);
                        MI_pc{k, 1} = [MI_pc_class1, MI_pc_class2];
                    end
                else
                    %% ONE SINGLE NON MATCHED_LENGTH COMPUTATION
                    if(k == PCs(1))
                        [~, MI_pc_class1, ~] = MI_copnorm(eeg_class1,[],pc_class1, true, false, false);
                        [MI_env_class1,~, ~] = MI_copnorm(eeg_class1,env_class1,[], false, true, false);
                        [~, MI_pc_class2, ~] = MI_copnorm(eeg_class2,[],pc_class2, true, false, false);
                        [MI_env_class2,~, ~] = MI_copnorm(eeg_class2,env_class2,[], false, true, false);
                        MI_env = [MI_env_class1, MI_env_class2];
                        MI_pc{k, 1} = [MI_pc_class1, MI_pc_class2];
                    else
                        [~, MI_pc_class1, ~] = MI_copnorm(eeg_class1,[],pc_class1, true, false, false);
                        [~, MI_pc_class2, ~] = MI_copnorm(eeg_class2,[],pc_class2, true, false, false);
                        MI_pc{k, 1} = [MI_pc_class1, MI_pc_class2];
                    end
                end               
            end
            if(MATCH_LENGTH)
                MI_env_subs{n}{mi_iter} = MI_env;
                MI_pc_subs{n}{mi_iter} = MI_pc;
            else
                MI_env_subs{n} = MI_env;
                MI_pc_subs{n} = MI_pc;          
            end
        end
    end   
%     if(computeMI_difficulty)
%         figure
%         plot(lengths_easyTrials_subjects/1000, '-o')
%         hold on
%         plot(lengths_hardTrials_subjects/1000, '-o')
%         legend("Easy", "Hard")
%         hold off
%         savefig("dataSplitDifficulty.fig");
%     end

%     if(computeMI_performance)
%         figure
%         plot(lengths_goodPerformanceTrials_subjects/1000, '-o')
%         hold on
%         plot(lengths_badPerformanceTrials_subjects/1000, '-o')
%         legend("Good", "Bad")
%         hold off
%         savefig("dataSplitPerformance.fig");
%     end   
    %% SAVING RESULTS
    cd 'results/MI/'
    if(COMPUTE_MI_DIFFICULTY)
        fname = "MI_ENVELOPE_" + fname_offset + fname_difficulty + fname_bandp;
        Env = MI_env_subs;
        save(fname,'Env');
        fname = "MI_PCs_" + fname_offset + fname_difficulty + fname_bandp;
        PCs = MI_pc_subs;
        save(fname,'PCs');
    end
    cd '../../'
end

