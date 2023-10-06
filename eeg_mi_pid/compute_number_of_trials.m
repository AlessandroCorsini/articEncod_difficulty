clc;
clear;
close all;

%% ADDING PATHS
addpath '../fieldtrip-20210311'
addpath '../gcmi-master/matlab'
addpath '../analisi'

%% PARAMETERS
rng(1000)
rng(1000)
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
MATCH_LENGTH = false;
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
SINGLE_TRIAL_LENGTH = 'max';
% Number of iteration of PID. In combination with match_length = true at
% every computation a different segment is cutted for matching
MI_ITERATIONS = 1;
flows = [5];
fhighs = [7];
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

num_trials = {};
parfor n = 1:22
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
                t0 = data_epochs.time{epoch}(1);
                data_epochs.time{epoch}=round(data_epochs.time{epoch}(1,start_cutTime:end-end_cutTime)+ t0 - data_epochs.time{epoch}(1,start_cutTime), 5);
                data_epochs.trial{epoch}=data_epochs.trial{epoch}(:,start_cutTime:end-end_cutTime);
            elseif(goodTrials_original(trial) == 1)
                % Remove the invalid trials (phrase 34 and the NaN answers)
                % that were instead considered as goodTrials in the
                % original version
                epoch = sum(find(goodTrials_original) <= trial);
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
    for k=1:length(all_pcas_ordered)
        if length(all_pcas_ordered{k}) < 2000
            % The index in corrupted_pca_trials is referred to the index of
            % the trial in the validTrials only
            corrupted_pca_trials = [corrupted_pca_trials k];
        end
    end
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
        dataSplit.easyTrials(idxs_valid(corrupted_pca_trials)) = 0;
        dataSplit.hardTrials(idxs_valid(corrupted_pca_trials)) = 0;
    end

    class1 = find(dataSplit.easyTrials(validTrials == 1));
    class2 = find(dataSplit.hardTrials(validTrials == 1));

    num_trials{n} = length(class1) + length(class2);
end

