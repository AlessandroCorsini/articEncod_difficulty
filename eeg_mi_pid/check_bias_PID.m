clear
clc
cd 'D:\Entrainment__\final_functions_behavior\'

%% ADDING PATHS AND SETTING PARAMS
rng(1000)

addpath '..\partial-info-decomp-master'
addpath '..\fieldtrip-20210311'
addpath '..\gcmi-master\matlab'
addpath '..\analisi'

% DEFINITION OF PID CALCULATION
REALIGN_INPUTS=true;
REALIGN_EEG=true;
ITERATIONS = [1,2,4,8,16,32,64,128,256] ;

flows=[0.5];
fhighs=[10];
% The lag is between the Stimuli (speech
% envelope, PCs) and the Response (EEG channel) which is supposed to be
% delayed
Fs=400;
delay=0.2;%s
lag=80; %delay in sample -> i.e. 200ms at 400hz
start_time = 0.5;
PCs = [1]; % PCs to use
pca_dir="../pca_new_nocut";
env_dir="../envelope";

load('cut_times')
cut_times=round(cut_times*Fs*0.001); %from ms to samples

%% LOADING DATA NAMES EEG

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
N_SUBJ = length(subs);

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
    if REALIGN_INPUTS
        % read the time from where to cut files to eliminate silences present
        % at the beginning of sentences
        start_cutTime=cut_times(i,1);
        end_cutTime=cut_times(i,2);
        envelopes{1,i}=envelopes{1,i}(start_cutTime:end-end_cutTime);
        pcas{i}=pcas{i}(:,start_cutTime:end-end_cutTime);
    end
end

for freq=1:length(flows)
    flow = flows(freq);
    fhigh = fhighs(freq);
    disp("**flow: " + flow + ", fhigh: " + fhigh + "**");
    fname_redundancy = "REDUNDANCY_";
    fname_unq1 = "UNIQUE1_";
    fname_unq2 = "UNIQUE2_";
    fname_synergy = "SYNERGY_";
    fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_";
    fname_checkbias = "alltrials_checkbias_";
    fname_bandp = "trials_concat_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs";
    Rdn_subs ={};
    Unq1_subs ={};
    Unq2_subs ={};
    Syn_subs ={};
    for n=1:length(subs)
        disp("*Processing subject " + subs(n) + "*");
        Rdn={};
        Unq1={};
        Unq2={};
        Syn={};
        cfg=[];
        file=subs_files{n};
        cfg.dataset=file;
        data_epochs=ft_preprocessing(cfg)
        behavior = load("..\behavioral_\results\subs\" + subs{n} + "_behavior").behavior;
        goodTrials_original = zeros(1, length(behavior.goodTrials));
        goodTrials_original(load("..\analisi\" + subs{n} + "_goodTrials.mat").good_trials + 1) = 1;
        validTrials_original = goodTrials_original == 1 & (behavior.answers == 0 | behavior.answers == 1);
        validTrials = behavior.goodTrials == 1 & (behavior.answers == 0 | behavior.answers == 1);
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
        %REMOVING TRIALS WHERE THE ANSWER TO THE ATTENTIVE QUESTION WAS
        %INCORRECT OR NOT VALID
        data_epochs.time(epochs_invalid) = [];
        data_epochs.trial(epochs_invalid) = [];
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
        for k=PCs
            % (length(validTrials == 1) x 1) cell array. Cell i contains the
            % (samples x 1) double matrix representing the k-th PC.
            pcas_timeShifted_k={};
            for epoch=1:length(pcas_timeShifted)
                pcas_timeShifted_k=[pcas_timeShifted_k; pcas_timeShifted{epoch}(:,k)];
            end
            data_eeg=data_epochs.trial;

            for trial=1:length(pcas_timeShifted_k)
                if trial==1
                    eeg=data_eeg{trial};
                    env=envelopes_timeShifted{trial};
                    pc=pcas_timeShifted_k{trial};
                else
                    eeg=[eeg data_eeg{trial}];
                    env=[env; envelopes_timeShifted{trial}];
                    pc=[pc; pcas_timeShifted_k{trial}];
                end
            end
            disp("START PID SUBJECT: " + subs_n(n));
            index = 1;
            for iter=ITERATIONS
                disp("Dividing the total length into " + iter + " pieces");
                num_samples = floor(length(eeg)/iter);
                for split=1:iter
                    disp("Computing the PID for piece: " + iter + "_" + split);
                    start_sample = (split - 1)*num_samples + 1;
                    end_sample = split*num_samples;
                    [red_all_ch,unq1_all_ch,unq2_all_ch,syn_all_ch] = PID_copnorm(eeg(:,start_sample:end_sample),env(start_sample:end_sample),pc(start_sample:end_sample));
                    if(split ==1)
                        Rdn{k, index} = [red_all_ch];
                        Unq1{k, index} = [unq1_all_ch];
                        Unq2{k, index} = [unq2_all_ch];
                        Syn{k, index} = [syn_all_ch];
                    else
                        Rdn{k, index} = [Rdn{k, index}, red_all_ch];
                        Unq1{k, index} = [Unq1{k, index}, unq1_all_ch];
                        Unq2{k, index} = [Unq2{k, index}, unq2_all_ch];
                        Syn{k, index} = [Syn{k, index}, syn_all_ch];
                    end
                end
                index = index + 1;
            end
        end
        Rdn_subs{n}=Rdn;
        Unq1_subs{n}=Unq1;
        Unq2_subs{n}=Unq2;
        Syn_subs{n}=Syn;
    end
    %% SAVING RESULTS
    cd 'results/PID/'
    Rdn = Rdn_subs;
    fname = fname_redundancy + fname_offset + fname_checkbias + fname_bandp;
    save(fname + ".mat", "Rdn");
    Unq1 = Unq1_subs;
    fname = fname_unq1 + fname_offset + fname_checkbias + fname_bandp;
    save(fname + ".mat", "Unq1");
    Unq2 = Unq2_subs;
    fname = fname_unq2 + fname_offset + fname_checkbias + fname_bandp;
    save(fname + ".mat", "Unq2");
    Syn = Syn_subs;
    fname = fname_synergy + fname_offset + fname_checkbias + fname_bandp;
    save(fname + ".mat", "Syn");
    cd 'D:\Entrainment__\final_functions_behavior\'
end