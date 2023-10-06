clear
close all
clc
% cd 'D:\Entrainment__\final_functions_behavior\'

%% ADDING PATHS
addpath '../fieldtrip-20210311'
addpath './results'
addpath '../analisi'
addpath(pwd)

%% PARAMETERS
cmap = brewermap(1024,'RdBu');
cmap = flip(cmap);
PCs = [1];
flows = [0.5];
fhighs = [4];
N_SUBJ = 22;
COMPUTE_STATISTICS_WITHIN_SUBJECTS = true
start_time = 0.5
COMPUTE_STATISTICS_DIFFICULTY = true
COMPUTE_STATISTICS_PERFORMANCE = false
TRIALS_CONCAT_SINGLE_TRIAL_MATCHED_LENGTH = false;
SAVE_DATA = false;
TIME = true;
TIME_FREQ = false;
% Loading some templates
load 'chs_eeg_ok'
load('freq_template')
load('neighbours_ordered')
if(TIME)
    for freq=1:length(flows)
        %% FILE LOADING
        cd 'results/MI/'
        flow = flows(freq);
        fhigh = fhighs(freq);
        fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_";
        fname_difficulty = "difficultySplit_";
        fname_performance = "performanceSplit_";
        if(TRIALS_CONCAT_SINGLE_TRIAL_MATCHED_LENGTH)
            fname_bandp = "trials_concat_singleTrial_matchedLength_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";
        else
            fname_bandp = "trials_concat_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";
        end
        if(COMPUTE_STATISTICS_DIFFICULTY)
            fname_MI_env_difficulty = "MI_ENVELOPE_" + fname_offset + fname_difficulty + fname_bandp;
            fname_MI_PCs_difficulty = "MI_PCs_" + fname_offset + fname_difficulty + fname_bandp;
            MI_env = load(fname_MI_env_difficulty).Env;
            MI_pc = load(fname_MI_PCs_difficulty).PCs;
        end
        if(COMPUTE_STATISTICS_PERFORMANCE)
            fname_MI_env_performance = "MI_ENVELOPE_" + fname_offset + fname_performance + fname_bandp;
            fname_MI_PCs_performance ="MI_PCs_" + fname_offset + fname_performance + fname_bandp;
            lagged_MI_env_subjects_performance = load(fname_MI_env_performance).lagged_MI_env_subjects_performance;
            lagged_MI_pc_subjects_performance = load(fname_MI_PCs_performance).lagged_MI_pc_subjects_performance;
        end
        
        
        %% PREPARING MI DATA
        % MI difficulty
        if(COMPUTE_STATISTICS_DIFFICULTY)
            MI_env_class1 = zeros(N_SUBJ, 64);
            MI_env_class2 = zeros(N_SUBJ, 64);
            for sub=1:N_SUBJ
                MI_env_class1(sub, :) = MI_env{sub}(:, 1)';
                MI_env_class2(sub, :) = MI_env{sub}(:, 2)';
            end
            for pc=PCs
                MI_pc_class1{pc} = zeros(N_SUBJ, 64);
                MI_pc_class2{pc} = zeros(N_SUBJ, 64);
                for sub=1:N_SUBJ
                    MI_pc_class1{pc}(sub, :) = MI_pc{sub}{pc, 1}(:, 1)';
                    MI_pc_class2{pc}(sub, :) = MI_pc{sub}{pc, 1}(:, 2)';
                end
            end
        end
        % ENV performance
        if(COMPUTE_STATISTICS_PERFORMANCE)
            MI_env_subjects_good = zeros(N_SUBJ, 64);
            lagged_MI_env_subjects_bad = zeros(N_SUBJ, 64);
            for sub=1:N_SUBJ
                lagged_MI_env_subjects_good(sub, :) = lagged_MI_env_subjects_performance{sub}{1, 1}(:, 1)';
                lagged_MI_env_subjects_bad(sub, :) = lagged_MI_env_subjects_performance{sub}{1, 1}(:, 2)';
            end
        end
        % PC performance
        if(COMPUTE_STATISTICS_PERFORMANCE)
            lagged_MI_pc_subjects_good = zeros(N_SUBJ, 64);
            lagged_MI_pc_subjects_bad = zeros(N_SUBJ, 64);
            for sub=1:N_SUBJ
                lagged_MI_pc_subjects_good(sub, :) = lagged_MI_pc_subjects_performance{sub}{pc, 1}(:, 1)';
                lagged_MI_pc_subjects_bad(sub, :) = lagged_MI_pc_subjects_performance{sub}{pc, 1}(:, 2)';
            end
        end
        
        %% REMOVING NON EEG CHANNELS
        bad_chs=[22,28,32,41,46];
        freq_fourier_all.label = cellstr(chs_eeg_ok);
        
        % MI env difficulty
        MI_env_class1(:, bad_chs) = [];
        MI_env_class2(:, bad_chs) = [];
        for pc=PCs
            MI_pc_class1{pc}(:, bad_chs) = [];
            MI_pc_class2{pc}(:, bad_chs) = [];
        end
        
        %% CLUSTER-BASED STATISTICS WITHIN SUBJECTS
        cd '../../'
        if COMPUTE_STATISTICS_WITHIN_SUBJECTS
            
            cfg=[];
            cfg.method           = 'montecarlo';
            cfg.statistic        = 'depsamplesT';
            cfg.correctm         = 'cluster';
            cfg.clusteralpha     = 0.05;
            cfg.clusterstatistic = 'maxsum';
            cfg.minnbchan        = 2;
            cfg.neighbours       = neighbours_ordered;  % same as defined for the between-trials experiment
            cfg.tail             = 0;
            cfg.clustertail      = 0;
            cfg.alpha            = 0.05;
            cfg.correcttail = 'alpha';
            cfg.numrandomization = 5000;     % number of draws from the permutation distribution
            design = zeros(2, N_SUBJ*2);
            design(1,:) = [1:N_SUBJ 1:N_SUBJ];
            design(2,:) = [ones(1,N_SUBJ) ones(1,N_SUBJ)*2];
            cfg.design = design;
            cfg.uvar   = 1;
            cfg.ivar   = 2;
            
            allSubjEnvClass1.label  = chs_eeg_ok;
            allSubjEnvClass1.dimord = 'subj_chan_time';
            allSubjEnvClass1.time   = 1;
            allSubjEnvClass1.avg = MI_env_class1;
            allSubjEnvClass2.label  = chs_eeg_ok;
            allSubjEnvClass2.dimord = 'subj_chan_time';
            allSubjEnvClass2.time   = 1;
            allSubjEnvClass2.avg = MI_env_class2;
            for pc= PCs
                allSubjPCsClass1{pc}.label  = chs_eeg_ok;
                allSubjPCsClass1{pc}.dimord = 'subj_chan_time';
                allSubjPCsClass1{pc}.time   = 1;
                allSubjPCsClass1{pc}.avg = MI_pc_class1{pc};
                allSubjPCsClass2{pc}.label  = chs_eeg_ok;
                allSubjPCslass2{pc}.dimord = 'subj_chan_time';
                allSubjPCsClass2{pc}.time   = 1;
                allSubjPCsClass2{pc}.avg = MI_pc_class2{pc};
            end
            
            %% Statistics
            [stat_env] = ft_timelockstatistics(cfg, allSubjEnvClass1, allSubjEnvClass2);
            h_MI_env_np=stat_env.mask;
            for pc = PCs
                [stat_pc{pc}] = ft_timelockstatistics(cfg, allSubjPCsClass1{pc}, allSubjPCsClass2{pc});
                h_MI_pc_np{pc} = stat_pc{pc}.mask;
            end
            
            %% Topoplot
            cfg_plot=[];
            cfg_plot.colormap=cmap;
            cfg_plot.marker =  'off'
            cfg_plot.highlight='on'
            cfg_plot.highlightchannel=[]
            cfg_plot.highlightsymbol='.';
            cfg_plot.highlightcolor='k';
            cfg_plot.highlightsize=18;
            cfg_plot.zlim = [-0.0003 0.0003];
            freq_fourier_all.freq=[1];
            figure
            freq_fourier_all.powspctrm = mean(MI_env_class1) - mean(MI_env_class2);
            cfg_plot.highlightchannel = freq_fourier_all.label(h_MI_env_np);
            ft_topoplotER(cfg_plot, freq_fourier_all);
            if(SAVE_DATA)
                if(TRIALS_CONCAT_SINGLE_TRIAL_MATCHED_LENGTH)
                    savefig("figures/MI/statisticsCB/EEG_ENVELOPE_DIFFICULTY_OFFSET" + string(start_time * 1000) + "_TRIALS_CONCAT_SINGLETRIAL_MATCHEDLENGTH_" + flow + "TO" + fhigh + "HZ.fig");
                    save("results/MI/statisticsCB/stat_env_difficulty_offset" + string(start_time * 1000) + "_trials_concat_singleTrial_matchedLength_" + flow + "to" + fhigh + "hz.mat", "stat_env");
                else
                    savefig("figures/MI/statisticsCB/EEG_ENVELOPE_DIFFICULTY_OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
                    save("results/MI/statisticsCB/stat_env_difficulty_offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_env");
                end
            end
            figure
            cfg_plot.zlim = [-0.0003 0.0003];
            for pc=PCs
                cfg_plot.highlightchannel=[];
                freq_fourier_all.powspctrm = mean(MI_pc_class1{pc}) - mean(MI_pc_class2{pc});
                cfg_plot.highlightchannel=freq_fourier_all.label(h_MI_pc_np{pc});
                ft_topoplotER(cfg_plot, freq_fourier_all);
                if(SAVE_DATA)
                    if(TRIALS_CONCAT_SINGLE_TRIAL_MATCHED_LENGTH)
                        savefig("figures/MI/statisticsCB/EEG_PC"+ string(pc)+"_DIFFICULTY_OFFSET" + string(start_time * 1000) + "_TRIALS_CONCAT_SINGLETRIAL_MATCHEDLENGTH_" + flow + "TO" + fhigh + "HZ.fig");
                    else
                        savefig("figures/MI/statisticsCB/EEG_PC"+ string(pc)+"_DIFFICULTY_OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
                    end
                end
            end
            if(SAVE_DATA)
                if(TRIALS_CONCAT_SINGLE_TRIAL_MATCHED_LENGTH)
                    save("results/MI/statisticsCB/" + "stat_pcs_difficulty_offset" + string(start_time * 1000) + "_trials_concat_singleTrial_matchedLength_" + flow + "to" + fhigh + "hz.mat", "stat_pc");
                else
                    save("results/MI/statisticsCB/" + "stat_pcs_difficulty_offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_pc");
                end
            end
        end
    end
end

if(TIME_FREQ)
    cd 'results/MI/'
    Env_class1 = zeros(N_SUBJ, 64, length(flows));
    Env_class2 = zeros(N_SUBJ, 64, length(flows));
    PC_class1 = zeros(N_SUBJ, 64, length(flows));
    PC_class2 = zeros(N_SUBJ, 64, length(flows));
    pc = PCs(1);
    for freq=1:length(flows)
        flow = flows(freq);
        fhigh = fhighs(freq);
        fname_env = "MI_ENVELOPE_";
        fname_pcs = "MI_PCs_";
        fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_";
        fname_difficulty = "difficultySplit_";
        fname_bandp = "trials_concat_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";
        Env = load(fname_env + fname_offset + fname_difficulty + fname_bandp).Env;
        PCs = load(fname_pcs + fname_offset + fname_difficulty + fname_bandp).PCs;
        for sub=1:N_SUBJ
            Env_class1(sub, :, freq) = Env{sub}(:, 1)';
            Env_class2(sub, :, freq) = Env{sub}(:, 2)';
            PC_class1(sub, :, freq) = PCs{sub}{pc,1}(:, 1)';
            PC_class2(sub, :, freq) = PCs{sub}{pc,1}(:, 2)';
        end
    end
    %% REMOVING NON EEG CHANNELS
    bad_chs=[22,28,32,41,46];
    freq_fourier_all.label = cellstr(chs_eeg_ok);
    Env_class1(:, bad_chs, :) = [];
    Env_class2(:, bad_chs, :) = [];
    PC_class1(:, bad_chs, :) = [];
    PC_class2(:, bad_chs, :) = [];
    %% CLUSTER-BASED STATISTICS
    cfg=[];
    load('neighbours_ordered')
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan        = 2;
    cfg.neighbours       = neighbours_ordered;  % same as defined for the between-trials experiment
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.05;
    cfg.correcttail = 'alpha';
    cfg.numrandomization = 5000;     % number of draws from the permutation distribution
    design = zeros(2,2*N_SUBJ);
    for k = 1:N_SUBJ
        design(1,k) = k;
    end
    for k = 1:N_SUBJ
        design(1,N_SUBJ+k) = k;
    end
    design(2,1:N_SUBJ)        = 1;
    design(2,N_SUBJ+1:2*N_SUBJ) = 2;
    
    cfg.design = design;
    cfg.uvar   = 1;
    cfg.ivar   = 2;
    
    allSubjEnvClass1.label  = chs_eeg_ok;
    allSubjEnvClass1.dimord = 'subj_chan_freq_time';
    allSubjEnvClass1.time   = 1;
    allSubjEnvClass1.freq = (flows + fhighs) /2;
    allSubjEnvClass1.powspctrm = Env_class1;
    allSubjEnvClass2.label  = chs_eeg_ok;
    allSubjEnvClass2.dimord = 'subj_chan_freq_time';
    allSubjEnvClass2.time   = 1;
    allSubjEnvClass2.freq = (flows + fhighs) /2;
    allSubjEnvClass2.powspctrm = Env_class2;
    
    
    allSubjPCClass1.label  = chs_eeg_ok;
    allSubjPCClass1.dimord = 'subj_chan_freq_time';
    allSubjPCClass1.time   = 1;
    allSubjPCClass1.freq = (flows + fhighs) /2;
    allSubjPCClass1.powspctrm = PC_class1;
    allSubjPCClass2.label  = chs_eeg_ok;
    allSubjPCClass2.dimord = 'subj_chan_freq_time';
    allSubjPCClass2.time   = 1;
    allSubjPCClass2.freq = (flows + fhighs) /2;
    allSubjPCClass2.powspctrm = PC_class2;
    
    %% Statistics
    if(pc == 1)
        [stat_env] = ft_freqstatistics(cfg, allSubjEnvClass1, allSubjEnvClass2)
        h_env_np = stat_env.mask;
    end
    [stat_pc] = ft_freqstatistics(cfg, allSubjPCClass1, allSubjPCClass2)
    h_pc_np = stat_pc.mask;
end
