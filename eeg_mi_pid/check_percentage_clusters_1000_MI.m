clear
close all
% cd 'D:\Entrainment__\final_functions_behavior\'

%% ADDING PATHS AND SETTING PARAMS
addpath '../fieldtrip-20210311'
addpath './results'
addpath '../analisi'
cmap = brewermap(1024,'RdBu');
cmap = flip(cmap);


%%% PARAMETERS %%%
PC = [1]
flows=[0.5];
fhighs=[4];
N_SUBJ = 22
start_time = 0.5
NUM_ITERATIONS = 1000
COMPUTE_STATISTICS_DIFFICULTY = true
COMPUTE_STATISTICS_WITHIN_SUBJECTS = true
load 'chs_eeg_ok'
load('freq_template')
bad_chs=[22,28,32,41,46];
freq_fourier_all.label = cellstr(chs_eeg_ok);
load('neighbours_ordered')


% Performance split to use
for freq=1:length(flows)
    flow = flows(freq)
    fhigh = fhighs(freq)
    %% FILE LOADING
    cd 'results/MI/'
    fname_Env = "MI_ENVELOPE_";
    fname_PCs = "MI_PCs_";
    fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_";
    fname_difficulty = "difficultySplit_matchedLength_1000_reps_";
    fname_performance = "performanceSplit_";
    fname_bandp = "trials_concat_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";

    % Difficulty
    if(COMPUTE_STATISTICS_DIFFICULTY)
        fname_env_difficulty = fname_Env + fname_offset + fname_difficulty + fname_bandp;
        fname_pcs_difficulty = fname_PCs + fname_offset + fname_difficulty + fname_bandp;

        Env = {load(fname_env_difficulty).Env};
        PCs = {load(fname_pcs_difficulty).PCs};
    end

    %% PREPARING PID DATA
    pos_clusters_env = zeros(1, NUM_ITERATIONS);
    neg_clusters_env = zeros(1, NUM_ITERATIONS);
    pos_clusters_pcs = zeros(1, NUM_ITERATIONS);
    neg_clusters_pcs = zeros(1, NUM_ITERATIONS);

    parfor iter=1:NUM_ITERATIONS
        disp("Starting with iteration: " + iter);
        for i =1:length(Env)
            for pc=PC
                Js = [1];
                for j=Js
                    Env_class1 = {};
                    Env_class2 = {};
                    PCs_class1 = {};
                    PCs_class2 = {};

                    % Difficulty
                    Env_class1{i} = zeros(N_SUBJ, 64);
                    Env_class2{i} = zeros(N_SUBJ, 64);
                    PCs_class1{i} = zeros(N_SUBJ, 64);
                    PCs_class2{i} = zeros(N_SUBJ, 64);

                    for sub=1:N_SUBJ
                        Env_class1{i}(sub, :) = Env{i}{sub}{iter}(:, 1)';
                        Env_class2{i}(sub, :) = Env{i}{sub}{iter}(:, 2)';
                        PCs_class1{i}(sub, :) = PCs{i}{sub}{iter}{pc, j}(:, 1)';
                        PCs_class2{i}(sub, :) = PCs{i}{sub}{iter}{pc, j}(:, 2)';
                    end

                    %% REMOVING NON EEG CHANNELS
                    % Difficulty
                    Env_class1{i}(:, bad_chs) = [];
                    Env_class2{i}(:, bad_chs) = [];
                    PCs_class1{i}(:, bad_chs) = [];
                    PCs_class2{i}(:, bad_chs) = [];

                    %% CLUSTER-BASED STATISTICS WITHIN SUBJECTS
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
                        design = zeros(2, (N_SUBJ)*2);
                        design(1,:) = [1:N_SUBJ 1:N_SUBJ];
                        design(2,:) = [ones(1,N_SUBJ) ones(1,N_SUBJ)*2];

                        cfg.design = design;
                        cfg.uvar   = 1;
                        cfg.ivar   = 2;
                        
                        allSubjEnvClass1 = [];
                        allSubjEnvClass2 = [];
                        allSubjEnvClass1.label  = chs_eeg_ok;
                        allSubjEnvClass1.dimord = 'subj_chan_time';
                        allSubjEnvClass1.time   = 1;
                        allSubjEnvClass1.avg = Env_class1{i};
                        allSubjEnvClass2.label  = chs_eeg_ok;
                        allSubjEnvClass2.dimord = 'subj_chan_time';
                        allSubjEnvClass2.time   = 1;
                        allSubjEnvClass2.avg = Env_class2{i};

                        allSubjPCsClass1 = [];
                        allSubjPCsClass2 = [];
                        allSubjPCsClass1.label  = chs_eeg_ok;
                        allSubjPCsClass1.dimord = 'subj_chan_time';
                        allSubjPCsClass1.time   = 1;
                        allSubjPCsClass1.avg = PCs_class1{i};
                        allSubjPCsClass2.label  = chs_eeg_ok;
                        allSubjPCsClass2.dimord = 'subj_chan_time';
                        allSubjPCsClass2.time   = 1;
                        allSubjPCsClass2.avg = PCs_class2{i};


                        %% Statistics
                        [stat_env] = ft_timelockstatistics(cfg, allSubjEnvClass1, allSubjEnvClass2)
                        [stat_pcs] = ft_timelockstatistics(cfg, allSubjPCsClass1, allSubjPCsClass2)

                        if(isfield(stat_env, 'posclusters'))
                            if(~isempty(stat_env.posclusters))
                                for cluster=1:length(stat_env.posclusters)
                                    if(stat_env.posclusters(cluster).prob <= 0.025)
                                        pos_clusters_env(1,iter) = pos_clusters_env(1,iter) + 1;
                                    end
                                end
                            end
                        end
                        if(isfield(stat_env, 'negclusters'))
                            if(~isempty(stat_env.negclusters))
                                for cluster=1:length(stat_env.negclusters)
                                    if(stat_env.negclusters(cluster).prob <= 0.025)
                                        neg_clusters_env(1,iter) = neg_clusters_env(1,iter) + 1;
                                    end
                                end
                            end
                        end
                        if(isfield(stat_pcs, 'posclusters'))
                            if(~isempty(stat_pcs.posclusters))
                                for cluster=1:length(stat_pcs.posclusters)
                                    if(stat_pcs.posclusters(cluster).prob <= 0.025)
                                        pos_clusters_pcs(1,iter) = pos_clusters_pcs(1,iter) + 1;
                                    end
                                end
                            end
                        end
                        if(isfield(stat_pcs, 'negclusters'))
                            if(~isempty(stat_pcs.negclusters))
                                for cluster=1:length(stat_pcs.negclusters)
                                    if(stat_pcs.negclusters(cluster).prob <= 0.025)
                                        neg_clusters_pcs(1,iter) = neg_clusters_pcs(1,iter) + 1;
                                    end
                                end
                            end
                        end
                    end
                    if(not(pc == PC(end)))
                        cd 'results/MI/'
                    end
                end
            end
        end
    end
    cd 'statisticsCB/'
    percentages_pos_clusters = zeros(1,4);
    percentages_neg_clusters = zeros(1,4);
    percentages_pos_clusters(1) = sum(pos_clusters_env(1, :) > 0)/NUM_ITERATIONS;
    percentages_neg_clusters(1) = sum(neg_clusters_env(1,:) > 0)/NUM_ITERATIONS;
    percentages_pos_clusters(2) = sum(pos_clusters_pcs(1, :) > 0)/NUM_ITERATIONS;
    percentages_neg_clusters(2) = sum(neg_clusters_pcs(1,:) > 0)/NUM_ITERATIONS;

    save('negative_clusters_' + string(NUM_ITERATIONS) + '_reps_' + flows(freq) + '_' + fhighs(freq) + '.mat', 'percentages_neg_clusters')
    save('positive_clusters_' + string(NUM_ITERATIONS) + '_reps_' + flows(freq) + '_' + fhighs(freq) + '.mat', 'percentages_pos_clusters')
    cd '../../'
end