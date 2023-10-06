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
PCs = [1]
% flows = [0.5:0.5:9.5]
% fhighs = [1.5:0.5:10.5]
flows=[0.5];
fhighs=[4];
N_SUBJ = 22
start_time = 0.5
NUM_ITERATIONS = 1000
COMPUTE_STATISTICS_DIFFICULTY = true
COMPUTE_STATISTICS_PERFORMANCE = false
COMPUTE_STATISTICS_WITHIN_SUBJECTS = true
load 'chs_eeg_ok'
load('freq_template')
bad_chs=[22,28,32,41,46];
freq_fourier_all.label = cellstr(chs_eeg_ok);
load('neighbours_ordered')


NUM_PERFORMANCE = length(load("../behavioral_/results/subs/AM_08" + "_dataSplit").dataSplit.goodPerformanceTrials);
% Performance split to use
for freq=1:length(flows)
    flow = flows(freq)
    fhigh = fhighs(freq)
    %% FILE LOADING
    cd 'results/PID/'
    fname_redundancy = "REDUNDANCY_";
    fname_unq1 = "UNIQUE1_";
    fname_unq2 = "UNIQUE2_";
    fname_synergy = "SYNERGY_";
    fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_";
    fname_difficulty = "difficultySplit_matchedLength_1000_reps_";
    fname_performance = "performanceSplit_";
    fname_bandp = "trials_concat_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";

    % Difficulty
    if(COMPUTE_STATISTICS_DIFFICULTY)
        fname_red_difficulty = fname_redundancy + fname_offset + fname_difficulty + fname_bandp;
        fname_unq1_difficulty = fname_unq1 + fname_offset + fname_difficulty + fname_bandp;
        fname_unq2_difficulty = fname_unq2 + fname_offset + fname_difficulty + fname_bandp;
        fname_syn_difficulty = fname_synergy + fname_offset + fname_difficulty + fname_bandp;

        Rdn = {load(fname_red_difficulty).Rdn};
        Unq1 = {load(fname_unq1_difficulty).Unq1};
        Unq2 = {load(fname_unq2_difficulty).Unq2};
        Syn = {load(fname_syn_difficulty).Syn};
    end
    % Performance
    if(COMPUTE_STATISTICS_PERFORMANCE)
        for n_performance=1:NUM_PERFORMANCE
            if(performances(n_performance))
                fname_red_performance = fname_redundancy + fname_offset + fname_performance + string(n_performance) + "_" + fname_bandp;
                fname_unq1_performance = fname_unq1 + fname_offset + fname_performance + string(n_performance) + "_" + fname_bandp;
                fname_unq2_performance = fname_unq2 + fname_offset + fname_performance + string(n_performance) + "_" + fname_bandp;
                fname_syn_performance = fname_synergy + fname_offset + fname_performance + string(n_performance) + "_" + fname_bandp;

                Rdn{n_performance} = load(fname_red_performance).Rdn;
                Unq1{n_performance} = load(fname_unq1_performance).Unq1;
                Unq2{n_performance} = load(fname_unq2_performance).Unq2;
                Syn{n_performance} = load(fname_syn_performance).Syn;
            end
        end
    end
    %% PREPARING PID DATA
    pos_clusters_rdn = zeros(1, NUM_ITERATIONS);
    neg_clusters_rdn = zeros(1, NUM_ITERATIONS);
    pos_clusters_unq1 = zeros(1, NUM_ITERATIONS);
    neg_clusters_unq1 = zeros(1, NUM_ITERATIONS);
    pos_clusters_unq2 = zeros(1, NUM_ITERATIONS);
    neg_clusters_unq2 = zeros(1, NUM_ITERATIONS);
    pos_clusters_syn = zeros(1, NUM_ITERATIONS);
    neg_clusters_syn = zeros(1, NUM_ITERATIONS);
    parfor iter=1:NUM_ITERATIONS
        disp("Starting with iteration: " + iter);
        for i =1:length(Rdn)
            for pc=PCs
                Js = [1];
                for j=Js
                    Rdn_class1 = {};
                    Rdn_class2 = {};
                    Unq1_class1 = {};
                    Unq1_class2 = {};
                    Unq2_class1 = {};
                    Unq2_class2 = {};
                    Syn_class1 = {};
                    Syn_class2 = {};
                    % Difficulty
                    Rdn_class1{i} = zeros(N_SUBJ, 64);
                    Rdn_class2{i} = zeros(N_SUBJ, 64);
                    Unq1_class1{i} = zeros(N_SUBJ, 64);
                    Unq1_class2{i} = zeros(N_SUBJ, 64);
                    Unq2_class1{i} = zeros(N_SUBJ, 64);
                    Unq2_class2{i} = zeros(N_SUBJ, 64);
                    Syn_class1{i} = zeros(N_SUBJ, 64);
                    Syn_class2{i} = zeros(N_SUBJ, 64);
                    for sub=1:N_SUBJ
                        Rdn_class1{i}(sub, :) = Rdn{i}{sub}{iter}{pc, j}(:, 1)';
                        Rdn_class2{i}(sub, :) = Rdn{i}{sub}{iter}{pc, j}(:, 2)';
                        Unq1_class1{i}(sub, :) = Unq1{i}{sub}{iter}{pc, j}(:, 1)';
                        Unq1_class2{i}(sub, :) = Unq1{i}{sub}{iter}{pc, j}(:, 2)';
                        Unq2_class1{i}(sub, :) = Unq2{i}{sub}{iter}{pc, j}(:, 1)';
                        Unq2_class2{i}(sub, :) = Unq2{i}{sub}{iter}{pc, j}(:, 2)';
                        Syn_class1{i}(sub, :) = Syn{i}{sub}{iter}{pc, j}(:, 1)';
                        Syn_class2{i}(sub, :) = Syn{i}{sub}{iter}{pc, j}(:, 2)';
                    end

                    %% REMOVING NON EEG CHANNELS
                    % Difficulty
                    Rdn_class1{i}(:, bad_chs) = [];
                    Rdn_class2{i}(:, bad_chs) = [];
                    Unq1_class1{i}(:, bad_chs) = [];
                    Unq1_class2{i}(:, bad_chs) = [];
                    Unq2_class1{i}(:, bad_chs) = [];
                    Unq2_class2{i}(:, bad_chs) = [];
                    Syn_class1{i}(:, bad_chs) = [];
                    Syn_class2{i}(:, bad_chs) = [];
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
                        
                        allSubjRdnClass1 = [];
                        allSubjRdnClass2 = [];
                        allSubjRdnClass1.label  = chs_eeg_ok;
                        allSubjRdnClass1.dimord = 'subj_chan_time';
                        allSubjRdnClass1.time   = 1;
                        allSubjRdnClass1.avg = Rdn_class1{i};
                        allSubjRdnClass2.label  = chs_eeg_ok;
                        allSubjRdnClass2.dimord = 'subj_chan_time';
                        allSubjRdnClass2.time   = 1;
                        allSubjRdnClass2.avg = Rdn_class2{i};

                        allSubjUnq1Class1 = [];
                        allSubjUnq1Class2 = [];
                        allSubjUnq1Class1.label  = chs_eeg_ok;
                        allSubjUnq1Class1.dimord = 'subj_chan_time';
                        allSubjUnq1Class1.time   = 1;
                        allSubjUnq1Class1.avg = Unq1_class1{i};
                        allSubjUnq1Class2.label  = chs_eeg_ok;
                        allSubjUnq1Class2.dimord = 'subj_chan_time';
                        allSubjUnq1Class2.time   = 1;
                        allSubjUnq1Class2.avg = Unq1_class2{i};

                        allSubjUnq2Class1 = [];
                        allSubjUnq2Class2 = [];
                        allSubjUnq2Class1.label  = chs_eeg_ok;
                        allSubjUnq2Class1.dimord = 'subj_chan_time';
                        allSubjUnq2Class1.time   = 1;
                        allSubjUnq2Class1.avg = Unq2_class1{i};
                        allSubjUnq2Class2.label  = chs_eeg_ok;
                        allSubjUnq2Class2.dimord = 'subj_chan_time';
                        allSubjUnq2Class2.time   = 1;
                        allSubjUnq2Class2.avg = Unq2_class2{i};
                        
                        allSubjSynClass1 = [];
                        allSubjSynClass2 = [];
                        allSubjSynClass1.label  = chs_eeg_ok;
                        allSubjSynClass1.dimord = 'subj_chan_time';
                        allSubjSynClass1.time   = 1;
                        allSubjSynClass1.avg = Syn_class1{i};
                        allSubjSynClass2.label  = chs_eeg_ok;
                        allSubjSynClass2.dimord = 'subj_chan_time';
                        allSubjSynClass2.time   = 1;
                        allSubjSynClass2.avg = Syn_class2{i};
                        %% Statistics
                        [stat_rdn] = ft_timelockstatistics(cfg, allSubjRdnClass1, allSubjRdnClass2)
                        [stat_unq1] = ft_timelockstatistics(cfg, allSubjUnq1Class1, allSubjUnq1Class2)
                        [stat_unq2] = ft_timelockstatistics(cfg, allSubjUnq2Class1, allSubjUnq2Class2)
                        [stat_syn] = ft_timelockstatistics(cfg, allSubjSynClass1, allSubjSynClass2)
                        % update rdn clusters
                        if(isfield(stat_rdn, 'posclusters'))
                            if(~isempty(stat_rdn.posclusters))
                                for cluster=1:length(stat_rdn.posclusters)
                                    if(stat_rdn.posclusters(cluster).prob <= 0.025)
                                        pos_clusters_rdn(1,iter) = pos_clusters_rdn(1,iter) + 1;
                                    end
                                end
                            end
                        end
                        if(isfield(stat_rdn, 'negclusters'))
                            if(~isempty(stat_rdn.negclusters))
                                for cluster=1:length(stat_rdn.negclusters)
                                    if(stat_rdn.negclusters(cluster).prob <= 0.025)
                                        neg_clusters_rdn(1,iter) = neg_clusters_rdn(1,iter) + 1;
                                    end
                                end
                            end
                        end
                        if(isfield(stat_unq1, 'posclusters'))
                            if(~isempty(stat_unq1.posclusters))
                                for cluster=1:length(stat_unq1.posclusters)
                                    if(stat_unq1.posclusters(cluster).prob <= 0.025)
                                        pos_clusters_unq1(1,iter) = pos_clusters_unq1(1,iter) + 1;
                                    end
                                end
                            end
                        end
                        if(isfield(stat_unq1, 'negclusters'))
                            if(~isempty(stat_unq1.negclusters))
                                for cluster=1:length(stat_unq1.negclusters)
                                    if(stat_unq1.negclusters(cluster).prob <= 0.025)
                                        neg_clusters_unq1(1,iter) = neg_clusters_unq1(1,iter) + 1;
                                    end
                                end
                            end
                        end
                        if(isfield(stat_unq2, 'posclusters'))
                            if(~isempty(stat_unq2.posclusters))
                                for cluster=1:length(stat_unq2.posclusters)
                                    if(stat_unq2.posclusters(cluster).prob <= 0.025)
                                        pos_clusters_unq2(1,iter) = pos_clusters_unq2(1,iter) + 1;
                                    end
                                end
                            end
                        end
                        if(isfield(stat_unq2, 'negclusters'))
                            if(~isempty(stat_unq2.negclusters))
                                for cluster=1:length(stat_unq2.negclusters)
                                    if(stat_unq2.negclusters(cluster).prob <= 0.025)
                                        neg_clusters_unq2(1,iter) = neg_clusters_unq2(1,iter) + 1;
                                    end
                                end
                            end
                        end
                        if(isfield(stat_syn, 'posclusters'))
                            if(~isempty(stat_syn.posclusters))
                                for cluster=1:length(stat_syn.posclusters)
                                    if(stat_syn.posclusters(cluster).prob <= 0.025)
                                        pos_clusters_syn(1,iter) = pos_clusters_syn(1,iter) + 1;
                                    end
                                end
                            end
                        end
                        if(isfield(stat_syn, 'negclusters'))
                            if(~isempty(stat_syn.negclusters))
                                for cluster=1:length(stat_syn.negclusters)
                                    if(stat_syn.negclusters(cluster).prob <= 0.025)
                                        neg_clusters_syn(1,iter) = neg_clusters_syn(1,iter) + 1;
                                    end
                                end
                            end
                        end
                    end
                    if(not(pc == PCs(end)))
                        cd 'results/PID/'
                    end
                end
            end
        end
    end
    cd 'statisticsCB/PC1'
    percentages_pos_clusters = zeros(1,4);
    percentages_neg_clusters = zeros(1,4);
    percentages_pos_clusters(1) = sum(pos_clusters_rdn(1, :) > 0)/NUM_ITERATIONS;
    percentages_neg_clusters(1) = sum(neg_clusters_rdn(1,:) > 0)/NUM_ITERATIONS;
    percentages_pos_clusters(2) = sum(pos_clusters_unq1(1, :) > 0)/NUM_ITERATIONS;
    percentages_neg_clusters(2) = sum(neg_clusters_unq1(1,:) > 0)/NUM_ITERATIONS;
    percentages_pos_clusters(3) = sum(pos_clusters_unq2(1, :) > 0)/NUM_ITERATIONS;
    percentages_neg_clusters(3) = sum(neg_clusters_unq2(1,:) > 0)/NUM_ITERATIONS;
    percentages_pos_clusters(4) = sum(pos_clusters_syn(1, :) > 0)/NUM_ITERATIONS;
    percentages_neg_clusters(4) = sum(neg_clusters_syn(1,:) > 0)/NUM_ITERATIONS;
    save('negative_clusters_' + string(NUM_ITERATIONS) + '_reps_' + flows(freq) + '_' + fhighs(freq) + '.mat', 'percentages_neg_clusters')
    save('positive_clusters_' + string(NUM_ITERATIONS) + '_reps_' + flows(freq) + '_' + fhighs(freq) + '.mat', 'percentages_pos_clusters')
    cd '../../'
end