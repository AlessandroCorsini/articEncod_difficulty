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
flows = [0.5:0.5:9.5];
fhighs = [1.5:0.5:10.5];
flows = [0.5];
fhighs = [4];
lags = [0.2];
lag_eeg = 0.3;
lag_pc = 0.1;
lags_env = [0.2];
N_SUBJ = 22
start_time = 0.5
COMPUTE_STATISTICS_DIFFICULTY = true
COMPUTE_STATISTICS_PERFORMANCE = false
COMPUTE_STATISTICS_WITHIN_SUBJECTS = true
COMPUTE_STATISTICS_BETWEEN_PCS = false
TIME = true;
TIME_FREQ = false;
load 'chs_eeg_ok'
load('freq_template')

NUM_PERFORMANCE = length(load("../behavioral_/results/subs/AM_08" + "_dataSplit").dataSplit.goodPerformanceTrials);
% Performance split to use
if(TIME)
    for freq=1:length(flows)
        flow = flows(freq)
        fhigh = fhighs(freq)
                cd 'results/PID/'
        for lag=lags_env
        %% FILE LOADING

        if(COMPUTE_STATISTICS_BETWEEN_PCS)
            fname_redundancy = "REDUNDANCY_PCs_";
            fname_unq1 = "UNIQUE1_PCs_";
            fname_unq2 = "UNIQUE2_PCs_";
            fname_synergy = "SYNERGY_PCs_";
        else
            fname_redundancy = "REDUNDANCY_";
            fname_unq1 = "UNIQUE1_";
            fname_unq2 = "UNIQUE2_";
            fname_synergy = "SYNERGY_";
        end
        % Levare matchedLength_ per non usare i matched
        fname_offset = "offset" + string(start_time * 1000) + "ms_lag" + string(lag*1000) + "ms_allFreqs_";
        %fname_offset = "offset" + string(start_time * 1000) + "ms_lag_eeg" + string(lag_eeg*1000) + "ms_lag_pc" + string(lag_pc*1000) + "ms_lag_env" + string(lags_env(1)*1000) + "ms_allFreqs_";
        %fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_";
        fname_difficulty = "difficultySplit_";
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
        % levare il for
        %for iter=145
        for i =1:length(Rdn)
            for pc=PCs
                Js = [1];
                if(COMPUTE_STATISTICS_BETWEEN_PCS)
                    Js = PCs(PCs > pc);
                end
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
% Levare {iter} per non usare i matchedLength
                    for sub=1:N_SUBJ
                        Rdn_class1{i}(sub, :) = Rdn{i}{sub}{pc, j}(:, 1)';
                        Rdn_class2{i}(sub, :) = Rdn{i}{sub}{pc, j}(:, 2)';
                        Unq1_class1{i}(sub, :) = Unq1{i}{sub}{pc, j}(:, 1)';
                        Unq1_class2{i}(sub, :) = Unq1{i}{sub}{pc, j}(:, 2)';
                        Unq2_class1{i}(sub, :) = Unq2{i}{sub}{pc, j}(:, 1)';
                        Unq2_class2{i}(sub, :) = Unq2{i}{sub}{pc, j}(:, 2)';
                        Syn_class1{i}(sub, :) = Syn{i}{sub}{pc, j}(:, 1)';
                        Syn_class2{i}(sub, :) = Syn{i}{sub}{pc, j}(:, 2)';
                    end

                    %% REMOVING NON EEG CHANNELS
                    bad_chs=[22,28,32,41,46];
                    freq_fourier_all.label = cellstr(chs_eeg_ok);
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
                        load('neighbours_ordered')
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

                        allSubjRdnClass1.label  = chs_eeg_ok;
                        allSubjRdnClass1.dimord = 'subj_chan_time';
                        allSubjRdnClass1.time   = 1;
                        allSubjRdnClass1.avg = Rdn_class1{i};
                        allSubjRdnClass2.label  = chs_eeg_ok;
                        allSubjRdnClass2.dimord = 'subj_chan_time';
                        allSubjRdnClass2.time   = 1;
                        allSubjRdnClass2.avg = Rdn_class2{i};


                        allSubjUnq1Class1.label  = chs_eeg_ok;
                        allSubjUnq1Class1.dimord = 'subj_chan_time';
                        allSubjUnq1Class1.time   = 1;
                        allSubjUnq1Class1.avg = Unq1_class1{i};
                        allSubjUnq1Class2.label  = chs_eeg_ok;
                        allSubjUnq1Class2.dimord = 'subj_chan_time';
                        allSubjUnq1Class2.time   = 1;
                        allSubjUnq1Class2.avg = Unq1_class2{i};


                        allSubjUnq2Class1.label  = chs_eeg_ok;
                        allSubjUnq2Class1.dimord = 'subj_chan_time';
                        allSubjUnq2Class1.time   = 1;
                        allSubjUnq2Class1.avg = Unq2_class1{i};
                        allSubjUnq2Class2.label  = chs_eeg_ok;
                        allSubjUnq2Class2.dimord = 'subj_chan_time';
                        allSubjUnq2Class2.time   = 1;
                        allSubjUnq2Class2.avg = Unq2_class2{i};

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
                        h_rdn_np=stat_rdn.mask;
                        [stat_unq1] = ft_timelockstatistics(cfg, allSubjUnq1Class1, allSubjUnq1Class2)
                        h_unq1_np=stat_unq1.mask;
                        [stat_unq2] = ft_timelockstatistics(cfg, allSubjUnq2Class1, allSubjUnq2Class2)
                        h_unq2_np=stat_unq2.mask;
                        [stat_syn] = ft_timelockstatistics(cfg, allSubjSynClass1, allSubjSynClass2)
                        h_syn_np=stat_syn.mask;
                        %% Topoplot
%                         try
                            cfg_plot=[];
                            cfg_plot.colormap=cmap;
                            cfg_plot.marker =  'off'
                            cfg_plot.highlight='on'
                            cfg_plot.highlightchannel=[]
                            cfg_plot.highlightsymbol='.';
                            cfg_plot.highlightcolor='k';
                            cfg_plot.highlightsize=18;
                            freq_fourier_all.freq=[1];
%                             cd D:/Entrainment__/final_functions_behavior\
% levare matchedlenght da fname fig
                            fname_fig = "";%"matchedLength_" + iter + "_";
                            if(COMPUTE_STATISTICS_BETWEEN_PCS)
                                fname_fig = "PC" + string(pc) + "_" + string(j) + "_";
                            end
                            cfg_plot.zlim = [-0.00015, 0.00015];
                            freq_fourier_all.powspctrm = mean(Rdn_class1{i}) - mean(Rdn_class2{i});
%                             ft_topoplotER(cfg_plot, freq_fourier_all);
%                             title("REDUNDANCY DIFFICULTY " + pc + "-" + j)
                            cfg_plot.highlightchannel = freq_fourier_all.label(h_rdn_np)
                            ft_topoplotER(cfg_plot, freq_fourier_all);
%                            savefig("../../figures/PID/statisticsCB/PC" + string(pc) + "/REDUNDANCY_DIFFICULTY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
%                            save("../../results/PID/statisticsCB/PC" + string(pc) + "/stat_rdn_difficulty_" + fname_fig + "offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_rdn");
                            
                            cfg_plot.zlim = [-0.0005, 0.0005];
                            figure
                            cfg_plot.highlightchannel=[]
                            freq_fourier_all.powspctrm = mean(Unq1_class1{i}) - mean(Unq1_class2{i});
                            ft_topoplotER(cfg_plot, freq_fourier_all);
                            title("UNIQUE1 SE")
                            cfg_plot.highlightchannel = freq_fourier_all.label(h_unq1_np)
                            ft_topoplotER(cfg_plot, freq_fourier_all);
                         % savefig("../../figures/PID/statisticsCB/PC" + string(pc) + "/UNIQUE1_DIFFICULTY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
                            % save("../../results/PID/statisticsCB/PC" + string(pc) + "/stat_unq1_difficulty_" + fname_fig + "offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_unq1");
                            
                            cfg_plot.zlim = [-0.0005, 0.0005];
                            figure
                            cfg_plot.highlightchannel=[]
                            freq_fourier_all.powspctrm = mean(Unq2_class1{i}) - mean(Unq2_class2{i});
                            ft_topoplotER(cfg_plot, freq_fourier_all);
                            title("UNIQUE PC1")
                            cfg_plot.highlightchannel = freq_fourier_all.label(h_unq2_np)
                            ft_topoplotER(cfg_plot, freq_fourier_all);
                            % savefig("../../figures/PID/statisticsCB/PC" + string(pc) + "/UNIQUE2_DIFFICULTY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
                           % save("../../results/PID/statisticsCB/PC" + string(pc) + "/stat_unq2_difficulty_" + fname_fig + "offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_unq2");

                            cfg_plot.zlim = [-0.0001, 0.0001];
                            figure
                            cfg_plot.highlightchannel=[]
                            freq_fourier_all.powspctrm = mean(Syn_class1{i}) - mean(Syn_class2{i});
%                             ft_topoplotER(cfg_plot, freq_fourier_all);
                            title("SYNERGY DIFFICULTY " + pc + "-" + j)
                            cfg_plot.highlightchannel = freq_fourier_all.label(h_syn_np)
                            ft_topoplotER(cfg_plot,freq_fourier_all);
%                             savefig("../../figures/PID/statisticsCB/PC" + string(pc) + "/SYNERGY_DIFFICULTY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
%                             save("../../results/PID/statisticsCB/PC" + string(pc) + "/stat_syn_difficulty_" + fname_fig + "offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_syn");
%                         end
                    end
                    if(not(pc == PCs(end)))
%                         cd 'results/PID/'
                    end
                end
            end
        end
        end
    end
end

NUM_SPLIT = 1;
if(TIME_FREQ)
    %% PREPARING PID DATA
    for i =1:NUM_SPLIT
        for pc=PCs
            cd 'results/PID/'
            Rdn_class1 = {};
            Rdn_class2 = {};
            Unq1_class1 = {};
            Unq1_class2 = {};
            Unq2_class1 = {};
            Unq2_class2 = {};
            Syn_class1 = {};
            Syn_class2 = {};
            % Difficulty
            Rdn_class1{i} = zeros(N_SUBJ, 64, length(flows));
            Rdn_class2{i} = zeros(N_SUBJ, 64, length(flows));
            Unq1_class1{i} = zeros(N_SUBJ, 64, length(flows));
            Unq1_class2{i} = zeros(N_SUBJ, 64, length(flows));
            Unq2_class1{i} = zeros(N_SUBJ, 64, length(flows));
            Unq2_class2{i} = zeros(N_SUBJ, 64, length(flows));
            Syn_class1{i} = zeros(N_SUBJ, 64, length(flows));
            Syn_class2{i} = zeros(N_SUBJ, 64, length(flows));
            Js = [1];
            if(COMPUTE_STATISTICS_BETWEEN_PCS)
                Js = PCs(PCs > pc);
            end
            for j=Js
                for freq=1:length(flows)
                    flow = flows(freq)
                    fhigh = fhighs(freq)
                    %% FILE LOADING
                    if(COMPUTE_STATISTICS_BETWEEN_PCS)
                        fname_redundancy = "REDUNDANCY_PCs_";
                        fname_unq1 = "UNIQUE1_PCs_";
                        fname_unq2 = "UNIQUE2_PCs_";
                        fname_synergy = "SYNERGY_PCs_";
                    else
                        fname_redundancy = "REDUNDANCY_";
                        fname_unq1 = "UNIQUE1_";
                        fname_unq2 = "UNIQUE2_";
                        fname_synergy = "SYNERGY_";
                    end
                    fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_";
                    fname_difficulty = "difficultySplit_";
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

                    for sub=1:N_SUBJ
                        Rdn_class1{i}(sub, :, freq) = Rdn{i}{sub}{pc, j}(:, 1)';
                        Rdn_class2{i}(sub, :, freq) = Rdn{i}{sub}{pc, j}(:, 2)';
                        Unq1_class1{i}(sub, :, freq) = Unq1{i}{sub}{pc, j}(:, 1)';
                        Unq1_class2{i}(sub, :, freq) = Unq1{i}{sub}{pc, j}(:, 2)';
                        Unq2_class1{i}(sub, :, freq) = Unq2{i}{sub}{pc, j}(:, 1)';
                        Unq2_class2{i}(sub, :, freq) = Unq2{i}{sub}{pc, j}(:, 2)';
                        Syn_class1{i}(sub, :, freq) = Syn{i}{sub}{pc, j}(:, 1)';
                        Syn_class2{i}(sub, :, freq) = Syn{i}{sub}{pc, j}(:, 2)';
                    end
                end
                %% REMOVING NON EEG CHANNELS
                bad_chs=[22,28,32,41,46];
                freq_fourier_all.label = cellstr(chs_eeg_ok);
                % Difficulty
                Rdn_class1{i}(:, bad_chs, :) = [];
                Rdn_class2{i}(:, bad_chs, :) = [];
                Unq1_class1{i}(:, bad_chs, :) = [];
                Unq1_class2{i}(:, bad_chs, :) = [];
                Unq2_class1{i}(:, bad_chs, :) = [];
                Unq2_class2{i}(:, bad_chs, :) = [];
                Syn_class1{i}(:, bad_chs, :) = [];
                Syn_class2{i}(:, bad_chs, :) = [];
                %% CLUSTER-BASED STATISTICS WITHIN SUBJECTS

                if COMPUTE_STATISTICS_WITHIN_SUBJECTS

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

                    allSubjRdnClass1.label  = chs_eeg_ok;
                    allSubjRdnClass1.dimord = 'subj_chan_freq_time';
                    allSubjRdnClass1.time   = 1;
                    allSubjRdnClass1.freq = (flows + fhighs) /2;
                    allSubjRdnClass1.powspctrm = Rdn_class1{i};
                    allSubjRdnClass2.label  = chs_eeg_ok;
                    allSubjRdnClass2.dimord = 'subj_chan_freq_time';
                    allSubjRdnClass2.time   = 1;
                    allSubjRdnClass2.freq = (flows + fhighs) /2;
                    allSubjRdnClass2.powspctrm = Rdn_class2{i};


                    allSubjUnq1Class1.label  = chs_eeg_ok;
                    allSubjUnq1Class1.dimord = 'subj_chan_freq_time';
                    allSubjUnq1Class1.time   = 1;
                    allSubjUnq1Class1.freq = (flows + fhighs) /2;
                    allSubjUnq1Class1.powspctrm = Unq1_class1{i};
                    allSubjUnq1Class2.label  = chs_eeg_ok;
                    allSubjUnq1Class2.dimord = 'subj_chan_freq_time';
                    allSubjUnq1Class2.time   = 1;
                    allSubjUnq1Class2.freq = (flows + fhighs) /2;
                    allSubjUnq1Class2.powspctrm = Unq1_class2{i};


                    allSubjUnq2Class1.label  = chs_eeg_ok;
                    allSubjUnq2Class1.dimord = 'subj_chan_freq_time';
                    allSubjUnq2Class1.time   = 1;
                    allSubjUnq2Class1.freq = (flows + fhighs) /2;
                    allSubjUnq2Class1.powspctrm = Unq2_class1{i};
                    allSubjUnq2Class2.label  = chs_eeg_ok;
                    allSubjUnq2Class2.dimord = 'subj_chan_freq_time';
                    allSubjUnq2Class2.time   = 1;
                    allSubjUnq2Class2.freq = (flows + fhighs) /2;
                    allSubjUnq2Class2.powspctrm = Unq2_class2{i};

                    allSubjSynClass1.label  = chs_eeg_ok;
                    allSubjSynClass1.dimord = 'subj_chan_freq_time';
                    allSubjSynClass1.time   = 1;
                    allSubjSynClass1.freq = (flows + fhighs) /2;
                    allSubjSynClass1.powspctrm = Syn_class1{i};
                    allSubjSynClass2.label  = chs_eeg_ok;
                    allSubjSynClass2.dimord = 'subj_chan_freq_time';
                    allSubjSynClass2.time   = 1;
                    allSubjSynClass2.freq = (flows + fhighs) /2;
                    allSubjSynClass2.powspctrm = Syn_class2{i};
                    %% Statistics
                    [stat_rdn] = ft_freqstatistics(cfg, allSubjRdnClass1, allSubjRdnClass2)
                    h_rdn_np=stat_rdn.mask;
                    [stat_unq1] = ft_freqstatistics(cfg, allSubjUnq1Class1, allSubjUnq1Class2)
                    h_unq1_np=stat_unq1.mask;
                    [stat_unq2] = ft_freqstatistics(cfg, allSubjUnq2Class1, allSubjUnq2Class2)
                    h_unq2_np=stat_unq2.mask;
                    [stat_syn] = ft_freqstatistics(cfg, allSubjSynClass1, allSubjSynClass2)
                    h_syn_np=stat_syn.mask;
                    %% Topoplot
                    %try
                        figure
                        cfg_plot = [];
                        cfg_plot.alpha  = 0.025;
                        cfg_plot.colormap=cmap;
                        %cfg_plot.marker =  'off'
%                         cfg_plot.highlight='on'
%                         %cfg_plot.highlightchannel=[]
%                         cfg_plot.highlightsymbol='.';
%                         cfg_plot.highlightcolor='k';
%                         cfg_plot.highlightsize=18;

                        %cd D:/Entrainment__/final_functions_behavior/

                        fname_fig = "";
                        if(COMPUTE_STATISTICS_BETWEEN_PCS)
                            fname_fig = "PC" + string(pc) + "_" + string(j) + "_";
                        end
                        cfg_plot.zlim = [-0.00015, 0.00015];
                        %                         stat_rdn.raweffect = mean(Rdn_class1{i},1) - mean(Rdn_class2{i},1);
                        cfg_plot.parameter = 'stat_rdn';
                        title("REDUNDANCY DIFFICULTY " + pc + "-" + j)
                        %cfg_plot.highlightchannel = stat_rdn.label(h_rdn_np)
                        ft_clusterplot(cfg_plot, stat_rdn);
                        %                         savefig("figures\PID\statisticsCB\PC" + string(pc) + "\REDUNDANCY_DIFFICULTY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_all_freqs_HZ.fig");
%                         save("results\PID\statisticsCB\PC" + string(pc) + "\stat_rdn_difficulty_" + fname_fig + "offset" + string(start_time * 1000) + "_all_freqs_2Hz.mat", "stat_rdn");

                        cfg_plot.zlim = [-0.0005, 0.0005];
                        figure
                        %cfg_plot.highlightchannel=[]
                        cfg_plot.parameter = 'stat_unq1';

                        title("UNIQUE1 DIFFICULTY " + pc + "-" + j)
                        %cfg_plot.highlightchannel = stat_unq1.label(h_unq1_np)
                        ft_clusterplot(cfg_plot, stat_unq1);
                        %                         savefig("figures\PID\statisticsCB\PC" + string(pc) + "\UNIQUE1_DIFFICULTY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_all_freqs_HZ.fig");
%                         save("results\PID\statisticsCB\PC" + string(pc) + "\stat_unq1_difficulty_" + fname_fig + "offset" + string(start_time * 1000) + "_all_freqs_2Hz.mat", "stat_unq1");

                        cfg_plot.zlim = [-0.0005, 0.000];
                        figure
                        %cfg_plot.highlightchannel=[]
                        %                         stat_unq2.raweffect = mean(Unq2_class1{i},1) - mean(Unq2_class2{i},1);
                        cfg_plot.parameter = 'stat_unq2';

                        title("UNIQUE2 DIFFICULTY " + pc + "-" + j)
                        %cfg_plot.highlightchannel = stat_unq2.label(h_unq2_np)
                        ft_clusterplot(cfg_plot, stat_unq2);
                        %                         savefig("figures\PID\statisticsCB\PC" + string(pc) + "\UNIQUE2_DIFFICULTY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_all_freqs_HZ.fig");
%                         save("results\PID\statisticsCB\PC" + string(pc) + "\stat_unq2_difficulty_" + fname_fig + "offset" + string(start_time * 1000) + "_all_freqs_2Hz.mat", "stat_unq2");

                        cfg_plot.zlim = [-0.0001, 0.0001];
                        figure
                        %cfg_plot.highlightchannel=[]
                        %                         stat_syn.raweffect = mean(Syn_class1{i},1) - mean(Syn_class2{i},1);
                        cfg_plot.parameter = 'stat_syn';

                        title("SYNERGY DIFFICULTY " + pc + "-" + j)
                        %cfg_plot.highlightchannel = stat_syn.label(h_syn_np)
                        ft_clusterplot(cfg_plot, stat_syn);
                        %                         savefig("figures\PID\statisticsCB\PC" + string(pc) + "\SYNERGY_DIFFICULTY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_all_freqs_HZ.fig");
%                         save("results\PID\statisticsCB\PC" + string(pc) + "\stat_syn_difficulty_" + fname_fig + "offset" + string(start_time * 1000) + "_all_freqs_2Hz.mat", "stat_syn");
                    %end
                end
                if(not(pc == PCs(end)))
                    cd 'results/PID/'
                end
            end
        end
    end
end