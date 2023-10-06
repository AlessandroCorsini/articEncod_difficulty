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

flows=[5];
fhighs=[7];
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


    Rdn_mean_all_subj_class1 = cell(1,N_SUBJ);
    Unq1_mean_all_subj_class1 = cell(1,N_SUBJ);
    Unq2_mean_all_subj_class1 = cell(1,N_SUBJ);
    Syn_mean_all_subj_class1 = cell(1,N_SUBJ);

    Rdn_mean_all_subj_class2 = cell(1,N_SUBJ);
    Unq1_mean_all_subj_class2 = cell(1,N_SUBJ);
    Unq2_mean_all_subj_class2 = cell(1,N_SUBJ);
    Syn_mean_all_subj_class2 = cell(1,N_SUBJ);

    for subj=1:length(Rdn{1,1})

        Rdn_mean = zeros(64,2); % easy and hard
        Unq1_mean = zeros(64,2);
        Unq2_mean = zeros(64,2);
        Syn_mean = zeros(64,2);

        for iter=1:length(Rdn{1,1}{1,subj})

            Rdn_mean = Rdn_mean + Rdn{1,1}{1,subj}{1,iter}{1,1};
            Unq1_mean = Unq1_mean + Unq1{1,1}{1,subj}{1,iter}{1,1};
            Unq2_mean = Unq2_mean + Unq2{1,1}{1,subj}{1,iter}{1,1};
            Syn_mean = Syn_mean + Syn{1,1}{1,subj}{1,iter}{1,1};

        end

        Rdn_mean = Rdn_mean/length(Rdn{1,1}{1,subj});
        Unq1_mean = Unq1_mean/length(Unq1{1,1}{1,subj});
        Unq2_mean = Unq2_mean/length(Unq2{1,1}{1,subj});
        Syn_mean = Syn_mean/length(Syn{1,1}{1,subj});

        Rdn_mean_all_subj_class1{1,subj} = Rdn_mean(:,1);
        Unq1_mean_all_subj_class1{1,subj} = Unq1_mean(:,1);
        Unq2_mean_all_subj_class1{1,subj} = Unq2_mean(:,1);
        Syn_mean_all_subj_class1{1,subj} = Syn_mean(:,1);
        Rdn_mean_all_subj_class2{1,subj} = Rdn_mean(:,2);
        Unq1_mean_all_subj_class2{1,subj} = Unq1_mean(:,2);
        Unq2_mean_all_subj_class2{1,subj} = Unq2_mean(:,2);
        Syn_mean_all_subj_class2{1,subj} = Syn_mean(:,2);
    end


    % Difficulty
    Rdn_class1 = zeros(N_SUBJ, 64);
    Rdn_class2 = zeros(N_SUBJ, 64);
    Unq1_class1 = zeros(N_SUBJ, 64);
    Unq1_class2 = zeros(N_SUBJ, 64);
    Unq2_class1 = zeros(N_SUBJ, 64);
    Unq2_class2 = zeros(N_SUBJ, 64);
    Syn_class1 = zeros(N_SUBJ, 64);
    Syn_class2 = zeros(N_SUBJ, 64);
    for sub=1:N_SUBJ
        Rdn_class1(sub, :) = Rdn_mean_all_subj_class1{1,sub}';
        Rdn_class2(sub, :) = Rdn_mean_all_subj_class2{1,sub}';
        Unq1_class1(sub, :) = Unq1_mean_all_subj_class1{1,sub}';
        Unq1_class2(sub, :) = Unq1_mean_all_subj_class2{1,sub}';
        Unq2_class1(sub, :) = Unq2_mean_all_subj_class1{1,sub}';
        Unq2_class2(sub, :) = Unq2_mean_all_subj_class2{1,sub}';
        Syn_class1(sub, :) = Syn_mean_all_subj_class1{1,sub}';
        Syn_class2(sub, :) = Syn_mean_all_subj_class2{1,sub}';
    end

    %% REMOVING NON EEG CHANNELS
    % Difficulty
    Rdn_class1(:, bad_chs) = [];
    Rdn_class2(:, bad_chs) = [];
    Unq1_class1(:, bad_chs) = [];
    Unq1_class2(:, bad_chs) = [];
    Unq2_class1(:, bad_chs) = [];
    Unq2_class2(:, bad_chs) = [];
    Syn_class1(:, bad_chs) = [];
    Syn_class2(:, bad_chs) = [];
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
        allSubjRdnClass1.avg = Rdn_class1;
        allSubjRdnClass2.label  = chs_eeg_ok;
        allSubjRdnClass2.dimord = 'subj_chan_time';
        allSubjRdnClass2.time   = 1;
        allSubjRdnClass2.avg = Rdn_class2;

        allSubjUnq1Class1 = [];
        allSubjUnq1Class2 = [];
        allSubjUnq1Class1.label  = chs_eeg_ok;
        allSubjUnq1Class1.dimord = 'subj_chan_time';
        allSubjUnq1Class1.time   = 1;
        allSubjUnq1Class1.avg = Unq1_class1;
        allSubjUnq1Class2.label  = chs_eeg_ok;
        allSubjUnq1Class2.dimord = 'subj_chan_time';
        allSubjUnq1Class2.time   = 1;
        allSubjUnq1Class2.avg = Unq1_class2;

        allSubjUnq2Class1 = [];
        allSubjUnq2Class2 = [];
        allSubjUnq2Class1.label  = chs_eeg_ok;
        allSubjUnq2Class1.dimord = 'subj_chan_time';
        allSubjUnq2Class1.time   = 1;
        allSubjUnq2Class1.avg = Unq2_class1;
        allSubjUnq2Class2.label  = chs_eeg_ok;
        allSubjUnq2Class2.dimord = 'subj_chan_time';
        allSubjUnq2Class2.time   = 1;
        allSubjUnq2Class2.avg = Unq2_class2;

        allSubjSynClass1 = [];
        allSubjSynClass2 = [];
        allSubjSynClass1.label  = chs_eeg_ok;
        allSubjSynClass1.dimord = 'subj_chan_time';
        allSubjSynClass1.time   = 1;
        allSubjSynClass1.avg = Syn_class1;
        allSubjSynClass2.label  = chs_eeg_ok;
        allSubjSynClass2.dimord = 'subj_chan_time';
        allSubjSynClass2.time   = 1;
        allSubjSynClass2.avg = Syn_class2;

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
        cfg_plot.zlim = [-0.0002, 0.0002];
        freq_fourier_all.powspctrm = mean(Rdn_class1) - mean(Rdn_class2);
        %                             ft_topoplotER(cfg_plot, freq_fourier_all);
        %                             title("REDUNDANCY DIFFICULTY " + pc + "-" + j)
        cfg_plot.highlightchannel = freq_fourier_all.label(h_rdn_np)
        ft_topoplotER(cfg_plot, freq_fourier_all);
        %                            savefig("../../figures/PID/statisticsCB/PC" + string(pc) + "/REDUNDANCY_DIFFICULTY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
        %                            save("../../results/PID/statisticsCB/PC" + string(pc) + "/stat_rdn_difficulty_" + fname_fig + "offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_rdn");

        cfg_plot.zlim = [-0.0003, 0.0003];
        figure
        cfg_plot.highlightchannel=[]
        freq_fourier_all.powspctrm = mean(Unq1_class1) - mean(Unq1_class2);
        ft_topoplotER(cfg_plot, freq_fourier_all);
        title("UNIQUE1 SE")
        cfg_plot.highlightchannel = freq_fourier_all.label(h_unq1_np)
        ft_topoplotER(cfg_plot, freq_fourier_all);
        % savefig("../../figures/PID/statisticsCB/PC" + string(pc) + "/UNIQUE1_DIFFICULTY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
        % save("../../results/PID/statisticsCB/PC" + string(pc) + "/stat_unq1_difficulty_" + fname_fig + "offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_unq1");

        cfg_plot.zlim = [-0.0003, 0.0003];
        figure
        cfg_plot.highlightchannel=[]
        freq_fourier_all.powspctrm = mean(Unq2_class1) - mean(Unq2_class2);
        ft_topoplotER(cfg_plot, freq_fourier_all);
        title("UNIQUE PC1")
        cfg_plot.highlightchannel = freq_fourier_all.label(h_unq2_np)
        ft_topoplotER(cfg_plot, freq_fourier_all);
        % savefig("../../figures/PID/statisticsCB/PC" + string(pc) + "/UNIQUE2_DIFFICULTY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
        % save("../../results/PID/statisticsCB/PC" + string(pc) + "/stat_unq2_difficulty_" + fname_fig + "offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_unq2");

        cfg_plot.zlim = [-0.0001, 0.0001];
        figure
        cfg_plot.highlightchannel=[]
        freq_fourier_all.powspctrm = mean(Syn_class1) - mean(Syn_class2);
        %                             ft_topoplotER(cfg_plot, freq_fourier_all);
        cfg_plot.highlightchannel = freq_fourier_all.label(h_syn_np)
        ft_topoplotER(cfg_plot,freq_fourier_all);
        %                             savefig("../../figures/PID/statisticsCB/PC" + string(pc) + "/SYNERGY_DIFFICULTY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
        %                             save("../../results/PID/statisticsCB/PC" + string(pc) + "/stat_syn_difficulty_" + fname_fig + "offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_syn");
    end
end
