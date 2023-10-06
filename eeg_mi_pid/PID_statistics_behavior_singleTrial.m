clear
clc
close all

%% ADDING PATHS
addpath '../fieldtrip-20210311';
addpath './results';
addpath '../analisi';

%% PARAMETERS
cmap = brewermap(1024,'RdBu');
cmap = flip(cmap);
PCs = 1;
flows = 0.5;
fhighs = 10;
N_SUBJ = 22;
start_time = 0.5;
% How to split the single trials:
% 1. 'difficulty': easy vs hard phrases
% 2. 'answers': correct vs wrong answers
% 3. 'performance': good vs bad performance. Also specify the split number
% (from 1 to 5) which stands for different definitions of performance
% metric. 
split_method = 'difficulty'; %'answers'; %'performance';
% In case split_method is 'performance', the type of performance metric to
% use (valid from 1 to 5).
NUM_PERFORMANCE = 1;
% if true the data are saved otherwise not.
SAVE_DATA = false;
trigger_name = "startSound";
% Loading some templates
load 'chs_eeg_ok';
load('freq_template');
load('neighbours_ordered');

%% LOADING SUBS NAMES
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

%% STATISTICAL ANALYSIS
for freq=1:length(flows)
    flow = flows(freq);
    fhigh = fhighs(freq);
    %% FILE LOADING
    cd 'results/PID/';
    fname_redundancy = "REDUNDANCY_";
    fname_unq1 = "UNIQUE1_";
    fname_unq2 = "UNIQUE2_";
    fname_synergy = "SYNERGY_";
    fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_noSplit_";
    fname_bandp = "singleTrial_matchedLength_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";
    
    fname_rdn = fname_redundancy + fname_offset + fname_bandp;
    fname_unq1 = fname_unq1 + fname_offset + fname_bandp;
    fname_unq2 = fname_unq2 + fname_offset + fname_bandp;
    fname_syn = fname_synergy + fname_offset + fname_bandp;
    
    Rdn = load(fname_rdn).Rdn;
    Unq1 = load(fname_unq1).Unq1;
    Unq2 = load(fname_unq2).Unq2;
    Syn = load(fname_syn).Syn;
    
    %% STATISTICAL ANALYSIS FOR EACH PC
    for pc=PCs
        %% PID SPLIT DEPENDING ON THE METHOD SELECTED
        [PID_class1, PID_class2] = splitPID(split_method, subs, pc, Rdn, Unq1, Unq2, Syn);
        
        %% FORMATTING FOR CLUSTER BASED STATISTICS
        Rdn_class1 = zeros(N_SUBJ, 64);
        Rdn_class2 = zeros(N_SUBJ, 64);
        Unq1_class1 = zeros(N_SUBJ, 64);
        Unq1_class2 = zeros(N_SUBJ, 64);
        Unq2_class1 = zeros(N_SUBJ, 64);
        Unq2_class2 = zeros(N_SUBJ, 64);
        Syn_class1 = zeros(N_SUBJ, 64);
        Syn_class2 = zeros(N_SUBJ, 64);
        
        for sub=1:N_SUBJ
            Rdn_class1(sub, :) = mean(PID_class1.Rdn{sub}, 2)';
            Rdn_class2(sub, :) = mean(PID_class2.Rdn{sub}, 2)';
            Unq1_class1(sub, :) = mean(PID_class1.Unq1{sub}, 2)';
            Unq1_class2(sub, :) = mean(PID_class2.Unq1{sub}, 2)';
            Unq2_class1(sub, :) = mean(PID_class1.Unq2{sub}, 2)';
            Unq2_class2(sub, :) = mean(PID_class2.Unq2{sub}, 2)';
            Syn_class1(sub, :) = mean(PID_class1.Syn{sub}, 2)';
            Syn_class2(sub, :) = mean(PID_class2.Syn{sub}, 2)';
        end
        
        %% REMOVING NON EEG CHANNELS
        bad_chs=[22,28,32,41,46];       
        
        Rdn_class1(:, bad_chs) = [];
        Rdn_class2(:, bad_chs) = [];
        Unq1_class1(:, bad_chs) = [];
        Unq1_class2(:, bad_chs) = [];
        Unq2_class1(:, bad_chs) = [];
        Unq2_class2(:, bad_chs) = [];
        Syn_class1(:, bad_chs) = [];
        Syn_class2(:, bad_chs) = [];
        
        %% CLUSTER-BASED STATISTICS PARAMETER SETTINGS
        cfg=[];
        freq_fourier_all.label = cellstr(chs_eeg_ok);  
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
        
        allSubjRdnClass1.label  = chs_eeg_ok;
        allSubjRdnClass1.dimord = 'subj_chan_time';
        allSubjRdnClass1.time   = 1;
        allSubjRdnClass1.avg = Rdn_class1;
        allSubjRdnClass2.label  = chs_eeg_ok;
        allSubjRdnClass2.dimord = 'subj_chan_time';
        allSubjRdnClass2.time   = 1;
        allSubjRdnClass2.avg = Rdn_class2;
        
        
        allSubjUnq1Class1.label  = chs_eeg_ok;
        allSubjUnq1Class1.dimord = 'subj_chan_time';
        allSubjUnq1Class1.time   = 1;
        allSubjUnq1Class1.avg = Unq1_class1;
        allSubjUnq1Class2.label  = chs_eeg_ok;
        allSubjUnq1Class2.dimord = 'subj_chan_time';
        allSubjUnq1Class2.time   = 1;
        allSubjUnq1Class2.avg = Unq1_class2;
        
        
        allSubjUnq2Class1.label  = chs_eeg_ok;
        allSubjUnq2Class1.dimord = 'subj_chan_time';
        allSubjUnq2Class1.time   = 1;
        allSubjUnq2Class1.avg = Unq2_class1;
        allSubjUnq2Class2.label  = chs_eeg_ok;
        allSubjUnq2Class2.dimord = 'subj_chan_time';
        allSubjUnq2Class2.time   = 1;
        allSubjUnq2Class2.avg = Unq2_class2;
        
        allSubjSynClass1.label  = chs_eeg_ok;
        allSubjSynClass1.dimord = 'subj_chan_time';
        allSubjSynClass1.time   = 1;
        allSubjSynClass1.avg = Syn_class1;
        allSubjSynClass2.label  = chs_eeg_ok;
        allSubjSynClass2.dimord = 'subj_chan_time';
        allSubjSynClass2.time   = 1;
        allSubjSynClass2.avg = Syn_class2;
        
        %% CLUSTER-BASED STATISTICS COMPUTATION
        [stat_rdn] = ft_timelockstatistics(cfg, allSubjRdnClass1, allSubjRdnClass2)
        h_rdn_np=stat_rdn.mask;
        [stat_unq1] = ft_timelockstatistics(cfg, allSubjUnq1Class1, allSubjUnq1Class2)
        h_unq1_np=stat_unq1.mask;
        [stat_unq2] = ft_timelockstatistics(cfg, allSubjUnq2Class1, allSubjUnq2Class2)
        h_unq2_np=stat_unq2.mask;
        [stat_syn] = ft_timelockstatistics(cfg, allSubjSynClass1, allSubjSynClass2)
        h_syn_np=stat_syn.mask;
        
        %% PLOTTING RESULTS AS SCALP MAP AND DATA SAVING
        cfg_plot=[];
        cfg_plot.colormap=cmap;
        cfg_plot.marker =  'off'
        cfg_plot.highlight='on'
        cfg_plot.highlightchannel=[]
        cfg_plot.highlightsymbol='.';
        cfg_plot.highlightcolor='k';
        cfg_plot.highlightsize=18;
        fname_fig = "_";
        
        figure
        cfg_plot.zlim = [-0.001, 0.001];
        freq_fourier_all.powspctrm = mean(Rdn_class1) - mean(Rdn_class2);
        freq_fourier_all.freq=[1];
        title("REDUNDANCY")
        cfg_plot.highlightchannel = freq_fourier_all.label(h_rdn_np)
        ft_topoplotER(cfg_plot, freq_fourier_all);
        if(SAVE_DATA)
            savefig("../../figures/PID/singleTrial/statisticsCB/PC" + string(pc) + "/REDUNDANCY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
            save("../../results/PID/statisticsCB/PC" + string(pc) + "/stat_rdn_" + fname_fig + "offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_rdn");
        end
        
        figure
        cfg_plot.zlim = [-0.002, 0.002];
        cfg_plot.highlightchannel=[]
        freq_fourier_all.powspctrm = mean(Unq1_class1) - mean(Unq1_class2);
        title("UNIQUE SE")
        cfg_plot.highlightchannel = freq_fourier_all.label(h_unq1_np)
        ft_topoplotER(cfg_plot, freq_fourier_all);
        if(SAVE_DATA)
            savefig("../../figures/PID/singleTrial/statisticsCB/PC" + string(pc) + "/UNIQUE1_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
            save("../../results/PID/statisticsCB/PC" + string(pc) + "/stat_unq1_" + fname_fig + "offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_unq1");
        end
        
        figure
        cfg_plot.highlightchannel=[]
        freq_fourier_all.powspctrm = mean(Unq2_class1) - mean(Unq2_class2);
        title("UNIQUE PC")
        cfg_plot.highlightchannel = freq_fourier_all.label(h_unq2_np)
        ft_topoplotER(cfg_plot, freq_fourier_all);
        if(SAVE_DATA)
            savefig("../../figures/PID/singleTrial/statisticsCB/PC" + string(pc) + "/UNIQUE2_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
            save("../../results/PID/statisticsCB/PC" + string(pc) + "/stat_unq2_" + fname_fig + "offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_unq2");
        end
        
        figure
        cfg_plot.highlightchannel=[]
        cfg_plot.zlim = [-0.001, 0.001];
        freq_fourier_all.powspctrm = mean(Syn_class1) - mean(Syn_class2);
        title("SYNERGY")
        cfg_plot.highlightchannel = freq_fourier_all.label(h_syn_np)
        ft_topoplotER(cfg_plot, freq_fourier_all);
        if(SAVE_DATA)
            savefig("../../figures/PID/singleTrial/statisticsCB/PC" + string(pc) + "/SYNERGY_" + fname_fig + "OFFSET" + string(start_time * 1000) + "_" + flow + "TO" + fhigh + "HZ.fig");
            save("../../results/PID/statisticsCB/PC" + string(pc) + "/stat_syn" + fname_fig + "offset" + string(start_time * 1000) + "_" + flow + "to" + fhigh + "hz.mat", "stat_syn");
        end
    end
end

%% PID split. It taskes as input the input method and the PID atoms for each trials and returns well formatted PID trials splitted in two class depending on the selected splt_method
function [PID_class1, PID_class2] = splitPID(split_method, subs, pc, Rdn, Unq1, Unq2, Syn)
Rdn_class1_concat = {};
Rdn_class2_concat = {};
Unq1_class1_concat = {};
Unq1_class2_concat = {};
Unq2_class1_concat = {};
Unq2_class2_concat = {};
Syn_class1_concat = {};
Syn_class2_concat = {};
for sub=1:length(subs)
    for trial=1:length(Rdn{sub}{pc, 1})
        if(strcmp('answers', split_method))
            load("../../../behavioral_/results/subs/" + subs(sub) + "_behavior.mat");
            split_trials = behavior.answers;
            fname_fig = "answer_";
        elseif(strcmp('difficulty', split_method))
            load("../../../behavioral_/results/subs/" + subs(sub) + "_datasplit.mat");
            split_trials = dataSplit.easyTrials;
            fname_fig = "difficulty_";          
        elseif(strcmp('performance', split_method))
            load("../../../behavioral_/results/subs/" + subs(sub) + "_datasplit.mat");
            split_trials = dataSplit.goodPerformanceTrials{NUM_PERFORMANCE};
            fname_fig = "performance" + NUM_PERFORMANCE + "_";        
        end
        if(length(Rdn_class1_concat) < sub && split_trials(trial) == 1 && not(isempty(Rdn{sub}{pc, 1}{trial})))
            Rdn_class1_concat{sub} = Rdn{sub}{pc, 1}{1,trial};
            Unq1_class1_concat{sub} = Unq1{sub}{pc, 1}{1,trial};
            Unq2_class1_concat{sub}  = Unq2{sub}{pc, 1}{1,trial};
            Syn_class1_concat{sub} = Syn{sub}{pc, 1}{1,trial};
        elseif(length(Rdn_class2_concat) < sub && not(isempty(Rdn{sub}{pc, 1}{trial})))
            Rdn_class2_concat{sub} = Rdn{sub}{pc, 1}{1,trial};
            Unq1_class2_concat{sub} = Unq1{sub}{pc, 1}{1,trial};
            Unq2_class2_concat{sub}  = Unq2{sub}{pc, 1}{1,trial};
            Syn_class2_concat{sub} = Syn{sub}{pc, 1}{1,trial};
        else
            if(split_trials(trial) == 1 && not(isempty(Rdn{sub}{pc, 1}{trial})))
                Rdn_class1_concat{sub} = [Rdn_class1_concat{sub},Rdn{sub}{pc, 1}{1,trial}];
                Unq1_class1_concat{sub} = [Unq1_class1_concat{sub},Unq1{sub}{pc, 1}{1,trial}];
                Unq2_class1_concat{sub}  = [Unq2_class1_concat{sub},Unq2{sub}{pc, 1}{1,trial}];
                Syn_class1_concat{sub} = [Syn_class1_concat{sub}, Syn{sub}{pc, 1}{1,trial}];
            elseif(not(isempty(Rdn{sub}{pc, 1}{trial})))
                Rdn_class2_concat{sub} = [Rdn_class2_concat{sub}, Rdn{sub}{pc, 1}{1,trial}];
                Unq1_class2_concat{sub} = [Unq1_class2_concat{sub}, Unq1{sub}{pc, 1}{1,trial}];
                Unq2_class2_concat{sub}  = [Unq2_class2_concat{sub}, Unq2{sub}{pc, 1}{1,trial}];
                Syn_class2_concat{sub} = [Syn_class2_concat{sub}, Syn{sub}{pc, 1}{1,trial}];
            end
        end
    end
end
PID_class1.Rdn = Rdn_class1_concat;
PID_class1.Unq1 = Unq1_class1_concat;
PID_class1.Unq2 = Unq2_class1_concat;
PID_class1.Syn = Syn_class1_concat;
PID_class2.Rdn = Rdn_class2_concat;
PID_class2.Unq1 = Unq1_class2_concat;
PID_class2.Unq2 = Unq2_class2_concat;
PID_class2.Syn = Syn_class2_concat;
end