clear
clc
cd 'D:\Entrainment__\final_functions_behavior\'
addpath '..\partial-info-decomp-master'
addpath '..\fieldtrip-20210311'
addpath '..\gcmi-master\matlab'
addpath '..\analisi'

load 'channels_list';
chs_eeg=chs;
load('freq_template')
cmap = brewermap(1024,'RdBu');
cmap=flip(cmap);
n_pca = 4;

bad_chs=[22,28,32,41,46];

chs_eeg_ok=chs_eeg;
chs_eeg_ok(bad_chs)=[]; %non eeg
freq_fourier_all.label=cellstr(chs_eeg_ok);

cd 'results\behavior_results\PID\'

%% PID DIFFICULTY
for freq=0.5%linspace(0.5, 3.5, 7)
    fname =  "REDUNDANCY_offset1500ms_allFreqs_difficultySplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_std.mat";
    load (fname);
    fname = "REDUNDANCY_offset1500ms_allFreqs_difficultySplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_grandAverage.mat";
    load (fname);
    fname =  "UNIQUE1_offset1500ms_allFreqs_difficultySplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_std.mat";
    load (fname);
    fname = "UNIQUE1_offset1500ms_allFreqs_difficultySplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_grandAverage.mat";
    load (fname);
    fname =  "UNIQUE2_offset1500ms_allFreqs_difficultySplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_std.mat";
    load (fname);
    fname = "UNIQUE2_offset1500ms_allFreqs_difficultySplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_grandAverage.mat";
    load (fname);
    fname =  "SYNERGY_offset1500ms_allFreqs_difficultySplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_std.mat";
    load (fname);
    fname = "SYNERGY_offset1500ms_allFreqs_difficultySplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_grandAverage.mat";
    load (fname);
    
    figure
    for pc=1:n_pca
        lagged_redundancy_difficulty_grandAverage_difference = lagged_redundancy_subjects_difficulty_grandAverage{pc}(:, 1) - lagged_redundancy_subjects_difficulty_grandAverage{pc}(:, 2);
        lagged_unq1_difficulty_grandAverage_difference = lagged_unq1_subjects_difficulty_grandAverage{pc}(:, 1) - lagged_unq1_subjects_difficulty_grandAverage{pc}(:, 2);
        lagged_unq2_difficulty_grandAverage_difference = lagged_unq2_subjects_difficulty_grandAverage{pc}(:, 1) - lagged_unq2_subjects_difficulty_grandAverage{pc}(:, 2);
        lagged_synergy_difficulty_grandAverage_difference = lagged_synergy_subjects_difficulty_grandAverage{pc}(:, 1) - lagged_synergy_subjects_difficulty_grandAverage{pc}(:, 2);
        
        % Removing non eeg chs
        % grandAverage
        lagged_redundancy_subjects_difficulty_grandAverage{pc}(bad_chs, :) = [];
        lagged_unq1_subjects_difficulty_grandAverage{pc}(bad_chs, :) = [];
        lagged_unq2_subjects_difficulty_grandAverage{pc}(bad_chs, :) = [];
        lagged_synergy_subjects_difficulty_grandAverage{pc}(bad_chs, :) = [];
        
        % std
        lagged_redundancy_subjects_difficulty_std{pc}(bad_chs, :) = [];
        lagged_unq1_subjects_difficulty_std{pc}(bad_chs, :) = [];
        lagged_unq2_subjects_difficulty_std{pc}(bad_chs, :) = [];
        lagged_synergy_subjects_difficulty_std{pc}(bad_chs, :) = [];
        
        % difference
        lagged_redundancy_difficulty_grandAverage_difference(bad_chs, :) = [];
        lagged_unq1_difficulty_grandAverage_difference(bad_chs, :) = [];
        lagged_unq2_difficulty_grandAverage_difference(bad_chs, :) = [];
        lagged_synergy_difficulty_grandAverage_difference(bad_chs, :) = [];

        hlim_difficulty = max([max(lagged_redundancy_subjects_difficulty_grandAverage{pc}), max(lagged_unq1_subjects_difficulty_grandAverage{pc}), max(lagged_unq2_subjects_difficulty_grandAverage{pc}), max(lagged_synergy_subjects_difficulty_grandAverage{pc})]);
        hlim_difficulty_difference = max([abs(lagged_redundancy_difficulty_grandAverage_difference); abs(lagged_unq1_difficulty_grandAverage_difference); abs(lagged_unq2_difficulty_grandAverage_difference); abs(lagged_synergy_difficulty_grandAverage_difference)]);

        cfg=[];
        cfg.marker = 'off';%, 'numbers', 'off'
        colorb=true;
        cfg.colormap=cmap;
        freq_fourier_all.freq=[1];
        
        % Adjusting variables
        % grandAverage
        topo_data_red_easy_grandAverage = lagged_redundancy_subjects_difficulty_grandAverage{pc}(:, 1)';
        topo_data_red_hard_grandAverage = lagged_redundancy_subjects_difficulty_grandAverage{pc}(:, 2)';
        topo_data_unq1_easy_grandAverage = lagged_unq1_subjects_difficulty_grandAverage{pc}(:, 1)';
        topo_data_unq1_hard_grandAverage = lagged_unq1_subjects_difficulty_grandAverage{pc}(:, 2)';
        topo_data_unq2_easy_grandAverage = lagged_unq2_subjects_difficulty_grandAverage{pc}(:, 1)';
        topo_data_unq2_hard_grandAverage = lagged_unq2_subjects_difficulty_grandAverage{pc}(:, 2)';
        topo_data_syn_easy_grandAverage = lagged_synergy_subjects_difficulty_grandAverage{pc}(:, 1)';
        topo_data_syn_hard_grandAverage = lagged_synergy_subjects_difficulty_grandAverage{pc}(:, 2)';

        % std
        topo_data_red_easy_std = lagged_redundancy_subjects_difficulty_std{pc}(:, 1)';
        topo_data_red_hard_std = lagged_redundancy_subjects_difficulty_std{pc}(:, 2)';
        topo_data_unq1_easy_std = lagged_unq1_subjects_difficulty_std{pc}(:, 1)';
        topo_data_unq1_hard_std = lagged_unq1_subjects_difficulty_std{pc}(:, 2)';
        topo_data_unq2_easy_std = lagged_unq2_subjects_difficulty_std{pc}(:, 1)';
        topo_data_unq2_hard_std = lagged_unq2_subjects_difficulty_std{pc}(:, 2)';
        topo_data_syn_easy_std = lagged_synergy_subjects_difficulty_std{pc}(:, 1)';
        topo_data_syn_hard_std = lagged_synergy_subjects_difficulty_std{pc}(:, 2)';

        % difference
        topo_data_red_difficulty_grandAverage = lagged_redundancy_difficulty_grandAverage_difference';
        topo_data_unq1_difficulty_grandAverage = lagged_unq1_difficulty_grandAverage_difference';
        topo_data_unq2_difficulty_grandAverage = lagged_unq2_difficulty_grandAverage_difference';
        topo_data_syn_difficulty_grandAverage = lagged_synergy_difficulty_grandAverage_difference';
        
        % Plotting
        cfg.zlim=[0, hlim_difficulty];
        
        % grandAverage
%         freq_fourier_all.powspctrm=topo_data_red_easy_grandAverage;
%         title("REDUNDANCY easy " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\REDUNDANCY_DIFFICULTY_EASY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_red_hard_grandAverage;
%         title("REDUNDANCY hard " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\REDUNDANCY_DIFFICULTY_HARD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_unq1_easy_grandAverage;
%         title("UNIQUE1 easy " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE1_DIFFICULTY_EASY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_unq1_hard_grandAverage;
%         title("UNIQUE1 hard " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE1_DIFFICULTY_HARD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
% 
%         freq_fourier_all.powspctrm=topo_data_unq2_easy_grandAverage;
%         title("UNIQUE2 easy " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE2_DIFFICULTY_EASY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_unq2_hard_grandAverage;
%         title("UNIQUE2 hard " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE2_DIFFICULTY_HARD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%         
%         freq_fourier_all.powspctrm=topo_data_syn_easy_grandAverage;
%         title("SYNERGY easy " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\SYNERGY_DIFFICULTY_EASY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_syn_hard_grandAverage;
%         title("SYNERGY hard " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\SYNERGY_DIFFICULTY_HARD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%         
%         % std
%         freq_fourier_all.powspctrm=topo_data_red_easy_std;
%         title("REDUNDANCY easy " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\REDUNDANCY_DIFFICULTY_EASY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_red_hard_std;
%         title("REDUNDANCY hard " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\REDUNDANCY_DIFFICULTY_HARD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_unq1_easy_std;
%         title("UNIQUE1 easy " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE1_DIFFICULTY_EASY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_unq1_hard_std;
%         title("UNIQUE1 hard " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE1_DIFFICULTY_HARD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
% 
%         freq_fourier_all.powspctrm=topo_data_unq2_easy_std;
%         title("UNIQUE2 easy " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE2_DIFFICULTY_EASY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_unq2_hard_std;
%         title("UNIQUE2 hard " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE2_DIFFICULTY_HARD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%         
%         freq_fourier_all.powspctrm=topo_data_syn_easy_std;
%         title("SYNERGY easy " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\SYNERGY_DIFFICULTY_EASY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_syn_hard_std;
%         title("SYNERGY hard " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\SYNERGY_DIFFICULTY_HARD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
        
        % difference
        cfg.zlim = [-hlim_difficulty_difference, hlim_difficulty_difference];

        freq_fourier_all.powspctrm = topo_data_red_difficulty_grandAverage;
        title("REDUNDANCY DIFFICULTY " + string(freq) + "to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\PC" + string(pc) + "\REDUNDANCY_DIFFICULTY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")

        freq_fourier_all.powspctrm = topo_data_unq1_difficulty_grandAverage;
        title("UNIQUE1 DIFFICULTY " + string(freq) + "to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\PC" + string(pc) + "\UNIQUE1_DIFFICULTY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
        
        freq_fourier_all.powspctrm = topo_data_unq2_difficulty_grandAverage;
        title("UNIQUE2 DIFFICULTY " + string(freq) + "to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\PC" + string(pc) + "\UNIQUE2_DIFFICULTY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")

        freq_fourier_all.powspctrm = topo_data_syn_difficulty_grandAverage;
        title("SYNERGY DIFFICULTY " + string(freq) + "to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\PC" + string(pc) + "\SYNERGY_DIFFICULTY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
    end
end

%% PID PERFORMANCE
for freq=0.5%linspace(0.5, 3.5, 7)
    fname =  "REDUNDANCY_offset1500ms_allFreqs_performanceSplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_std.mat";
    load (fname);
    fname = "REDUNDANCY_offset1500ms_allFreqs_performanceSplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_grandAverage.mat";
    load (fname);
    fname =  "UNIQUE1_offset1500ms_allFreqs_performanceSplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_std.mat";
    load (fname);
    fname = "UNIQUE1_offset1500ms_allFreqs_performanceSplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_grandAverage.mat";
    load (fname);
    fname =  "UNIQUE2_offset1500ms_allFreqs_performanceSplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_std.mat";
    load (fname);
    fname = "UNIQUE2_offset1500ms_allFreqs_performanceSplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_grandAverage.mat";
    load (fname);
    fname =  "SYNERGY_offset1500ms_allFreqs_performanceSplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_std.mat";
    load (fname);
    fname = "SYNERGY_offset1500ms_allFreqs_performanceSplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_grandAverage.mat";
    load (fname);
    
    figure
    for pc=1:n_pca
        lagged_redundancy_performance_grandAverage_difference = lagged_redundancy_subjects_performance_grandAverage{pc}(:, 1) - lagged_redundancy_subjects_performance_grandAverage{pc}(:, 2);
        lagged_unq1_performance_grandAverage_difference = lagged_unq1_subjects_performance_grandAverage{pc}(:, 1) - lagged_unq1_subjects_performance_grandAverage{pc}(:, 2);
        lagged_unq2_performance_grandAverage_difference = lagged_unq2_subjects_performance_grandAverage{pc}(:, 1) - lagged_unq2_subjects_performance_grandAverage{pc}(:, 2);
        lagged_synergy_performance_grandAverage_difference = lagged_synergy_subjects_performance_grandAverage{pc}(:, 1) - lagged_synergy_subjects_performance_grandAverage{pc}(:, 2);
        
        % Removing non eeg chs
        % grandAverage
        lagged_redundancy_subjects_performance_grandAverage{pc}(bad_chs, :) = [];
        lagged_unq1_subjects_performance_grandAverage{pc}(bad_chs, :) = [];
        lagged_unq2_subjects_performance_grandAverage{pc}(bad_chs, :) = [];
        lagged_synergy_subjects_performance_grandAverage{pc}(bad_chs, :) = [];
        
        % std
        lagged_redundancy_subjects_performance_std{pc}(bad_chs, :) = [];
        lagged_unq1_subjects_performance_std{pc}(bad_chs, :) = [];
        lagged_unq2_subjects_performance_std{pc}(bad_chs, :) = [];
        lagged_synergy_subjects_performance_std{pc}(bad_chs, :) = [];
        
        % difference
        lagged_redundancy_performance_grandAverage_difference(bad_chs, :) = [];
        lagged_unq1_performance_grandAverage_difference(bad_chs, :) = [];
        lagged_unq2_performance_grandAverage_difference(bad_chs, :) = [];
        lagged_synergy_performance_grandAverage_difference(bad_chs, :) = [];

        hlim_performance = max([max(lagged_redundancy_subjects_performance_grandAverage{pc}), max(lagged_unq1_subjects_performance_grandAverage{pc}), max(lagged_unq2_subjects_performance_grandAverage{pc}), max(lagged_synergy_subjects_performance_grandAverage{pc})]);
        hlim_performance_difference = max([abs(lagged_redundancy_performance_grandAverage_difference); abs(lagged_unq1_performance_grandAverage_difference); abs(lagged_unq2_performance_grandAverage_difference); abs(lagged_synergy_performance_grandAverage_difference)]);

        cfg=[];
        cfg.marker = 'off';%, 'numbers', 'off'
        colorb=true;
        cfg.colormap=cmap;
        freq_fourier_all.freq=[1];
        
        % Adjusting variables
        % grandAverage
        topo_data_red_good_grandAverage = lagged_redundancy_subjects_performance_grandAverage{pc}(:, 1)';
        topo_data_red_bad_grandAverage = lagged_redundancy_subjects_performance_grandAverage{pc}(:, 2)';
        topo_data_unq1_good_grandAverage = lagged_unq1_subjects_performance_grandAverage{pc}(:, 1)';
        topo_data_unq1_bad_grandAverage = lagged_unq1_subjects_performance_grandAverage{pc}(:, 2)';
        topo_data_unq2_good_grandAverage = lagged_unq2_subjects_performance_grandAverage{pc}(:, 1)';
        topo_data_unq2_bad_grandAverage = lagged_unq2_subjects_performance_grandAverage{pc}(:, 2)';
        topo_data_syn_good_grandAverage = lagged_synergy_subjects_performance_grandAverage{pc}(:, 1)';
        topo_data_syn_bad_grandAverage = lagged_synergy_subjects_performance_grandAverage{pc}(:, 2)';

        % std
        topo_data_red_good_std = lagged_redundancy_subjects_performance_std{pc}(:, 1)';
        topo_data_red_bad_std = lagged_redundancy_subjects_performance_std{pc}(:, 2)';
        topo_data_unq1_good_std = lagged_unq1_subjects_performance_std{pc}(:, 1)';
        topo_data_unq1_bad_std = lagged_unq1_subjects_performance_std{pc}(:, 2)';
        topo_data_unq2_good_std = lagged_unq2_subjects_performance_std{pc}(:, 1)';
        topo_data_unq2_bad_std = lagged_unq2_subjects_performance_std{pc}(:, 2)';
        topo_data_syn_good_std = lagged_synergy_subjects_performance_std{pc}(:, 1)';
        topo_data_syn_bad_std = lagged_synergy_subjects_performance_std{pc}(:, 2)';

        % difference
        topo_data_red_performance_grandAverage = lagged_redundancy_performance_grandAverage_difference';
        topo_data_unq1_performance_grandAverage = lagged_unq1_performance_grandAverage_difference';
        topo_data_unq2_performance_grandAverage = lagged_unq2_performance_grandAverage_difference';
        topo_data_syn_performance_grandAverage = lagged_synergy_performance_grandAverage_difference';
        
        % Plotting
        cfg.zlim=[0, hlim_performance];
        
%         % grandAverage
%         freq_fourier_all.powspctrm=topo_data_red_good_grandAverage;
%         title("REDUNDANCY good " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\REDUNDANCY_PERFORMANCE_GOOD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_red_bad_grandAverage;
%         title("REDUNDANCY bad " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\REDUNDANCY_PERFORMANCE_BAD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_unq1_good_grandAverage;
%         title("UNIQUE1 good " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE1_PERFORMANCE_GOOD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_unq1_bad_grandAverage;
%         title("UNIQUE1 bad " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE1_PERFORMANCE_BAD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
% 
%         freq_fourier_all.powspctrm=topo_data_unq2_good_grandAverage;
%         title("UNIQUE2 good " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE2_PERFORMANCE_GOOD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_unq2_bad_grandAverage;
%         title("UNIQUE2 bad " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE2_PERFORMANCE_BAD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%         
%         freq_fourier_all.powspctrm=topo_data_syn_good_grandAverage;
%         title("SYNERGY good " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\SYNERGY_PERFORMANCE_GOOD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_syn_bad_grandAverage;
%         title("SYNERGY bad " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\SYNERGY_PERFORMANCE_BAD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
%         
%         % std
%         freq_fourier_all.powspctrm=topo_data_red_good_std;
%         title("REDUNDANCY good " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\REDUNDANCY_PERFORMANCE_GOOD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_red_bad_std;
%         title("REDUNDANCY bad " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\REDUNDANCY_PERFORMANCE_BAD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_unq1_good_std;
%         title("UNIQUE1 good " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE1_PERFORMANCE_GOOD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_unq1_bad_std;
%         title("UNIQUE1 bad " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE1_PERFORMANCE_BAD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
% 
%         freq_fourier_all.powspctrm=topo_data_unq2_good_std;
%         title("UNIQUE2 good " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE2_PERFORMANCE_GOOD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_unq2_bad_std;
%         title("UNIQUE2 bad " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\UNIQUE2_PERFORMANCE_BAD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%         
%         freq_fourier_all.powspctrm=topo_data_syn_good_std;
%         title("SYNERGY good " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\SYNERGY_PERFORMANCE_GOOD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%     
%         freq_fourier_all.powspctrm=topo_data_syn_bad_std;
%         title("SYNERGY bad " + string(freq) + " to " + string(freq + 9.5));
%         ft_topoplotER(cfg, freq_fourier_all);
%         savefig("figures\PC" + string(pc) + "\SYNERGY_PERFORMANCE_BAD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
%         
        % difference
        cfg.zlim = [-hlim_performance_difference, hlim_performance_difference];

        freq_fourier_all.powspctrm = topo_data_red_performance_grandAverage;
        title("REDUNDANCY PERFORMANCE " + string(freq) + "to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\PC" + string(pc) + "\REDUNDANCY_PERFORMANCE_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")

        freq_fourier_all.powspctrm = topo_data_unq1_performance_grandAverage;
        title("UNIQUE1 PERFORMANCE " + string(freq) + "to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\PC" + string(pc) + "\UNIQUE1_PERFORMANCE_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
        
        freq_fourier_all.powspctrm = topo_data_unq2_performance_grandAverage;
        title("UNIQUE2 PERFORMANCE " + string(freq) + "to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\PC" + string(pc) + "\UNIQUE2_PERFORMANCE_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")

        freq_fourier_all.powspctrm = topo_data_syn_performance_grandAverage;
        title("SYNERGY PERFORMANCE " + string(freq) + "to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\PC" + string(pc) + "\SYNERGY_PERFORMANCE_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
    end
end



