clear
clc
cd 'D:\Entrainment__\final_functions\'
addpath '..\partial-info-decomp-master'
addpath '..\fieldtrip-20210311'
addpath '..\gcmi-master\matlab'
addpath '..\analisi'

load 'channels_list';
chs_eeg=chs;
load('freq_template')
cmap = brewermap(1024,'RdBu');
cmap=flip(cmap);
n_pca = 1;

bad_chs=[22,28,32,41,46];

chs_eeg_ok=chs_eeg;
chs_eeg_ok(bad_chs)=[]; %non eeg
freq_fourier_all.label=cellstr(chs_eeg_ok);

cd 'results\behavior_results\MI\'
for freq=0.5
    fname =  "MI_ENVELOPE_offset1500ms_allFreqs_difficultySplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_std.mat";
    load (fname);
    fname = "MI_ENVELOPE_offset1500ms_allFreqs_performanceSplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_std.mat";
    load (fname);
    fname =  "MI_ENVELOPE_offset1500ms_allFreqs_difficultySplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_grandAverage.mat";
    load (fname);
    fname = "MI_ENVELOPE_offset1500ms_allFreqs_performanceSplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_grandAverage.mat";
    load (fname);
    lagged_MI_env_difficulty_grandAverage_difference = lagged_MI_env_difficulty_grandAverage(:, 1) - lagged_MI_env_difficulty_grandAverage(:, 2); 
    lagged_MI_env_performance_grandAverage_difference = lagged_MI_env_performance_grandAverage(:, 1) - lagged_MI_env_performance_grandAverage(:, 2);

    figure
    lagged_MI_env_difficulty_std(bad_chs, :) = [];
    lagged_MI_env_performance_std(bad_chs, :) = [];
    lagged_MI_env_difficulty_grandAverage(bad_chs, :) = [];
    lagged_MI_env_performance_grandAverage(bad_chs, :) = [];
    lagged_MI_env_difficulty_grandAverage_difference(bad_chs, :) = [];
    lagged_MI_env_performance_grandAverage_difference(bad_chs, :) = [];


    llim_difficulty = min(min(lagged_MI_env_difficulty_grandAverage));
    hlim_difficulty = max(max(lagged_MI_env_difficulty_grandAverage)); 
    llim_performance = min(min(lagged_MI_env_performance_grandAverage));
    hlim_performance = max(max(lagged_MI_env_performance_grandAverage));
    hlim_difficulty_difference = max(abs(lagged_MI_env_difficulty_grandAverage_difference));
    hlim_performance_difference = max(abs(lagged_MI_env_performance_grandAverage_difference));

    
    cfg=[];
    cfg.marker = 'off';%, 'numbers', 'off'
    colorb=true;
    cfg.colormap=cmap;
    freq_fourier_all.freq=[1];
    
    topo_data_easy_std = lagged_MI_env_difficulty_std(:, 1)';
    topo_data_hard_std = lagged_MI_env_difficulty_std(:, 2)';
    topo_data_good_std = lagged_MI_env_performance_std(:, 1)';
    topo_data_bad_std = lagged_MI_env_performance_std(:, 2)';
    
    topo_data_easy_grandAverage = lagged_MI_env_difficulty_grandAverage(:, 1)';
    topo_data_hard_grandAverage = lagged_MI_env_difficulty_grandAverage(:, 2)';
    topo_data_good_grandAverage = lagged_MI_env_performance_grandAverage(:, 1)';
    topo_data_bad_grandAverage = lagged_MI_env_performance_grandAverage(:, 2)';
    topo_data_difficulty_grandAverage = lagged_MI_env_difficulty_grandAverage_difference';
    topo_data_performance_grandAverage = lagged_MI_env_performance_grandAverage_difference';
    
    cfg.zlim = [llim_difficulty, hlim_difficulty];
    freq_fourier_all.powspctrm=topo_data_easy_std;
    title("ENVELOPE easy " + string(freq) + "to " + string(freq + 9.5));
    ft_topoplotER(cfg, freq_fourier_all);
    savefig("figures\MI_ENVELOPE_DIFFICULTY_EASY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")

    freq_fourier_all.powspctrm=topo_data_easy_grandAverage;
    title("ENVELOPE easy " + string(freq) + "to " + string(freq + 9.5));
    ft_topoplotER(cfg, freq_fourier_all);
    savefig("figures\MI_ENVELOPE_DIFFICULTY_EASY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")

    freq_fourier_all.powspctrm=topo_data_hard_std;
    title("ENVELOPE hard " + string(freq) + "to " + string(freq + 9.5));
    ft_topoplotER(cfg, freq_fourier_all);
    savefig("figures\MI_ENVELOPE_DIFFICULTY_HARD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
    
    freq_fourier_all.powspctrm=topo_data_hard_grandAverage;
    title("ENVELOPE hard " + string(freq) + "to " + string(freq + 9.5));
    ft_topoplotER(cfg, freq_fourier_all);
    savefig("figures\MI_ENVELOPE_DIFFICULTY_HARD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")

    cfg.zlim = [-hlim_difficulty_difference, hlim_difficulty_difference];
    freq_fourier_all.powspctrm = topo_data_difficulty_grandAverage;
    title("ENVELOPE hard " + string(freq) + "to " + string(freq + 9.5));
    ft_topoplotER(cfg, freq_fourier_all);
    savefig("figures\MI_ENVELOPE_DIFFICULTY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")

    cfg.zlim = [llim_performance, hlim_performance];
    freq_fourier_all.powspctrm=topo_data_good_std;
    title("ENVELOPE good " + string(freq) + "to " + string(freq + 9.5));
    ft_topoplotER(cfg, freq_fourier_all);
    savefig("figures\MI_ENVELOPE_PERFORMANCE_GOOD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")

    freq_fourier_all.powspctrm=topo_data_good_grandAverage;
    title("ENVELOPE good " + string(freq) + "to " + string(freq + 9.5));
    ft_topoplotER(cfg, freq_fourier_all);
    savefig("figures\MI_ENVELOPE_PERFORMANCE_GOOD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")

    freq_fourier_all.powspctrm=topo_data_bad_std;
    title("ENVELOPE bad " + string(freq) + "to " + string(freq + 9.5));
    ft_topoplotER(cfg, freq_fourier_all);
    savefig("figures\MI_ENVELOPE_PERFORMANCE_BAD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
    
    freq_fourier_all.powspctrm=topo_data_bad_grandAverage;
    title("ENVELOPE bad " + string(freq) + "to " + string(freq + 9.5));
    ft_topoplotER(cfg, freq_fourier_all);
    savefig("figures\MI_ENVELOPE_PERFORMANCE_BAD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")

    cfg.zlim = [-hlim_performance_difference, hlim_performance_difference];
    freq_fourier_all.powspctrm = topo_data_performance_grandAverage;
    title("ENVELOPE bad " + string(freq) + "to " + string(freq + 9.5));
    ft_topoplotER(cfg, freq_fourier_all);
    savefig("figures\MI_ENVELOPE_PERFORMANCE_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")

end

for freq=0.5%linspace(0.5, 3.5, 7)
    fname =  "MI_PCs_offset1500ms_allFreqs_difficultySplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_std.mat";
    load (fname);
    fname = "MI_PCs_offset1500ms_allFreqs_performanceSplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_std.mat";
    load (fname);
    fname =  "MI_PCs_offset1500ms_allFreqs_difficultySplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_grandAverage.mat";
    load (fname);
    fname = "MI_PCs_offset1500ms_allFreqs_performanceSplit_trials_concat_eegBandp_" + string(freq) + "to" + string(freq + 9.5) + "hz_filteredinputs_grandAverage.mat";
    load (fname);
    
    figure
    for pc=1:n_pca
        lagged_MI_pc_difficulty_grandAverage_difference = lagged_MI_pc_difficulty_grandAverage{pc}(:, 1) - lagged_MI_pc_difficulty_grandAverage{pc}(:, 2); 
        lagged_MI_pc_performance_grandAverage_difference = lagged_MI_pc_performance_grandAverage{pc}(:, 1) - lagged_MI_pc_performance_grandAverage{pc}(:, 2);
    
        lagged_MI_pc_difficulty_std{pc}(bad_chs, :) = [];
        lagged_MI_pc_performance_std{pc}(bad_chs, :) = [];
        lagged_MI_pc_difficulty_grandAverage{pc}(bad_chs, :) = [];
        lagged_MI_pc_performance_grandAverage{pc}(bad_chs, :) = [];
        lagged_MI_pc_difficulty_grandAverage_difference(bad_chs, :) = [];
        lagged_MI_pc_performance_grandAverage_difference(bad_chs, :) = [];

        llim_difficulty = min(min(lagged_MI_pc_difficulty_grandAverage{pc})); 
        hlim_difficulty = max(max(lagged_MI_pc_difficulty_grandAverage{pc}));
        llim_performance = min(min(lagged_MI_pc_performance_grandAverage{pc}));
        hlim_performance = max(max(lagged_MI_pc_performance_grandAverage{pc}));
        hlim_difficulty_difference = max(abs(lagged_MI_pc_difficulty_grandAverage_difference));
        hlim_performance_difference = max(abs(lagged_MI_pc_performance_grandAverage_difference));

        cfg=[];
        cfg.marker = 'off';%, 'numbers', 'off'
        colorb=true;
        cfg.colormap=cmap;
        freq_fourier_all.freq=[1];
        
        topo_data_easy_std = lagged_MI_pc_difficulty_std{pc}(:, 1)';
        topo_data_hard_std = lagged_MI_pc_difficulty_std{pc}(:, 2)';
        topo_data_good_std = lagged_MI_pc_performance_std{pc}(:, 1)';
        topo_data_bad_std = lagged_MI_pc_performance_std{pc}(:, 2)';

        topo_data_easy_grandAverage = lagged_MI_pc_difficulty_grandAverage{pc}(:, 1)';
        topo_data_hard_grandAverage = lagged_MI_pc_difficulty_grandAverage{pc}(:, 2)';
        topo_data_good_grandAverage = lagged_MI_pc_performance_grandAverage{pc}(:, 1)';
        topo_data_bad_grandAverage = lagged_MI_pc_performance_grandAverage{pc}(:, 2)';
        topo_data_difficulty_grandAverage = lagged_MI_pc_difficulty_grandAverage_difference';
        topo_data_performance_grandAverage = lagged_MI_pc_performance_grandAverage_difference';
        
        cfg.zlim=[llim_difficulty, hlim_difficulty];
        freq_fourier_all.powspctrm=topo_data_easy_std;
        title("PC" + pc +  "easy " + string(freq) + " to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\MI_PC" + pc + "_DIFFICULTY_EASY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
        
        freq_fourier_all.powspctrm=topo_data_easy_grandAverage;
        title("PC" + pc +  "easy " + string(freq) + " to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\MI_PC" + pc + "_DIFFICULTY_EASY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
    
        freq_fourier_all.powspctrm=topo_data_hard_std;
        title("PC" + pc +  "hard " + string(freq) + " to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\MI_PC" + pc + "_DIFFICULTY_HARD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")

        freq_fourier_all.powspctrm=topo_data_hard_grandAverage;
        title("PC" + pc +  "hard " + string(freq) + " to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\MI_PC" + pc + "_DIFFICULTY_HARD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
    
        cfg.zlim = [-hlim_difficulty_difference, hlim_difficulty_difference];
        freq_fourier_all.powspctrm = topo_data_difficulty_grandAverage;
        title("PC" + pc + " hard " + string(freq) + "to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\MI_PC" + pc + "_DIFFICULTY_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")

        cfg.zlim=[llim_performance, hlim_performance];
        freq_fourier_all.powspctrm=topo_data_good_std;
        title("PC" + pc +  "good " + string(freq) + " to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\MI_PC" + pc + "_PERFORMANCE_GOOD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")
    
        freq_fourier_all.powspctrm=topo_data_good_grandAverage;
        title("PC" + pc +  "good " + string(freq) + " to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\MI_PC" + pc + "_PERFORMANCE_GOOD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")

        freq_fourier_all.powspctrm=topo_data_bad_std;
        title("PC" + pc +  "bad " + string(freq) + " to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\MI_PC" + pc + "_PERFORMANCE_BAD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ_STD.fig")

        freq_fourier_all.powspctrm=topo_data_bad_grandAverage;
        title("PC" + pc +  "bad " + string(freq) + " to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\MI_PC" + pc + "_PERFORMANCE_BAD_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")
        
        cfg.zlim = [-hlim_performance_difference, hlim_performance_difference];
        freq_fourier_all.powspctrm = topo_data_performance_grandAverage;
        title("PC" + pc + " bad " + string(freq) + "to " + string(freq + 9.5));
        ft_topoplotER(cfg, freq_fourier_all);
        savefig("figures\MI_PC" + pc + "_PERFORMANCE_" + string(freq) + "TO" + string(freq+ 9.5) + "HZ.fig")    
    end
end

