clear
close all
clc

addpath '../fieldtrip-20210311'
addpath './results'
addpath '../analisi'
cmap = brewermap(1024,'RdBu');
cmap = flip(cmap);
load 'chs_eeg_ok'
load('freq_template')
pc= 1;
addpath 'results/PID/'
Rdn_concat = load('REDUNDANCY_offset500ms_allFreqs_noSplit_trials_concat_singleTrial_matchedLength_eegBandp_0.5to10hz_filteredInputs.mat').Rdn;
Unq1_concat = load('UNIQUE1_offset500ms_allFreqs_noSplit_trials_concat_singleTrial_matchedLength_eegBandp_0.5to10hz_filteredInputs.mat').Unq1;
Unq2_concat = load('UNIQUE2_offset500ms_allFreqs_noSplit_trials_concat_singleTrial_matchedLength_eegBandp_0.5to10hz_filteredInputs.mat').Unq2;
Syn_concat= load('SYNERGY_offset500ms_allFreqs_noSplit_trials_concat_singleTrial_matchedLength_eegBandp_0.5to10hz_filteredInputs.mat').Syn;

Rdn_singleTrial = load('REDUNDANCY_offset500ms_allFreqs_noSplit_singleTrial_matchedLength_eegBandp_0.5to10hz_filteredInputs.mat').Rdn;
Unq1_singleTrial = load('UNIQUE1_offset500ms_allFreqs_noSplit_singleTrial_matchedLength_eegBandp_0.5to10hz_filteredInputs.mat').Unq1;
Unq2_singleTrial = load('UNIQUE2_offset500ms_allFreqs_noSplit_singleTrial_matchedLength_eegBandp_0.5to10hz_filteredInputs.mat').Unq2;
Syn_singleTrial = load('SYNERGY_offset500ms_allFreqs_noSplit_singleTrial_matchedLength_eegBandp_0.5to10hz_filteredInputs.mat').Syn;

for sub=9
    %for sample_size=1:7%length(Rdn{sub})
    %% REMOVING NON EEG CHANNELS
    bad_chs=[22,28,32,41,46];
    freq_fourier_all.label = cellstr(chs_eeg_ok);
    freq_fourier_all.freq=[1];
    cfg_plot=[];
    cfg_plot.highlightchannel=[];
    
    cfg_plot.colormap=cmap;
    
    % Difficulty
    Rdn_concat{1, sub}{1,1}(bad_chs) = [];
    Unq1_concat{1, sub}{1,1}(bad_chs) = [];
    Unq2_concat{1, sub}{1,1}(bad_chs) = [];
    Syn_concat{1, sub}{1,1}(bad_chs) = [];
    
    figure 
    freq_fourier_all.powspctrm = Rdn_concat{1, sub}{1,1};
    ft_topoplotER(cfg_plot, freq_fourier_all);
    
    figure
    freq_fourier_all.powspctrm = Unq1_concat{1, sub}{1,1};
    ft_topoplotER(cfg_plot, freq_fourier_all);
    
    figure
    freq_fourier_all.powspctrm = Unq2_concat{1, sub}{1,1};
    ft_topoplotER(cfg_plot, freq_fourier_all);
    
    figure
    freq_fourier_all.powspctrm = Syn_concat{1, sub}{1,1};
    ft_topoplotER(cfg_plot, freq_fourier_all);
    
    Rdn = {};
    Unq1 = {};
    Unq2 = {};
    Syn = {};
    for trial=1:length(Rdn_singleTrial{sub}{pc, 1})
        if(length(Rdn) < sub && not(isempty(Rdn_singleTrial{sub}{pc, 1}{trial})))
            Rdn{sub} = Rdn_singleTrial{sub}{pc, 1}{1,trial};
            Unq1{sub} = Unq1_singleTrial{sub}{pc, 1}{1,trial};
            Unq2{sub}  = Unq2_singleTrial{sub}{pc, 1}{1,trial};
            Syn{sub} = Syn_singleTrial{sub}{pc, 1}{1,trial};
        elseif(not(isempty(Rdn_singleTrial{sub}{pc, 1}{trial})))
            Rdn{sub} = [Rdn{sub},Rdn_singleTrial{sub}{pc, 1}{1,trial}];
            Unq1{sub} = [Unq1{sub},Unq1_singleTrial{sub}{pc, 1}{1,trial}];
            Unq2{sub}  = [Unq2{sub},Unq2_singleTrial{sub}{pc, 1}{1,trial}];
            Syn{sub} = [Syn{sub}, Syn_singleTrial{sub}{pc, 1}{1,trial}];
        end
    end
    
    Rdn{sub}(bad_chs, :) = [];
    Unq1{sub}(bad_chs, :) = [];
    Unq2{sub}(bad_chs, :) = [];
    Syn{sub}(bad_chs, :) = [];
    
    figure
    freq_fourier_all.powspctrm = mean(Rdn{sub}, 2);
    ft_topoplotER(cfg_plot, freq_fourier_all);
    
    figure
    freq_fourier_all.powspctrm = mean(Unq1{sub}, 2);
    ft_topoplotER(cfg_plot, freq_fourier_all);
    
    figure
    freq_fourier_all.powspctrm = mean(Unq2{sub}, 2);
    ft_topoplotER(cfg_plot, freq_fourier_all);
    
    figure
    freq_fourier_all.powspctrm = mean(Syn{sub}, 2);
    ft_topoplotER(cfg_plot, freq_fourier_all);
end

