clear
clc
close all
addpath(pwd)
cd 'results/PID'

load("REDUNDANCY_offset500ms_allFreqs_alltrials_checkbias_500reps_trials_concat_eegBandp_0.5to10hz_filteredInputs.mat")
load("UNIQUE1_offset500ms_allFreqs_alltrials_checkbias_500reps_trials_concat_eegBandp_0.5to10hz_filteredInputs.mat")
load("UNIQUE2_offset500ms_allFreqs_alltrials_checkbias_500reps_trials_concat_eegBandp_0.5to10hz_filteredInputs.mat")
load("SYNERGY_offset500ms_allFreqs_alltrials_checkbias_500reps_trials_concat_eegBandp_0.5to10hz_filteredInputs.mat")

topoplot=false;
cmap = brewermap(1024,'RdBu');
cmap = flip(cmap);
bad_chs=[22,28,32,41,46];
addpath '..\fieldtrip-20210311'
addpath '.\results'
addpath '..\analisi'
load 'chs_eeg_ok'
load('freq_template')

if(topoplot)
    cfg_plot=[];
    cfg_plot.colormap=cmap;
    cfg_plot.marker =  'off'
    freq_fourier_all.label = cellstr(chs_eeg_ok);

    for sub=1:1
        for len=1:3%length(Rdn{sub})
            if(not(isempty(Rdn{1, sub}{1, len})))
                Rdn_len = [];
                Syn_len = [];
                Unq1_len = [];
                Unq2_len = [];
                for rep=1:500
                    Rdn{1,sub}{1,len}{1,rep}(bad_chs) = [];
                    Unq1{1,sub}{1,len}{1,rep}(bad_chs) = [];
                    Unq2{1,sub}{1,len}{1,rep}(bad_chs) = [];
                    Syn{1,sub}{1,len}{1,rep}(bad_chs) = [];
                    if(rep ==1)
                        Rdn_len = [Rdn{1,sub}{1,len}{1,rep}];
                        Syn_len = [Syn{1,sub}{1,len}{1,rep}];
                        Unq1_len = [Unq1{1,sub}{1,len}{1,rep}];
                        Unq2_len = [Unq2{1,sub}{1,len}{1,rep}];
                    else
                        Rdn_len = [Rdn_len, Rdn{1,sub}{1,len}{1,rep}];
                        Syn_len = [Syn_len, Syn{1,sub}{1,len}{1,rep}];
                        Unq1_len = [Unq1_len, Unq1{1,sub}{1,len}{1,rep}];
                        Unq2_len = [Unq2_len, Unq2{1,sub}{1,len}{1,rep}];
                    end
                end
            end
            figure
            cfg_plot.zlim = [0, 0.1];
            freq_fourier_all.powspctrm = mean(Rdn_len, 2);
            freq_fourier_all.freq=[1];
            ft_topoplotER(cfg_plot, freq_fourier_all);
            cfg_plot.zlim = [0, 0.001];
            figure
            freq_fourier_all.powspctrm = std(Rdn_len,0, 2);
            freq_fourier_all.freq=[1];
            ft_topoplotER(cfg_plot, freq_fourier_all);
            cfg_plot.zlim = [0, 0.1];
            figure
            freq_fourier_all.powspctrm = mean(Unq1_len, 2);
            freq_fourier_all.freq=[1];
            ft_topoplotER(cfg_plot, freq_fourier_all);
            cfg_plot.zlim = [0, 0.001];
            figure
            freq_fourier_all.powspctrm = std(Unq1_len,0, 2);
            freq_fourier_all.freq=[1];
            ft_topoplotER(cfg_plot, freq_fourier_all);
            cfg_plot.zlim = [0, 0.1];
            figure
            freq_fourier_all.powspctrm = mean(Unq2_len, 2);
            freq_fourier_all.freq=[1];
            ft_topoplotER(cfg_plot, freq_fourier_all);
            cfg_plot.zlim = [0, 0.001];
            figure
            freq_fourier_all.powspctrm = std(Unq2_len,0, 2);
            freq_fourier_all.freq=[1];
            ft_topoplotER(cfg_plot, freq_fourier_all);
            cfg_plot.zlim = [0, 0.1];
            figure
            freq_fourier_all.powspctrm = mean(Syn_len, 2);
            freq_fourier_all.freq=[1];
            ft_topoplotER(cfg_plot, freq_fourier_all);
            cfg_plot.zlim = [0, 0.001];
            figure
            freq_fourier_all.powspctrm = std(Syn_len,0, 2);
            freq_fourier_all.freq=[1];
            ft_topoplotER(cfg_plot, freq_fourier_all);
        end
    end
end

mean_subs = [];
for sub=1:22
    mean_sub = [];
    std_sub = [];
    for len=1:length(Rdn{sub})
        if(not(isempty(Rdn{1, sub}{1, len})))
            Rdn_mean = [];
            Syn_mean = [];
            Unq1_mean = [];
            Unq2_mean = [];
            for rep=1:500
                Rdn{1,sub}{1,len}{1,rep}(bad_chs) = [];
                Unq1{1,sub}{1,len}{1,rep}(bad_chs) = [];
                Unq2{1,sub}{1,len}{1,rep}(bad_chs) = [];
                Syn{1,sub}{1,len}{1,rep}(bad_chs) = [];
                if(rep ==1)
                    Rdn_mean = [mean(Rdn{1,sub}{1,len}{1,rep})];
                    Syn_mean = [mean(Syn{1,sub}{1,len}{1,rep})];
                    Unq1_mean = [mean(Unq1{1,sub}{1,len}{1,rep})];
                    Unq2_mean = [mean(Unq2{1,sub}{1,len}{1,rep})];
                else
                    Rdn_mean = [Rdn_mean, mean(Rdn{1,sub}{1,len}{1,rep})];
                    Syn_mean = [Syn_mean, mean(Syn{1,sub}{1,len}{1,rep})];
                    Unq1_mean = [Unq1_mean, mean(Unq1{1,sub}{1,len}{1,rep})];
                    Unq2_mean = [Unq2_mean, mean(Unq2{1,sub}{1,len}{1,rep})];
                end
            end
        end
        if(len==1)
            mean_sub = [mean(Rdn_mean), mean(Unq1_mean), mean(Unq2_mean), mean(Syn_mean)];
            std_sub = [std(Rdn_mean), std(Unq1_mean), std(Unq2_mean), std(Syn_mean)];
        else
            mean_sub = [mean_sub; mean(Rdn_mean), mean(Unq1_mean), mean(Unq2_mean), mean(Syn_mean)];
            std_sub = [std_sub; std(Rdn_mean), std(Unq1_mean), std(Unq2_mean), std(Syn_mean)];
        end
    end
    figure
    errorbar(mean_sub(3:end, 1)/max(mean_sub(3:end, 1)), std_sub(3:end, 1)/max(mean_sub(3:end, 1)));
    hold on
    errorbar(mean_sub(3:end, 2)/max(mean_sub(3:end, 2)), std_sub(3:end, 2)/max(mean_sub(3:end, 2)));
    errorbar(mean_sub(3:end, 3)/max(mean_sub(3:end, 3)), std_sub(3:end, 3)/max(mean_sub(3:end, 3)));
    errorbar(mean_sub(3:end, 4)/max(mean_sub(3:end, 4)), std_sub(3:end, 4)/max(mean_sub(3:end, 4)));
    legend('Rdn', 'Unq1', 'Unq2', 'Syn');
    xticks([1, 2, 3, 4, 5]);
    xticklabels({'2^1^0', '2^1^2', '2^1^4', '2^1^6', '2^1^8'});
    ylim([0, 1]);
    saveas(gcf, "PID_reps_sub_" + sub + ".pdf");
%     saveas(gcf, "PID_reps_sub_" + sub + ".jpg");
%     saveas(gcf, "PID_reps_sub_" + sub + ".tif");

    if(sub == 1)
        mean_subs_rdn = [mean_sub(3:end, 1)/max(mean_sub(3:end, 1))];
        mean_subs_unq1 = [mean_sub(3:end, 2)/max(mean_sub(3:end, 2))];
        mean_subs_unq2 = [mean_sub(3:end, 3)/max(mean_sub(3:end, 3))];
        mean_subs_syn = [mean_sub(3:end, 4)/max(mean_sub(3:end, 4))];
    else
        if(sub ~= 22)
            mean_subs_rdn = [mean_subs_rdn, mean_sub(3:end, 1)/max(mean_sub(3:end, 1))];
            mean_subs_unq1 = [mean_subs_unq1, mean_sub(3:end, 2)/max(mean_sub(3:end, 2))];
            mean_subs_unq2 = [mean_subs_unq2, mean_sub(3:end, 3)/max(mean_sub(3:end, 3))];
            mean_subs_syn = [mean_subs_syn, mean_sub(3:end, 4)/max(mean_sub(3:end, 4))];
        else
            mean_subs_rdn = [mean_subs_rdn, [mean_sub(3:end, 1)/max(mean_sub(3:end, 1));nan]];
            mean_subs_unq1 = [mean_subs_unq1, [mean_sub(3:end, 2)/max(mean_sub(3:end, 2));nan]];
            mean_subs_unq2 = [mean_subs_unq2, [mean_sub(3:end, 3)/max(mean_sub(3:end, 3));nan]];
            mean_subs_syn = [mean_subs_syn, [mean_sub(3:end, 4)/max(mean_sub(3:end, 4));nan]];
        end
    end
end
figure
errorbar(mean(mean_subs_rdn, 2, 'omitnan'), [std(mean_subs_rdn(1:end-1, :), 0, 2, 'omitnan')/sqrt(22);std(mean_subs_rdn(end, :), 0, 2, 'omitnan')/sqrt(21)]);
hold on
errorbar(mean(mean_subs_unq1, 2, 'omitnan'), [std(mean_subs_unq1(1:end-1, :), 0, 2, 'omitnan')/sqrt(22);std(mean_subs_unq1(end, :), 0, 2, 'omitnan')/sqrt(21)]);
errorbar(mean(mean_subs_unq2, 2, 'omitnan'), [std(mean_subs_unq2(1:end-1, :), 0, 2, 'omitnan')/sqrt(22);std(mean_subs_unq2(end, :), 0, 2, 'omitnan')/sqrt(21)]);
errorbar(mean(mean_subs_syn, 2, 'omitnan'), [std(mean_subs_syn(1:end-1, :), 0, 2, 'omitnan')/sqrt(22);std(mean_subs_syn(end, :), 0, 2, 'omitnan')/sqrt(21)]);
legend('Rdn', 'Unq1', 'Unq2', 'Syn');
xticks([1,2, 3,4, 5]);
xticklabels({'2^1^0', '2^1^2', '2^1^4', '2^1^6', '2^1^8'});
ylim([0, 1]);