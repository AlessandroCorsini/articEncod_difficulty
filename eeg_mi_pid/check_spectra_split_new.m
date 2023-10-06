clear
clc
close all

%% ADDING PATHS AND SETTING PARAMS
addpath '../fieldtrip-20210311'
addpath './results'
addpath '../analisi'
addpath(pwd)
Fs=400;
INVALID_PHRASE = 34;
pca_dir="../pca_new_nocut";
env_dir="../envelope";
load('cut_times')
cut_times=round(cut_times*Fs*0.001);
PC = 1;

%% LOADING PCA NAMES AND DATA
d_pca_em=dir(pca_dir);
pca_names={};
for i=4:length(d_pca_em)
    pca_names=[pca_names; d_pca_em(i).name];
end
pca_names=natsort(pca_names);
pcas={};
for i=1:length(pca_names)
    pca_tmp=load(pca_dir+"/"+pca_names{i});
    pcas=[pcas pca_tmp.pca_tmp];
end

%% LOADING ENVELOPE NAMES AND DATA
d_env=dir(env_dir);
envelope_names={};
for i=3:length(d_env)
    envelope_names=[envelope_names; d_env(i).name];
end
envelope_names=natsort(envelope_names);
envelopes={};
for i=1:length(envelope_names)
    envelopes=[envelopes readNPY(env_dir+"/"+envelope_names{i})];
end

%% MATCH THE LENGTH OF ENVELOPES WITH PCAS
for i=1:length(envelopes)
    diff = length(envelopes{i})-length(pcas{i});
    envelopes{i}=envelopes{i}(1:end-diff);
end

%% CUTTING TIMES
for i=1:length(cut_times)
    start_cutTime=cut_times(i,1);
    end_cutTime=cut_times(i,2);
    pcas{i}=pcas{i}(:,start_cutTime:end-end_cutTime);
    envelopes{1,i}=envelopes{1,i}(start_cutTime:end-end_cutTime);
end

%% Exclude invalid pcas and envelopes
corrupted_pca=[];
for i=1:length(pcas)
    if length(pcas{i}) < 2000
        corrupted_pca = [corrupted_pca i];
    end
end
pcas{corrupted_pca} = NaN;
pcas{INVALID_PHRASE} = NaN;
envelopes{corrupted_pca} = NaN;
envelopes{INVALID_PHRASE} = NaN;

%% Compute max length for padding
L_max = 0;
L_min = length(pcas{1});
for i=1:length(pcas)
    if(length(pcas{i}) > L_max)
        L_max = length(pcas{i});
    end
    if(not(isnan(pcas{i})) & length(pcas{i}) < L_min)
        L_min = length(pcas{i});
    end
end

%% Split the pcas in easy and hard
load("../behavioral_/results/phrases/behaviorPhrases.mat");
[accuracy_ordered, idxs_phrase] = sort(behaviorPhrases.accuracy);
pcas_easy = pcas(idxs_phrase(1:25));
pcas_hard = pcas(idxs_phrase(26:end));
envs_easy = envelopes(idxs_phrase(1:25));
envs_hard = envelopes(idxs_phrase(26:end));
pc=1;
while pc<= size(pcas_easy, 2)
    if(isnan(pcas_easy{1,pc}))
        pcas_easy(pc) = [];
        envs_easy(pc) = [];
    end
    pc = pc+1;
end
pc=1;
while pc<= size(pcas_hard, 2)
    if(isnan(pcas_hard{1,pc}))
        pcas_hard(pc) = [];
        envs_hard(pc) = [];
    end
    pc=pc+1;
end

ts = 1/Fs;
time_easy = {};
time_hard = {};
for stim=1:length(pcas_easy)
    data_easy.time{stim} = [ts:ts:ts*length(pcas_easy{stim})];
    tmp_pc = pcas_easy{stim};
    tmp_env = envs_easy{stim};
    data_easy.trial{stim} = [tmp_env';tmp_pc(1:4,:)];   
end

%% EASY
data_easy.hdr.chantype = ["misc";"misc"];
data_easy.hdr.chanunit = ["unknown";"unknown"];
data_easy.label = ["envelope";"pc1";"pc2";"pc3";"pc4"];
data_easy.fsample = Fs;
data_easy.sampleinfo = [];
data_easy.elec = [];
data_easy.cfg = [];

%% HARD
data_hard = data_easy;
for stim=1:length(pcas_hard)
    data_hard.time{stim} = [ts:ts:ts*length(pcas_hard{stim})];
    tmp_pc = pcas_hard{stim};
    tmp_env = envs_hard{stim};
    data_hard.trial{stim} = [tmp_env';tmp_pc(1:4,:)];
end

%% PLOT
cfg = [];
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';
cfg.taper = 'hanning';
cfg.foi = [0.5:Fs/L_min:20];
cfg.pad = ceil(L_max*ts);
freq_fourier_all_easy = ft_freqanalysis(cfg,data_easy);
freq_fourier_all_hard = ft_freqanalysis(cfg,data_hard);

envs_easy = squeeze(freq_fourier_all_easy.powspctrm(:,1,:));
envs_hard = squeeze(freq_fourier_all_hard.powspctrm(:,1,:));
pc1_easy = squeeze(freq_fourier_all_easy.powspctrm(:,2,:));
pc1_hard = squeeze(freq_fourier_all_hard.powspctrm(:,2,:));
pc2_easy = squeeze(freq_fourier_all_easy.powspctrm(:,3,:));
pc2_hard = squeeze(freq_fourier_all_hard.powspctrm(:,3,:));
pc3_easy = squeeze(freq_fourier_all_easy.powspctrm(:,4,:));
pc3_hard = squeeze(freq_fourier_all_hard.powspctrm(:,4,:));
pc4_easy = squeeze(freq_fourier_all_easy.powspctrm(:,5,:));
pc4_hard = squeeze(freq_fourier_all_hard.powspctrm(:,5,:));

subplot(3,2,[1,2])
shadedErrorBar(freq_fourier_all_easy.freq, mean(envs_easy).*freq_fourier_all_easy.freq, (std(envs_easy)/sqrt(size(envs_easy,1))).*freq_fourier_all_easy.freq, 'lineProps', '-r')
hold on
shadedErrorBar(freq_fourier_all_hard.freq, mean(envs_hard).*freq_fourier_all_easy.freq, (std(envs_hard)/sqrt(size(envs_hard,1))).*freq_fourier_all_easy.freq, 'lineProps', '-b')
legend("EASY", "HARD");
title("SE")
xlabel("frequency (Hz)")
ylabel("Power")
xlim([0.5,10])

subplot(3,2,3)
shadedErrorBar(freq_fourier_all_easy.freq, mean(pc1_easy).*freq_fourier_all_easy.freq, (std(pc1_easy)/sqrt(size(pc1_easy,1))).*freq_fourier_all_easy.freq, 'lineProps', '-r')
hold on
shadedErrorBar(freq_fourier_all_hard.freq, mean(pc1_hard).*freq_fourier_all_easy.freq, (std(pc1_hard)/sqrt(size(pc1_hard,1))).*freq_fourier_all_easy.freq, 'lineProps', '-b')
title("PC1")
xlabel("frequency (Hz)")
ylabel("Power")
xlim([0.5,10])

subplot(3,2,4)
shadedErrorBar(freq_fourier_all_easy.freq, mean(pc2_easy).*freq_fourier_all_easy.freq, (std(pc2_easy)/sqrt(size(pc2_easy,1))).*freq_fourier_all_easy.freq, 'lineProps', '-r')
hold on
shadedErrorBar(freq_fourier_all_hard.freq, mean(pc2_hard).*freq_fourier_all_easy.freq, (std(pc2_hard)/sqrt(size(pc2_hard,1))).*freq_fourier_all_easy.freq, 'lineProps', '-b')
title("PC2")
xlabel("frequency (Hz)")
ylabel("Power")
xlim([0.5,10])

subplot(3,2,5)
shadedErrorBar(freq_fourier_all_easy.freq, mean(pc3_easy).*freq_fourier_all_easy.freq, (std(pc3_easy)/sqrt(size(pc3_easy,1))).*freq_fourier_all_easy.freq, 'lineProps', '-r')
hold on
shadedErrorBar(freq_fourier_all_hard.freq, mean(pc3_hard).*freq_fourier_all_easy.freq, (std(pc3_hard)/sqrt(size(pc3_hard,1))).*freq_fourier_all_easy.freq, 'lineProps', '-b')
title("PC3")
xlabel("frequency (Hz)")
ylabel("Power")
xlim([0.5,10])

subplot(3,2,6)
shadedErrorBar(freq_fourier_all_easy.freq, mean(pc4_easy).*freq_fourier_all_easy.freq, (std(pc4_easy)/sqrt(size(pc4_easy,1))).*freq_fourier_all_easy.freq, 'lineProps', '-r')
hold on
shadedErrorBar(freq_fourier_all_hard.freq, mean(pc4_hard).*freq_fourier_all_easy.freq, (std(pc4_hard)/sqrt(size(pc4_hard,1))).*freq_fourier_all_easy.freq, 'lineProps', '-b')
title("PC4")
xlabel("frequency (Hz)")
ylabel("Power")
xlim([0.5,10])




