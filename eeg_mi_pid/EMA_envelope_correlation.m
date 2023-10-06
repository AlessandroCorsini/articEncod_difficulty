clear
clc

addpath '../analisi'
addpath '../gcmi-master/matlab'

% 34, 29
d=dir('../analisi');
pca_dir = "../pca_new_nocut";
Fs = 400;
env_dir = "../envelope";
load('cut_times')
cut_times=round(cut_times*Fs*0.001); %from ms to samples

%% LOADING ENVELOPE
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

%% LOADING PCs
d_pca_em=dir(pca_dir);
pca_names={};
for i=4:length(d_pca_em)
    pca_names=[pca_names; d_pca_em(i).name];
end
pca_names=natsort(pca_names);
pcas={};
for i=1:length(pca_names);
    pca_tmp=load(pca_dir+"/"+pca_names{i});
    pcas=[pcas pca_tmp.pca_tmp];
end

%% CUTTING FEW SAMPLES (20to50) AT THE END OF ENVELOPES - they are little bit longer than PCs
for i=1:length(envelopes)
    %first compute the few samples in excess in speech respect to ema(or
    %pca) due to subsampling, and eliminate them from the end of envelopes
    diff= length(envelopes{i})-length(pcas{i});
    envelopes{i}=envelopes{i}(1:end-diff);
    % read the time from where to cut files to eliminate silences present
    % at the beginning of sentences
    start_cutTime=cut_times(i,1);
    end_cutTime=cut_times(i,2);
    envelopes{1,i}=envelopes{1,i}(start_cutTime:end-end_cutTime);
    pcas{i}=pcas{i}(:,start_cutTime:end-end_cutTime);
end

%% METAFILES LOADING
load("../behavioral_/results/phrases/behaviorPhrases.mat");
behaviorPhrases.accuracy([29,34]) = [];
[phrases_sorted, idx] = sort(behaviorPhrases.accuracy);

%% PREPARATION DATA
for i=1:length(pcas)
    pcas{1,i} = pcas{1,i}(1,:)';
end

% Define the filter parameters
flow = 0.5;     % Lower frequency cutoff
fhigh = 4;   % Higher frequency cutoff
fs = 300;    % Sampling rate of your data

% Design the filter
[b, a] = butter(2, [flow fhigh] / (fs/2), 'bandpass');

% Apply the filter to the data
for i=1:length(envelopes)
    envelopes{1,i} = filtfilt(b, a, envelopes{1,i});
    pcas{1,i} = filtfilt(b, a, pcas{1,i});
end

pcas([29,34]) = [];
envelopes([29,34]) = [];
env_sorted = envelopes(idx);
pcas_sorted = pcas(idx);
env_hard_concat = [];
pc_hard_concat = [];
env_easy_concat = [];
pc_easy_concat = [];
for i=2:length(env_sorted)/2
    env_hard_concat = [env_hard_concat; env_sorted{1,i}];
    pc_hard_concat = [pc_hard_concat; pcas_sorted{1,i}];
    env_easy_concat = [env_easy_concat; env_sorted{1,i+length(env_sorted)/2}];
    pc_easy_concat = [pc_easy_concat; pcas_sorted{1,i+length(env_sorted)/2}];
end
corrcoefs = zeros(1,48);
MI = zeros(1,48);
for i=1:48
    matrix_correlation = corrcoef(env_sorted{1,i}, pcas_sorted{1,i});
    corrcoefs(1,i) = matrix_correlation(1,2);
    MI(1,i) = gcmi_cc(copnorm(env_sorted{1,i}), copnorm(pcas_sorted{1,i}));
end

corrcoefs_easy = corrcoefs(idx(25:end));
corrcoefs_hard = corrcoefs(idx(1:24));
MI_easy = MI(idx(25:end));
MI_hard = MI(idx(1:24));

% Perform t-test
[h, p, ci, stats] = ttest2(corrcoefs_easy, corrcoefs_hard);

% Display the results
fprintf('t-statistic: %.4f\n', stats.tstat);
fprintf('p-value: %.4f\n', p);
fprintf('95%% confidence interval: [%.4f, %.4f]\n', ci(1), ci(2));

figure
plot(corrcoefs_easy, 'bo')
hold on
plot(corrcoefs_hard, 'ro')

% Perform t-test
[h, p, ci, stats] = ttest2(MI_easy, MI_hard);

% Display the results
fprintf('t-statistic: %.4f\n', stats.tstat);
fprintf('p-value: %.4f\n', p);
fprintf('95%% confidence interval: [%.4f, %.4f]\n', ci(1), ci(2));

figure
plot(MI_easy, 'bo')
hold on
plot(MI_hard, 'ro')
