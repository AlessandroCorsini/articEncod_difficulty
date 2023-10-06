%% ADDING PATHS
addpath '../fieldtrip-20210311'
addpath './results'
addpath '../analisi'
addpath(pwd)
load 'chs_eeg_ok'
load('freq_template')
load('neighbours_ordered')

%%
bad_chs=[22,28,32,41,46];
se_delta = zeros(22,64);
se_theta = zeros(22,64);
pc1_delta = zeros(22,64);
pc1_theta = zeros(22,64);

for subj=1:22
    se_delta(subj,:) = Unq1_delta{1,subj}{1,1}(:,1) - Unq1_delta{1,subj}{1,1}(:,2);
    se_theta(subj,:) = Unq1_theta{1,subj}{1,1}(:,1) - Unq1_theta{1,subj}{1,1}(:,2);
    pc1_delta(subj,:) = Unq2_delta{1,subj}{1,1}(:,1) - Unq2_delta{1,subj}{1,1}(:,2);
    pc1_theta(subj,:) = Unq2_theta{1,subj}{1,1}(:,1) - Unq2_theta{1,subj}{1,1}(:,2);

end

se_delta(:,bad_chs) = [];
se_theta(:,bad_chs) = [];
pc1_delta(:,bad_chs) = [];
pc1_theta(:,bad_chs) = [];

%%
figure
bar([mean(mean(se_delta,2)), mean(mean(se_theta,2))])
hold on
errorbar([1], [mean(mean(se_delta,2))], [std(mean(se_delta,2))/sqrt(22)])
errorbar([2], [mean(mean(se_theta,2))], [std(mean(se_theta,2))/sqrt(22)])

figure
bar([mean(mean(pc1_delta,2)), mean(mean(pc1_theta,2))])
hold on
errorbar([1], [mean(mean(pc1_delta,2))], [std(mean(pc1_delta,2))/sqrt(22)])
errorbar([2], [mean(mean(pc1_theta,2))], [std(mean(pc1_theta,2))/sqrt(22)])




