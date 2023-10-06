clear
close all

%% ADDING PATHS AND SETTING PARAMS
addpath './results'
addpath(pwd)

%%% PARAMETERS %%%
PCs = [1,2,3,4];
N_SUBJ = 22;
CLUSTERS = false;
start_time = 0.5;
bad_chs=[22,28,32,41,46];
flows = [0.5:0.5:9.5];
fhighs = [1.5:0.5:10.5];
load 'chs_eeg_ok'
chs_eeg = chs_eeg_ok;
load('freq_template')

%% AVERAGE PLOTS
% Easy
env_subjects_allFreqs_easy = zeros(length(fhighs),N_SUBJ,length(chs_eeg));
pc1_subjects_allFreqs_easy = zeros(length(fhighs),N_SUBJ,length(chs_eeg));
pc2_subjects_allFreqs_easy = zeros(length(fhighs),N_SUBJ,length(chs_eeg));
pc3_subjects_allFreqs_easy = zeros(length(fhighs),N_SUBJ,length(chs_eeg));
pc4_subjects_allFreqs_easy = zeros(length(fhighs),N_SUBJ,length(chs_eeg));
%Hard
env_subjects_allFreqs_hard = zeros(length(fhighs),N_SUBJ,length(chs_eeg));
pc1_subjects_allFreqs_hard = zeros(length(fhighs),N_SUBJ,length(chs_eeg));
pc2_subjects_allFreqs_hard = zeros(length(fhighs),N_SUBJ,length(chs_eeg));
pc3_subjects_allFreqs_hard = zeros(length(fhighs),N_SUBJ,length(chs_eeg));
pc4_subjects_allFreqs_hard = zeros(length(fhighs),N_SUBJ,length(chs_eeg));

for i=1:length(fhighs)
    fhigh=fhighs(i)
    flow=flows(i)
    %% FILE LOADING
    fname_env = "MI_ENVELOPE_";
    fname_pcs = "MI_PCs_";
    fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_";
    fname_difficulty = "difficultySplit_";
    fname_bandp = "trials_concat_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";
    clusters_env_subjects_difficulty = true(1,59);
    clusters_pc1_subjects_difficulty = true(1,59);
    clusters_pc2_subjects_difficulty = true(1,59);
    clusters_pc3_subjects_difficulty = true(1,59);
    clusters_pc4_subjects_difficulty = true(1,59);
    if CLUSTERS
        %%%
    end
    env_subjects_difficulty = load("results/MI/" + fname_env + fname_offset + fname_difficulty + fname_bandp).Env;
    pcs_subjects_difficulty = load("results/MI/" + fname_pcs + fname_offset + fname_difficulty + fname_bandp).PCs;
    %% Env
        env_subjects_easy=zeros(N_SUBJ,64);
        env_subjects_hard=zeros(N_SUBJ,64);  
        for sub=1:N_SUBJ
            env_subjects_easy(sub,:)=env_subjects_difficulty{sub}(:,1)';
            env_subjects_hard(sub,:)=env_subjects_difficulty{sub}(:,2)';
        end
        %% PC1
        pc1_subjects_easy=zeros(N_SUBJ,64);
        pc1_subjects_hard=zeros(N_SUBJ,64);
        
        for sub=1:N_SUBJ
            pc1_subjects_easy(sub,:)=pcs_subjects_difficulty{sub}{1,1}(:,1)';
            pc1_subjects_hard(sub,:)=pcs_subjects_difficulty{sub}{1,1}(:,2)';
        end
        %% PC2
        pc2_subjects_easy=zeros(N_SUBJ,64);
        pc2_subjects_hard=zeros(N_SUBJ,64);
        
        for sub=1:N_SUBJ
            pc2_subjects_easy(sub,:)=pcs_subjects_difficulty{sub}{2,1}(:,1)';
            pc2_subjects_hard(sub,:)=pcs_subjects_difficulty{sub}{2,1}(:,2)';
        end
        %% PC3
        pc3_subjects_easy=zeros(N_SUBJ,64);
        pc3_subjects_hard=zeros(N_SUBJ,64);
        
        for sub=1:N_SUBJ
            pc3_subjects_easy(sub,:)=pcs_subjects_difficulty{sub}{3,1}(:,1)';
            pc3_subjects_hard(sub,:)=pcs_subjects_difficulty{sub}{3,1}(:,2)';
        end
        %% PC4
        pc4_subjects_easy=zeros(N_SUBJ,64);
        pc4_subjects_hard=zeros(N_SUBJ,64);
        
        for sub=1:N_SUBJ
            pc4_subjects_easy(sub,:)=pcs_subjects_difficulty{sub}{4,1}(:,1)';
            pc4_subjects_hard(sub,:)=pcs_subjects_difficulty{sub}{4,1}(:,2)';
        end
        env_subjects_easy(:,bad_chs)=[];
        pc1_subjects_easy(:,bad_chs)=[];
        pc2_subjects_easy(:,bad_chs)=[];
        pc3_subjects_easy(:,bad_chs)=[];
        pc4_subjects_easy(:,bad_chs)=[];
        env_subjects_hard(:,bad_chs)=[];
        pc1_subjects_hard(:,bad_chs)=[];
        pc2_subjects_hard(:,bad_chs)=[];
        pc3_subjects_hard(:,bad_chs)=[];
        pc4_subjects_hard(:,bad_chs)=[];
        % Easy
        env_subjects_allFreqs_easy(i,:,:)=env_subjects_easy;
        pc1_subjects_allFreqs_easy(i,:,:)=pc1_subjects_easy;
        pc2_subjects_allFreqs_easy(i,:,:)=pc2_subjects_easy;
        pc3_subjects_allFreqs_easy(i,:,:)=pc3_subjects_easy;
        pc4_subjects_allFreqs_easy(i,:,:)=pc4_subjects_easy;
        % Hard
        env_subjects_allFreqs_hard(i,:,:)=env_subjects_hard;
        pc1_subjects_allFreqs_hard(i,:,:)=pc1_subjects_hard;
        pc2_subjects_allFreqs_hard(i,:,:)=pc2_subjects_hard;
        pc3_subjects_allFreqs_hard(i,:,:)=pc3_subjects_hard;
        pc4_subjects_allFreqs_hard(i,:,:)=pc4_subjects_hard;
end
%% PLOT
cd 'figures/MI/'
fname_fig = "";
% ENV
%mean over channels
mean_env_subjects_allFreqs_easy=mean(mean(env_subjects_allFreqs_easy(:, :, clusters_env_subjects_difficulty),3),2);
std_env_subjects_allFreqs_easy=std(mean(env_subjects_allFreqs_easy(:, :, clusters_env_subjects_difficulty), 3),0,2)/sqrt(N_SUBJ);
mean_env_subjects_allFreqs_hard=mean(mean(env_subjects_allFreqs_hard(:, :, clusters_env_subjects_difficulty),3),2);
std_env_subjects_allFreqs_hard=std(mean(env_subjects_allFreqs_hard(:, :, clusters_env_subjects_difficulty), 3),0,2)/sqrt(N_SUBJ);
figure
shadedErrorBar((fhighs + flows)/2, mean_env_subjects_allFreqs_easy ,std_env_subjects_allFreqs_easy,'lineProps','-r','transparent',1)
hold on
shadedErrorBar((fhighs + flows)/2, mean_env_subjects_allFreqs_hard ,std_env_subjects_allFreqs_hard,'lineProps','-b','transparent',1)
legend("Easy","Hard")
title("Env")
saveas(gcf, 'Env.svg')

% PC1
mean_pc1_subjects_allFreqs_easy=mean(mean(pc1_subjects_allFreqs_easy(:, :, clusters_pc1_subjects_difficulty),3),2);
std_pc1_subjects_allFreqs_easy=std(mean(pc1_subjects_allFreqs_easy(:, :, clusters_pc1_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
mean_pc1_subjects_allFreqs_hard=mean(mean(pc1_subjects_allFreqs_hard(:, :, clusters_pc1_subjects_difficulty),3),2);
std_pc1_subjects_allFreqs_hard=std(mean(pc1_subjects_allFreqs_hard(:, :, clusters_pc1_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
figure
shadedErrorBar((fhighs + flows)/2, mean_pc1_subjects_allFreqs_easy ,std_pc1_subjects_allFreqs_easy,'lineProps','-r','transparent',1)
hold on
shadedErrorBar((fhighs + flows)/2, mean_pc1_subjects_allFreqs_hard ,std_pc1_subjects_allFreqs_hard,'lineProps','-b','transparent',1)
legend("Easy","Hard")
title("PC1")
saveas(gcf, 'PC1.svg')

% PC2
mean_pc2_subjects_allFreqs_easy=mean(mean(pc2_subjects_allFreqs_easy(:, :, clusters_pc2_subjects_difficulty),3),2);
std_pc2_subjects_allFreqs_easy=std(mean(pc2_subjects_allFreqs_easy(:, :, clusters_pc2_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
mean_pc2_subjects_allFreqs_hard=mean(mean(pc2_subjects_allFreqs_hard(:, :, clusters_pc2_subjects_difficulty),3),2);
std_pc2_subjects_allFreqs_hard=std(mean(pc2_subjects_allFreqs_hard(:, :, clusters_pc2_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
figure
shadedErrorBar((fhighs + flows)/2, mean_pc2_subjects_allFreqs_easy ,std_pc2_subjects_allFreqs_easy,'lineProps','-r','transparent',1)
hold on
shadedErrorBar((fhighs + flows)/2, mean_pc2_subjects_allFreqs_hard ,std_pc2_subjects_allFreqs_hard,'lineProps','-b','transparent',1)
legend("Easy","Hard")
title("PC2")
saveas(gcf, 'PC2.svg')

%PC3
mean_pc3_subjects_allFreqs_easy=mean(mean(pc3_subjects_allFreqs_easy(:, :, clusters_pc3_subjects_difficulty),3),2);
std_pc3_subjects_allFreqs_easy=std(mean(pc3_subjects_allFreqs_easy(:, :, clusters_pc3_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
mean_pc3_subjects_allFreqs_hard=mean(mean(pc3_subjects_allFreqs_hard(:, :, clusters_pc3_subjects_difficulty),3),2);
std_pc3_subjects_allFreqs_hard=std(mean(pc3_subjects_allFreqs_hard(:, :, clusters_pc3_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
figure
shadedErrorBar((fhighs + flows)/2, mean_pc3_subjects_allFreqs_easy ,std_pc3_subjects_allFreqs_easy,'lineProps','-r','transparent',1)
hold on
shadedErrorBar((fhighs + flows)/2, mean_pc3_subjects_allFreqs_hard ,std_pc3_subjects_allFreqs_hard,'lineProps','-b','transparent',1)
legend("Easy","Hard")
title("PC3")
saveas(gcf, 'PC3.svg')

%PC4
mean_pc4_subjects_allFreqs_easy=mean(mean(pc4_subjects_allFreqs_easy(:, :, clusters_pc4_subjects_difficulty),3),2);
std_pc4_subjects_allFreqs_easy=std(mean(pc4_subjects_allFreqs_easy(:, :, clusters_pc4_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
mean_pc4_subjects_allFreqs_hard=mean(mean(pc4_subjects_allFreqs_hard(:, :, clusters_pc4_subjects_difficulty),3),2);
std_pc4_subjects_allFreqs_hard=std(mean(pc4_subjects_allFreqs_hard(:, :, clusters_pc4_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
figure
shadedErrorBar((fhighs + flows)/2, mean_pc4_subjects_allFreqs_easy ,std_pc4_subjects_allFreqs_easy,'lineProps','-r','transparent',1)
hold on
shadedErrorBar((fhighs + flows)/2, mean_pc4_subjects_allFreqs_hard ,std_pc4_subjects_allFreqs_hard,'lineProps','-b','transparent',1)
legend("Easy","Hard")
title("PC4")
saveas(gcf, 'PC4.svg')
