clear
close all
% cd 'D:/Entrainment__/final_functions_behavior/'

%% ADDING PATHS AND SETTING PARAMS
addpath '../fieldtrip-20210311'
addpath './results'
addpath '../analisi'
addpath(pwd)

%%% PARAMETERS %%%
PCs = [1];
TRIAL_CONCAT = true;
PLOT_DIFFICULTY = true;
PLOT_PERFORMANCE = false;
PLOT_BETWEEN_PCS = false;
N_SUBJ = 22;
CLUSTERS = false;
start_time = 0.5;
bad_chs=[22,28,32,41,46];
flows = [0.5:0.5:9.5];
fhighs = [1.5:0.5:10.5];
%flows=[0.5, 5];
%fhighs=[4, 7];
load 'chs_eeg_ok'
chs_eeg=chs_eeg_ok
load('freq_template')

%% AVERAGE PLOTS
if(PLOT_DIFFICULTY)
    unq1_subjects_allFreqs_easy=zeros(length(fhighs),N_SUBJ,length(chs_eeg));
    unq2_subjects_allFreqs_easy=zeros(length(fhighs),N_SUBJ,length(chs_eeg));
    rdn_subjects_allFreqs_easy=zeros(length(fhighs),N_SUBJ,length(chs_eeg));
    syn_subjects_allFreqs_easy=zeros(length(fhighs),N_SUBJ,length(chs_eeg));
    
    unq1_subjects_allFreqs_hard=zeros(length(fhighs),N_SUBJ,length(chs_eeg));
    unq2_subjects_allFreqs_hard=zeros(length(fhighs),N_SUBJ,length(chs_eeg));
    rdn_subjects_allFreqs_hard=zeros(length(fhighs),N_SUBJ,length(chs_eeg));
    syn_subjects_allFreqs_hard=zeros(length(fhighs),N_SUBJ,length(chs_eeg));
end

for pc=PCs
    Js = [1];
    if(PLOT_BETWEEN_PCS)
        Js = PCs(PCs > pc);
    end
    for j=Js
        for i=1:length(fhighs)
            
            fhigh=fhighs(i)
            flow=flows(i)
            
            %% FILE LOADING
            if(PLOT_BETWEEN_PCS)
                fname_redundancy = "REDUNDANCY_PCs_";
                fname_unq1 = "UNIQUE1_PCs_";
                fname_unq2 = "UNIQUE2_PCs_";
                fname_synergy = "SYNERGY_PCs_";
            else
                fname_redundancy = "REDUNDANCY_";
                fname_unq1 = "UNIQUE1_";
                fname_unq2 = "UNIQUE2_";
                fname_synergy = "SYNERGY_";
            end
            fname_offset = "offset" + string(start_time * 1000) + "ms_allFreqs_";
            fname_difficulty = "difficultySplit_";
            fname_performance = "performanceSplit_";
            fname_bandp = "trials_concat_eegBandp_" + string(flow) + "to" + string(fhigh) + "hz_filteredInputs.mat";
            if(PLOT_DIFFICULTY)
                clusters_rdn_subjects_difficulty = true(1,59);
                clusters_unq1_subjects_difficulty = true(1,59);
                clusters_unq2_subjects_difficulty = true(1,59);
                clusters_syn_subjects_difficulty = true(1,59);
                if CLUSTERS
                    %                     clusters_rdn_subjects_difficulty = (load("results\PID\statisticsCB\PC" + pc + "\stat_rdn_difficulty_offset500_" + string(flow) + "to" + string(fhigh) + "hz.mat").stat_rdn.mask);
                                        clusters_unq1_subjects_difficulty_delta = (load("results/PID/statisticsCB/PC" + pc + "/stat_unq1_difficulty_offset500_" + string(0.5) + "to" + string(4) + "hz.mat").stat_unq1.mask);
                    clusters_unq2_subjects_difficulty_delta = (load("results/PID/statisticsCB/PC" + pc + "/stat_unq2_difficulty_offset500_" + string(0.5) + "to" + string(4) + "hz.mat").stat_unq2.mask)';
                    clusters_unq2_subjects_difficulty_theta = (load("results/PID/statisticsCB/PC" + pc + "/stat_unq2_difficulty_offset500_" + string(5) + "to" + string(7) + "hz.mat").stat_unq2.mask)';
                    % clusters_syn_subjects_difficulty = (load("results\PID\statisticsCB\PC" + pc + "\stat_syn_difficulty_offset500_" + string(flow) + "to" + string(fhigh) + "hz.mat").stat_syn.mask);
                end
                rdn_subjects_difficulty = load("results/PID/" + fname_redundancy + fname_offset + fname_difficulty + fname_bandp).Rdn;
                unq1_subjects_difficulty = load("results/PID/" + fname_unq1 + fname_offset + fname_difficulty + fname_bandp).Unq1;
                unq2_subjects_difficulty = load("results/PID/" + fname_unq2 + fname_offset + fname_difficulty + fname_bandp).Unq2;
                syn_subjects_difficulty = load("results/PID/" + fname_synergy + fname_offset + fname_difficulty + fname_bandp).Syn;
            end
            %% UNQ1
            if(PLOT_DIFFICULTY)
                unq1_subjects_easy=zeros(N_SUBJ,64);
                unq1_subjects_hard=zeros(N_SUBJ,64);
                
                for sub=1:N_SUBJ
                    unq1_subjects_easy(sub,:)=unq1_subjects_difficulty{sub}{pc, j}(:,1)';
                    unq1_subjects_hard(sub,:)=unq1_subjects_difficulty{sub}{pc, j}(:,2)';
                end
                %% UNQ2
                unq2_subjects_easy=zeros(N_SUBJ,64);
                unq2_subjects_hard=zeros(N_SUBJ,64);
                
                for sub=1:N_SUBJ
                    unq2_subjects_easy(sub,:)=unq2_subjects_difficulty{sub}{pc, j}(:,1)';
                    unq2_subjects_hard(sub,:)=unq2_subjects_difficulty{sub}{pc, j}(:,2)';
                end
                %% SYNERGY
                syn_subjects_easy=zeros(N_SUBJ,64);
                syn_subjects_hard=zeros(N_SUBJ,64);
                
                for sub=1:N_SUBJ
                    syn_subjects_easy(sub,:)=syn_subjects_difficulty{sub}{pc, j}(:,1)';
                    syn_subjects_hard(sub,:)=syn_subjects_difficulty{sub}{pc, j}(:,2)';
                end
                %% REDUNDANCY
                rdn_subjects_easy=zeros(N_SUBJ,64);
                rdn_subjects_hard=zeros(N_SUBJ,64);
                
                for sub=1:N_SUBJ
                    rdn_subjects_easy(sub,:)=rdn_subjects_difficulty{sub}{pc, j}(:,1)';
                    rdn_subjects_hard(sub,:)=rdn_subjects_difficulty{sub}{pc, j}(:,2)';
                end
                
                rdn_subjects_easy(:,bad_chs)=[];
                rdn_subjects_hard(:,bad_chs)=[];
                syn_subjects_easy(:,bad_chs)=[];
                syn_subjects_hard(:,bad_chs)=[];
                unq1_subjects_easy(:,bad_chs)=[];
                unq1_subjects_hard(:,bad_chs)=[];
                unq2_subjects_easy(:,bad_chs)=[];
                unq2_subjects_hard(:,bad_chs)=[];
                
                unq1_subjects_allFreqs_easy(i,:,:)=unq1_subjects_easy;
                unq2_subjects_allFreqs_easy(i,:,:)=unq2_subjects_easy;
                rdn_subjects_allFreqs_easy(i,:,:)=rdn_subjects_easy;
                syn_subjects_allFreqs_easy(i,:,:)=syn_subjects_easy;
                
                unq1_subjects_allFreqs_hard(i,:,:)=unq1_subjects_hard;
                unq2_subjects_allFreqs_hard(i,:,:)=unq2_subjects_hard;
                rdn_subjects_allFreqs_hard(i,:,:)=rdn_subjects_hard;
                syn_subjects_allFreqs_hard(i,:,:)=syn_subjects_hard;
            end
        end
        %% PLOT
        cd 'figures/PID/'
        fname_fig = "";
        if(PLOT_BETWEEN_PCS)
            fname_fig = "PC" + string(pc) + "_" + string(j);
        end
        % UNIQUE1
        %mean over channels
        mean_unq1_subjects_allFreqs_easy=mean(mean(unq1_subjects_allFreqs_easy(:, :, clusters_unq1_subjects_difficulty),3),2);
        std_unq1_subjects_allFreqs_easy=std(mean(unq1_subjects_allFreqs_easy(:, :, clusters_unq1_subjects_difficulty), 3),0,2)/sqrt(N_SUBJ);
        mean_unq1_subjects_allFreqs_hard=mean(mean(unq1_subjects_allFreqs_hard(:, :, clusters_unq1_subjects_difficulty),3),2);
        std_unq1_subjects_allFreqs_hard=std(mean(unq1_subjects_allFreqs_hard(:, :, clusters_unq1_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
%         mean_unq1_subjects_theta_easy=mean(mean(unq1_subjects_allFreqs_easy(2, :, true(1, 59)),3),2);
%         std_unq1_subjects_theta_easy=std(mean(unq1_subjects_allFreqs_easy(2, :, true(1, 59)),3),0,2)/sqrt(N_SUBJ);
%         mean_unq1_subjects_theta_hard=mean(mean(unq1_subjects_allFreqs_hard(2, :, true(1, 59)),3),2);
%         std_unq1_subjects_theta_hard=std(mean(unq1_subjects_allFreqs_hard(2, :, true(1, 59)),3),0,2)/sqrt(N_SUBJ);
%         mean_unq1_subjects_delta_easy=mean(mean(unq1_subjects_allFreqs_easy(1, :, clusters_unq1_subjects_difficulty_delta),3),2);
%         std_unq1_subjects_delta_easy=std(mean(unq1_subjects_allFreqs_easy(1, :, clusters_unq1_subjects_difficulty_delta),3),0,2)/sqrt(N_SUBJ);
%         mean_unq1_subjects_delta_hard=mean(mean(unq1_subjects_allFreqs_hard(1, :, clusters_unq1_subjects_difficulty_delta),3),2);
%         std_unq1_subjects_delta_hard=std(mean(unq1_subjects_allFreqs_hard(1, :, clusters_unq1_subjects_difficulty_delta),3),0,2)/sqrt(N_SUBJ);
%         % max channel
        tmp=squeeze(mean(unq1_subjects_allFreqs_easy,2));
        [~,ch_max_easy]=find(tmp==max(max(tmp)));
        max_unq1_subjects_allFreqs_easy=mean(unq1_subjects_allFreqs_easy(:,:,ch_max_easy),2);
        std_max_unq1_subjects_allFreqs_easy=std(unq1_subjects_allFreqs_easy(:,:,ch_max_easy),0,2)/sqrt(N_SUBJ);
        tmp=squeeze(mean(unq1_subjects_allFreqs_hard,2));
        [~,ch_max_hard]=find(tmp==max(max(tmp)));
        max_unq1_subjects_allFreqs_hard=mean(unq1_subjects_allFreqs_hard(:,:,ch_max_hard),2);
        std_max_unq1_subjects_allFreqs_hard=std(unq1_subjects_allFreqs_hard(:,:,ch_max_hard),0,2)/sqrt(N_SUBJ);
        
%         figure
        % Plot easy
        %         shadedErrorBar((fhighs + flows)/2, mean_unq1_subjects_allFreqs_easy ,std_unq1_subjects_allFreqs_easy,'lineProps','-r','transparent',1)
        %         errorbar((fhighs + flows)/2, mean_unq1_subjects_allFreqs_easy ,std_unq1_subjects_allFreqs_easy)
        %         hold on
        %         shadedErrorBar((fhighs + flows)/2, mean_unq1_subjects_allFreqs_hard ,std_unq1_subjects_allFreqs_hard,'lineProps','-b','transparent',1)
        %         errorbar((fhighs + flows)/2, mean_unq1_subjects_allFreqs_hard ,std_unq1_subjects_allFreqs_hard)
        
        % errorbar((fhighs + flows)/2, max_PID_lagged_unq1_subjects_allFreqs_easy, std_max_PID_lagged_unq1_subjects_allFreqs_easy)
        %title("Mean and max across channels and subjects Unique 1 easy")
        ax = gca;
        ax.FontSize = 8;
        hold off
        %savefig("PC" + string(pc) + "\UNIQUE1_" + fname_fig + "EASY_ALLFREQS");
        
                figure
        %         % Plot hard
%         bar([mean_unq1_subjects_delta_easy-mean_unq1_subjects_delta_hard, mean_unq1_subjects_theta_easy-mean_unq1_subjects_theta_hard]);
% errorbar(1,mean_unq1_subjects_delta_easy-mean_unq1_subjects_delta_hard,std(mean(unq1_subjects_allFreqs_easy(1, :, clusters_unq1_subjects_difficulty_delta),3) - mean(unq1_subjects_allFreqs_hard(1, :, clusters_unq1_subjects_difficulty_delta),3))/sqrt(22));
%                 errorbar(2, mean_unq1_subjects_theta_easy-mean_unq1_subjects_theta_hard,std(mean(unq1_subjects_allFreqs_easy(2, :, clusters_unq1_subjects_difficulty),3) - mean(unq1_subjects_allFreqs_hard(2, :, clusters_unq1_subjects_difficulty),3))/sqrt(22))

errorbar((fhighs + flows)/2, [mean_unq1_subjects_delta_easy, mean_unq1_subjects_theta_easy] ,[std_unq1_subjects_delta_easy, std_unq1_subjects_theta_easy])
                hold on
                errorbar((fhighs + flows)/2, [mean_unq1_subjects_delta_hard, mean_unq1_subjects_theta_hard] ,[std_unq1_subjects_delta_hard, std_unq1_subjects_theta_hard])
        %         title("Mean and max across channels and subjects Unique 1 hard")
        %         ax = gca;
        %         ax.FontSize = 8;
        %         hold off
        %         %savefig("PC" + string(pc) + "\UNIQUE1_" + fname_fig + "HARD_ALLFREQS");
        
        % UNIQUE 2
        %mean over channels
        mean_unq2_subjects_allFreqs_easy=mean(mean(unq2_subjects_allFreqs_easy(:, :, clusters_unq2_subjects_difficulty),3),2);
        std_unq2_subjects_allFreqs_easy=std(mean(unq2_subjects_allFreqs_easy(:, :, clusters_unq2_subjects_difficulty), 3),0,2)/sqrt(N_SUBJ);
        mean_unq2_subjects_allFreqs_hard=mean(mean(unq2_subjects_allFreqs_hard(:, :, clusters_unq2_subjects_difficulty),3),2);
        std_unq2_subjects_allFreqs_hard=std(mean(unq2_subjects_allFreqs_hard(:, :, clusters_unq2_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
%         mean_unq2_subjects_theta_easy=mean(mean(unq2_subjects_allFreqs_easy(2, :, clusters_unq2_subjects_difficulty_theta),3),2);
%         std_unq2_subjects_theta_easy=std(mean(unq2_subjects_allFreqs_easy(2, :, clusters_unq2_subjects_difficulty_theta),3),0,2)/sqrt(N_SUBJ);
%         mean_unq2_subjects_theta_hard=mean(mean(unq2_subjects_allFreqs_hard(2, :, clusters_unq2_subjects_difficulty_theta),3),2);
%         std_unq2_subjects_theta_hard=std(mean(unq2_subjects_allFreqs_hard(2, :, clusters_unq2_subjects_difficulty_theta),3),0,2)/sqrt(N_SUBJ);
%         mean_unq2_subjects_delta_easy=mean(mean(unq2_subjects_allFreqs_easy(1, :, clusters_unq2_subjects_difficulty_delta),3),2);
%         std_unq2_subjects_delta_easy=std(mean(unq2_subjects_allFreqs_easy(1, :, clusters_unq2_subjects_difficulty_delta),3),0,2)/sqrt(N_SUBJ);
%         mean_unq2_subjects_delta_hard=mean(mean(unq2_subjects_allFreqs_hard(1, :, clusters_unq2_subjects_difficulty_delta),3),2);
%         std_unq2_subjects_delta_hard=std(mean(unq2_subjects_allFreqs_hard(1, :, clusters_unq2_subjects_difficulty_delta),3),0,2)/sqrt(N_SUBJ);
        %         % max channel
        tmp=squeeze(mean(unq2_subjects_allFreqs_easy,2));
        [~,ch_max_easy]=find(tmp==max(max(tmp)));
        max_unq2_subjects_allFreqs_easy=mean(unq2_subjects_allFreqs_easy(:,:,ch_max_easy),2);
        std_max_unq2_subjects_allFreqs_easy=std(unq2_subjects_allFreqs_easy(:,:,ch_max_easy),0,2)/sqrt(N_SUBJ);
        tmp=squeeze(mean(unq2_subjects_allFreqs_hard,2));
        [~,ch_max_hard]=find(tmp==max(max(tmp)));
        max_unq2_subjects_allFreqs_hard=mean(unq2_subjects_allFreqs_hard(:,:,ch_max_hard),2);
        std_max_unq2_subjects_allFreqs_hard=std(unq2_subjects_allFreqs_hard(:,:,ch_max_hard),0,2)/sqrt(N_SUBJ);
        
        figure
        % Plot easy
        %         shadedErrorBar((fhighs + flows)/2, mean_unq2_subjects_allFreqs_easy ,std_unq2_subjects_allFreqs_easy, 'lineProps','-r','transparent',1)
%                 bar([mean_unq2_subjects_delta_easy-mean_unq2_subjects_delta_hard, mean_unq2_subjects_theta_easy-mean_unq2_subjects_theta_hard]);
        %         shadedErrorBar((fhighs + flows)/2, mean_unq2_subjects_allFreqs_hard ,std_unq2_subjects_allFreqs_hard, 'lineProps','-b','transparent',1)

        errorbar((fhighs + flows)/2, [mean_unq2_subjects_delta_easy, mean_unq2_subjects_theta_easy] ,[std_unq2_subjects_delta_easy, std_unq2_subjects_theta_easy])
        hold on
        %         shadedErrorBar((fhighs + flows)/2, mean_unq2_subjects_allFreqs_hard ,std_unq2_subjects_allFreqs_hard, 'lineProps','-b','transparent',1)
        errorbar((fhighs + flows)/2, [mean_unq2_subjects_delta_hard, mean_unq2_subjects_theta_hard] ,[std_unq2_subjects_delta_hard, std_unq2_subjects_theta_hard])
        %errorbar((fhighs + flows)/2, max_PID_lagged_unq2_subjects_allFreqs_easy, std_max_PID_lagged_unq2_subjects_allFreqs_easy)
        title("Mean and max across channels and subjects Unique 2 easy")
        ax = gca;
        ax.FontSize = 8;
        savefig("PC" + string(pc) + "\UNIQUE2_" + fname_fig + "delta-theta");
        
        
        mean_unq1_delta = mean(unq1_subjects_allFreqs_easy(1,:,clusters_unq1_subjects_difficulty_delta), 3) - mean(unq1_subjects_allFreqs_hard(1,:,clusters_unq1_subjects_difficulty_delta), 3);
        mean_unq1_theta = mean(unq1_subjects_allFreqs_easy(2,:,true(1,59)), 3) - mean(unq1_subjects_allFreqs_hard(2,:,true(1,59)), 3);
        mean_unq2_delta = mean(unq2_subjects_allFreqs_easy(1,:,clusters_unq2_subjects_difficulty_delta), 3) - mean(unq2_subjects_allFreqs_hard(1,:,clusters_unq2_subjects_difficulty_delta), 3);
        mean_unq2_theta = mean(unq2_subjects_allFreqs_easy(2,:,clusters_unq2_subjects_difficulty_theta),3) - mean(unq2_subjects_allFreqs_hard(2,:,clusters_unq2_subjects_difficulty_theta), 3);

        std_unq1_delta = std(mean(unq1_subjects_allFreqs_easy(1,:,clusters_unq1_subjects_difficulty_delta), 3) - mean(unq1_subjects_allFreqs_hard(1,:,clusters_unq1_subjects_difficulty_delta), 3))/sqrt(22);
        std_unq1_theta = std(mean(unq1_subjects_allFreqs_easy(2,:,true(1,59)), 3) - mean(unq1_subjects_allFreqs_hard(2,:,true(1,59)), 3))/sqrt(22);
        std_unq2_delta = std(mean(unq2_subjects_allFreqs_easy(1,:,clusters_unq2_subjects_difficulty_delta), 3) - mean(unq2_subjects_allFreqs_hard(1,:,clusters_unq2_subjects_difficulty_delta), 3))/sqrt(22);
        std_unq2_theta = std(mean(unq2_subjects_allFreqs_easy(2,:,clusters_unq2_subjects_difficulty_theta), 3) - mean(unq2_subjects_allFreqs_hard(2,:,clusters_unq2_subjects_difficulty_theta), 3))/sqrt(22);

        %         figure
        %         % Plot hard
        %         errorbar((fhighs + flows)/2, mean_PID_lagged_unq2_subjects_allFreqs_hard ,std_PID_lagged_unq2_subjects_allFreqs_hard)
        %         hold on
        %         errorbar((fhighs + flows)/2, max_PID_lagged_unq2_subjects_allFreqs_hard, std_max_PID_lagged_unq2_subjects_allFreqs_hard)
        %         title("Mean and max across channels and subjects Unique 2 hard")
        %         ax = gca;
        %         ax.FontSize = 8;
        %         hold off
        %         %savefig("PC" + string(pc) + "\UNIQUE2_" + fname_fig + "HARD_ALLFREQS");
        
        % REDUNDANCY
        % mean over channels
        mean_rdn_subjects_allFreqs_easy=mean(mean(rdn_subjects_allFreqs_easy(:, :, clusters_rdn_subjects_difficulty),3),2);
        std_rdn_subjects_allFreqs_easy=std(mean(rdn_subjects_allFreqs_easy(:, :, clusters_rdn_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
        mean_rdn_subjects_allFreqs_hard=mean(mean(rdn_subjects_allFreqs_hard(:, :, clusters_rdn_subjects_difficulty),3),2);
        std_rdn_subjects_allFreqs_hard=std(mean(rdn_subjects_allFreqs_hard(:, :, clusters_rdn_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
        
        % max channel
        tmp=squeeze(mean(rdn_subjects_allFreqs_easy,2));
        [~,ch_max_easy]=find(tmp==max(max(tmp)));
        max_rdn_subjects_allFreqs_easy=mean(rdn_subjects_allFreqs_easy(:,:,ch_max_easy),2);
        std_max_rdn_subjects_allFreqs_easy=std(rdn_subjects_allFreqs_easy(:,:,ch_max_easy),0,2)/sqrt(N_SUBJ);
        tmp=squeeze(mean(rdn_subjects_allFreqs_hard,2));
        [~,ch_max_hard]=find(tmp==max(max(tmp)));
        max_rdn_subjects_allFreqs_hard=mean(rdn_subjects_allFreqs_hard(:,:,ch_max_hard),2);
        std_max_rdn_subjects_allFreqs_hard=std(rdn_subjects_allFreqs_hard(:,:,ch_max_hard),0,2)/sqrt(N_SUBJ);
        
        figure
        % Plot easy
        shadedErrorBar((fhighs + flows)/2, mean_rdn_subjects_allFreqs_easy ,std_rdn_subjects_allFreqs_easy,'lineProps','-r','transparent',1)
        %         errorbar((fhighs + flows)/2, mean_rdn_subjects_allFreqs_easy ,std_rdn_subjects_allFreqs_easy)
        hold on
        shadedErrorBar((fhighs + flows)/2, mean_rdn_subjects_allFreqs_hard ,std_rdn_subjects_allFreqs_hard, 'lineProps','-b','transparent',1)
        %         errorbar((fhighs + flows)/2, mean_rdn_subjects_allFreqs_hard ,std_rdn_subjects_allFreqs_hard)
        
        %errorbar((fhighs + flows)/2, max_PID_lagged_red_subjects_allFreqs_easy, std_max_PID_lagged_red_subjects_allFreqs_easy)
        title("Mean and max across channels and subjects Redundancy easy")
        ax = gca;
        ax.FontSize = 8;
        hold off
        %savefig("PC" + string(pc) + "\REDUNDANCY_" + fname_fig + "EASY_ALLFREQS");
        
        %         figure
        %         % Plot hard
        %         errorbar((fhighs + flows)/2, mean_PID_lagged_red_subjects_allFreqs_hard ,std_PID_lagged_red_subjects_allFreqs_hard)
        %         hold on
        %         errorbar((fhighs + flows)/2, max_PID_lagged_red_subjects_allFreqs_hard, std_max_PID_lagged_red_subjects_allFreqs_hard)
        %         title("Mean and max across channels and subjects Redundancy hard")
        %         ax = gca;
        %         ax.FontSize = 8;
        %         hold off
        %         %savefig("PC" + string(pc) + "\REDUNDANCY_" + fname_fig + "HARD_ALLFREQS");
        
        % SYNERGY
        %mean over channels
        mean_syn_subjects_allFreqs_easy=mean(mean(syn_subjects_allFreqs_easy(:, :, clusters_syn_subjects_difficulty),3),2);
        std_syn_subjects_allFreqs_easy=std(mean(syn_subjects_allFreqs_easy(:, :, clusters_syn_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
        mean_syn_subjects_allFreqs_hard=mean(mean(syn_subjects_allFreqs_hard(:, :, clusters_syn_subjects_difficulty),3),2);
        std_syn_subjects_allFreqs_hard=std(mean(syn_subjects_allFreqs_hard(:, :, clusters_syn_subjects_difficulty),3),0,2)/sqrt(N_SUBJ);
        
        % max channel
        tmp=squeeze(mean(syn_subjects_allFreqs_easy,2));
        [~,ch_max_easy]=find(tmp==max(max(tmp)));
        max_syn_subjects_allFreqs_easy=mean(syn_subjects_allFreqs_easy(:,:,ch_max_easy),2);
        std_max_syn_subjects_allFreqs_easy=std(syn_subjects_allFreqs_easy(:,:,ch_max_easy),0,2)/sqrt(N_SUBJ);
        tmp=squeeze(mean(syn_subjects_allFreqs_hard,2));
        [~,ch_max_hard]=find(tmp==max(max(tmp)));
        max_syn_subjects_allFreqs_hard=mean(syn_subjects_allFreqs_hard(:,:,ch_max_hard),2);
        std_max_syn_subjects_allFreqs_hard=std(syn_subjects_allFreqs_hard(:,:,ch_max_hard),0,2)/sqrt(N_SUBJ);
        
        figure
        % Plot easy
        shadedErrorBar((fhighs + flows)/2, mean_syn_subjects_allFreqs_easy ,std_syn_subjects_allFreqs_easy, 'lineProps','-r','transparent',1)
        %         errorbar((fhighs + flows)/2, mean_syn_subjects_allFreqs_easy ,std_syn_subjects_allFreqs_easy)
        hold on
        shadedErrorBar((fhighs + flows)/2, mean_syn_subjects_allFreqs_hard ,std_syn_subjects_allFreqs_hard, 'lineProps','-b','transparent',1)
        %         errorbar((fhighs + flows)/2, mean_syn_subjects_allFreqs_hard ,std_syn_subjects_allFreqs_hard)
        
        % errorbar((fhighs + flows)/2, max_PID_lagged_syn_subjects_allFreqs_easy, std_max_PID_lagged_syn_subjects_allFreqs_easy)
        title("Mean and max across channels and subjects Synergy easy")
        ax = gca;
        ax.FontSize = 8;
        hold on
        %savefig("PC" + string(pc) + "\SYNERGY_" + fname_fig + "EASY_ALLFREQS");
        
        %         figure
        %         % Plot hard
        %         errorbar((fhighs + flows)/2, mean_PID_lagged_syn_subjects_allFreqs_hard ,std_PID_lagged_syn_subjects_allFreqs_hard)
        %         hold on
        %         errorbar((fhighs + flows)/2, max_PID_lagged_syn_subjects_allFreqs_hard, std_max_PID_lagged_syn_subjects_allFreqs_hard)
        %         title("Mean and max across channels and subjects Synergy hard")
        %         ax = gca;
        %         ax.FontSize = 8;
        %         hold off
        %         %savefig("PC" + string(pc) + "\SYNERGY_" + fname_fig + "HARD_ALLFREQS");
        %         cd '..\..\'
    end
end

% %cd shuffled\all_circular
%
%
% try
%     load("ch_env_"+stat+string(pc))
%     c_idx_env=[];
%     for i=1:length(ch_env)
%         c_idx_env=[c_idx_env find(ch_env{i}==string(chs_eeg_ok))];
%     end
%
%     tmp=squeeze(mean(PID_lagged_unq1_subjects_chs_allFreqs(:,:,c_idx_env),2));
%     [~,ch_max_significative]=find (tmp==max(max(tmp)));
%     ch_max_significative=c_idx_env(ch_max_significative);
%     chs_eeg_ok(ch_max_significative)
%
%     figure
%     cluster_PID_lagged_env_subjects_chs_allFreqs=mean(mean(PID_lagged_unq1_subjects_chs_allFreqs(:,:,c_idx_env),3),2);
%     clusterstd_PID_lagged_env_subjects_chs_allFreqs=std(mean(PID_lagged_unq1_subjects_chs_allFreqs(:,:,c_idx_env),3),0,2)/sqrt(22);
%     errorbar(flow_all,cluster_PID_lagged_env_subjects_chs_allFreqs,clusterstd_PID_lagged_env_subjects_chs_allFreqs)
%     hold on
%     maxSignificative_PID_lagged_unq1_subjects_chs_allFreqs=mean(PID_lagged_unq1_subjects_chs_allFreqs(:,:,ch_max_significative),2);
%     maxSignificativestd_PID_lagged_unq1_subjects_chs_allFreqs=std(PID_lagged_unq1_subjects_chs_allFreqs(:,:,ch_max_significative),0,2)/sqrt(22);
%     errorbar(flow_all,maxSignificative_PID_lagged_unq1_subjects_chs_allFreqs,maxSignificativestd_PID_lagged_unq1_subjects_chs_allFreqs)
%     legend(["significative chs" "max ch"])
%     title ("PC"+string(pc)+" - envelope channels")
% end
%
% % figure
% % mean_PID_lagged_unq2_subjects_chs_allFreqs=mean(mean(PID_lagged_unq2_subjects_chs_allFreqs,3),2);
% % std_PID_lagged_unq2_subjects_chs_allFreqs=std(mean(PID_lagged_unq2_subjects_chs_allFreqs,3),0,2)/sqrt(22);
% % errorbar(flow_all,mean_PID_lagged_unq2_subjects_chs_allFreqs,std_PID_lagged_unq2_subjects_chs_allFreqs)
% % hold on
% % tmp=squeeze(mean(PID_lagged_unq2_subjects_chs_allFreqs,2));
% % [~,ch_max]=find(tmp==max(max(tmp)))
% % max_PID_lagged_unq2_subjects_chs_allFreqs=mean(PID_lagged_unq2_subjects_chs_allFreqs(:,:,ch_max),2);
% % maxstd_PID_lagged_unq2_subjects_chs_allFreqs=std(PID_lagged_unq2_subjects_chs_allFreqs(:,:,ch_max),0,2)/sqrt(22);
% % errorbar(flow_all,max_PID_lagged_unq2_subjects_chs_allFreqs,maxstd_PID_lagged_unq2_subjects_chs_allFreqs)
% % % cluster_PID_lagged_unq2_subjects_chs_allFreqs=mean(mean(PID_lagged_unq2_subjects_chs_allFreqs(:,:,chs_cluster),3),2)
% % % clusterstd_PID_lagged_unq2_subjects_chs_allFreqs=std(mean(PID_lagged_unq2_subjects_chs_allFreqs(:,:,chs_cluster),3),0,2)/sqrt(22)
% % % errorbar(flow_all,cluster_PID_lagged_unq2_subjects_chs_allFreqs,clusterstd_PID_lagged_unq2_subjects_chs_allFreqs)
% % legend("mean","max")%,"cluster mean")
% % title("Mean and max across channels and subjects PC"+string(pc))
% % % savefig("./figs/unique_pc"+string(pc)+"_allfreqs")
% % ax = gca;
% % ax.FontSize = 8;
%
% try
%     load("ch_unq_"+stat+string(pc))
%
%     c_idx_unq=[];
%     for i=1:length(ch_unq)
%         c_idx_unq=[c_idx_unq find(ch_unq{i}==string(chs_eeg_ok))];
%     end
%
%     tmp=squeeze(mean(PID_lagged_unq2_subjects_chs_allFreqs(:,:,c_idx_unq),2));
%     [~,ch_max_significative]=find (tmp==max(max(tmp)));
%     ch_max_significative=c_idx_unq(ch_max_significative);
%     chs_eeg_ok(ch_max_significative)
%
%     figure
%     cluster_PID_lagged_unq2_subjects_chs_allFreqs=mean(mean(PID_lagged_unq2_subjects_chs_allFreqs(:,:,c_idx_unq),3),2);
%     clusterstd_PID_lagged_unq2_subjects_chs_allFreqs=std(mean(PID_lagged_unq2_subjects_chs_allFreqs(:,:,c_idx_unq),3),0,2)/sqrt(22);
%     errorbar(flow_all,cluster_PID_lagged_unq2_subjects_chs_allFreqs,clusterstd_PID_lagged_unq2_subjects_chs_allFreqs)
%     hold on
%     maxSignificative_PID_lagged_unq2_subjects_chs_allFreqs=mean(PID_lagged_unq2_subjects_chs_allFreqs(:,:,ch_max_significative),2);
%     maxSignificativestd_PID_lagged_unq2_subjects_chs_allFreqs=std(PID_lagged_unq2_subjects_chs_allFreqs(:,:,ch_max_significative),0,2)/sqrt(22);
%     errorbar(flow_all,maxSignificative_PID_lagged_unq2_subjects_chs_allFreqs,maxSignificativestd_PID_lagged_unq2_subjects_chs_allFreqs)
%     legend(["significative chs" "max ch"])
%     title("PC"+string(pc)+" - unique channels")
% end
%
% % figure
% % mean_PID_lagged_syn_subjects_chs_allFreqs=mean(mean(PID_lagged_syn_subjects_chs_allFreqs,3),2);
% % std_PID_lagged_syn_subjects_chs_allFreqs=std(mean(PID_lagged_syn_subjects_chs_allFreqs,3),0,2)/sqrt(length(22));
% % errorbar(flow_all,mean_PID_lagged_syn_subjects_chs_allFreqs,std_PID_lagged_syn_subjects_chs_allFreqs)
% % hold on
% % tmp=squeeze(mean(PID_lagged_syn_subjects_chs_allFreqs,2));
% % [~,ch_max]=find(tmp==max(max(tmp)))
% % max_PID_lagged_syn_subjects_chs_allFreqs=mean(PID_lagged_syn_subjects_chs_allFreqs(:,:,ch_max),2);
% % maxstd_PID_lagged_syn_subjects_chs_allFreqs=std(PID_lagged_syn_subjects_chs_allFreqs(:,:,ch_max),0,2)/sqrt(22);
% % errorbar(flow_all,max_PID_lagged_syn_subjects_chs_allFreqs,maxstd_PID_lagged_syn_subjects_chs_allFreqs)
% % % cluster_PID_lagged_syn_subjects_chs_allFreqs=mean(mean(PID_lagged_syn_subjects_chs_allFreqs(:,:,chs_cluster),3),2)
% % % clusterstd_PID_lagged_syn_subjects_chs_allFreqs=std(mean(PID_lagged_syn_subjects_chs_allFreqs(:,:,chs_cluster),3),0,2)/sqrt(22)
% % % errorbar(flow_all,cluster_PID_lagged_syn_subjects_chs_allFreqs,clusterstd_PID_lagged_syn_subjects_chs_allFreqs)
% % legend("mean","max")%,"cluster mean")
% % title("Mean and max across channels and subjects synergy with PC"+string(pc))
% % % savefig("./figs/synergy_pc"+string(pc)+"_allfreqs")
% % ax = gca;
% % ax.FontSize = 8;
%
% try
%
%     load("ch_syn_"+stat+string(pc))
%
%     c_idx_syn=[];
%     for i=1:length(ch_syn)
%         c_idx_syn=[c_idx_syn find(ch_syn{i}==string(chs_eeg_ok))];
%     end
%
%     tmp=squeeze(mean(PID_lagged_syn_subjects_chs_allFreqs(:,:,c_idx_syn),2));
%     [~,ch_max_significative]=find (tmp==max(max(tmp)));
%     ch_max_significative=c_idx_syn(ch_max_significative);
%     chs_eeg_ok(ch_max_significative)
%
%     figure
%     cluster_PID_lagged_syn_subjects_chs_allFreqs=mean(mean(PID_lagged_syn_subjects_chs_allFreqs(:,:,c_idx_syn),3),2);
%     clusterstd_PID_lagged_syn_subjects_chs_allFreqs=std(mean(PID_lagged_syn_subjects_chs_allFreqs(:,:,c_idx_syn),3),0,2)/sqrt(22);
%     errorbar(flow_all,cluster_PID_lagged_syn_subjects_chs_allFreqs,clusterstd_PID_lagged_syn_subjects_chs_allFreqs)
%     hold on
%     maxSignificative_PID_lagged_syn_subjects_chs_allFreqs=mean(PID_lagged_syn_subjects_chs_allFreqs(:,:,ch_max_significative),2);
%     maxSignificativestd_PID_lagged_syn_subjects_chs_allFreqs=std(PID_lagged_syn_subjects_chs_allFreqs(:,:,ch_max_significative),0,2)/sqrt(22);
%     errorbar(flow_all,maxSignificative_PID_lagged_syn_subjects_chs_allFreqs,maxSignificativestd_PID_lagged_syn_subjects_chs_allFreqs)
%     legend(["significative chs" "max ch"])
%     title("PC"+string(pc)+" - synergistic channels")
% end
%
% % figure
% % mean_PID_lagged_red_subjects_chs_allFreqs=mean(mean(PID_lagged_red_subjects_chs_allFreqs,3),2);
% % std_PID_lagged_red_subjects_chs_allFreqs=std(mean(PID_lagged_red_subjects_chs_allFreqs,3),0,2)/sqrt(length(22));
% % errorbar(flow_all,mean_PID_lagged_red_subjects_chs_allFreqs,std_PID_lagged_red_subjects_chs_allFreqs)
% % hold on
% % tmp=squeeze(mean(PID_lagged_red_subjects_chs_allFreqs,2));
% % [~,ch_max]=find (tmp==max(max(tmp)))
% % max_PID_lagged_red_subjects_chs_allFreqs=mean(PID_lagged_red_subjects_chs_allFreqs(:,:,ch_max),2);
% % maxstd_PID_lagged_red_subjects_chs_allFreqs=std(PID_lagged_red_subjects_chs_allFreqs(:,:,ch_max),0,2)/sqrt(22);
% % errorbar(flow_all,max_PID_lagged_red_subjects_chs_allFreqs,maxstd_PID_lagged_red_subjects_chs_allFreqs)
% % % cluster_PID_lagged_red_subjects_chs_allFreqs=mean(mean(PID_lagged_red_subjects_chs_allFreqs(:,:,chs_cluster),3),2)
% % % clusterstd_PID_lagged_red_subjects_chs_allFreqs=std(mean(PID_lagged_red_subjects_chs_allFreqs(:,:,chs_cluster),3),0,2)/sqrt(22)
% % % errorbar(flow_all,cluster_PID_lagged_red_subjects_chs_allFreqs,clusterstd_PID_lagged_red_subjects_chs_allFreqs)
% % legend("mean","max")%,"cluster mean")
% % title("Mean and max across channels and subjects redundancy with PC"+string(pc))
% % % savefig("./figs/redundancy_pc"+string(pc)+"_allfreqs")
% % ax = gca;
% % ax.FontSize = 8;
%
% try
%     load("ch_red_"+stat+string(pc))
%
%     c_idx_red=[];
%     for i=1:length(ch_red)
%         c_idx_red=[c_idx_red find(ch_red{i}==string(chs_eeg_ok))];
%     end
%
%     tmp=squeeze(mean(PID_lagged_red_subjects_chs_allFreqs(:,:,c_idx_red),2));
%     [~,ch_max_significative]=find (tmp==max(max(tmp)))
%     ch_max_significative=c_idx_red(ch_max_significative);
%     chs_eeg_ok(ch_max_significative)
%
%     figure
%     cluster_PID_lagged_red_subjects_chs_allFreqs=mean(mean(PID_lagged_red_subjects_chs_allFreqs(:,:,c_idx_red),3),2);
%     clusterstd_PID_lagged_red_subjects_chs_allFreqs=std(mean(PID_lagged_red_subjects_chs_allFreqs(:,:,c_idx_red),3),0,2)/sqrt(22);
%     errorbar(flow_all,cluster_PID_lagged_red_subjects_chs_allFreqs,clusterstd_PID_lagged_red_subjects_chs_allFreqs)
%     hold on
%     maxSignificative_PID_lagged_red_subjects_chs_allFreqs=mean(PID_lagged_red_subjects_chs_allFreqs(:,:,ch_max_significative),2);
%     maxSignificativestd_PID_lagged_red_subjects_chs_allFreqs=std(PID_lagged_red_subjects_chs_allFreqs(:,:,ch_max_significative),0,2)/sqrt(22);
%     errorbar(flow_all,maxSignificative_PID_lagged_red_subjects_chs_allFreqs,maxSignificativestd_PID_lagged_red_subjects_chs_allFreqs)
%     legend(["significative chs" "max ch"])
%     title("PC"+string(pc)+" - redundant channels")
% end
