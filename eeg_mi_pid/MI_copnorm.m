function [MI_env, MI_pc, MI_env_pc] = MI_copnorm(eeg,envelope,pca_k,do_pc,do_env,do_coupled)
% This function returns the value of the MI computed using thr GCMI estimator 
% between the specified continous variables and all the recorded EEG
% channels.
%%% Parameters
% 1. data_eeg: {(channels x 1) x samples} cells array containing all the channels of the
%   EEG in the rows and all the samples measures in the columns
% 2. envelope: (samples x 1) matrix containing the value of the speech
%   envelope at each sample point
% 3. pca_k: (samples x 1) matrix containing the values of the kth PC of
%   the EMA at each sample point
% 4. do_pc: boolean, true means the MI between the kth PC and the EEG
%   signals is computed
% 5. do_env: boolean, true means the MI between the speech envelope and the
%   EEG signals is computed
% 6. do_coupled: boolean, true means the MI between the multivariate signal
%   [speech envelope kth PC] and the EEG signals is computed 

    % The number of sample points
%     n_trial=length(data_eeg);

    % From class cells array to default matrix constructor operator
    % Note that data_eeg, enevelope and PC must have the same length in sample
%     % units
%     for trial=1:n_trial
%         if trial==1
%             % data_eeg{trial} means the (channels x 1) array containing the
%             % values of all the channels at sample "trial".
%             eeg=data_eeg{trial};
%             env=envelope{trial};
%             pc=pca_k{trial};
%         else
%             eeg=[eeg data_eeg{trial}];
%             env=[env; envelope{trial}];
%             pc=[pc; pca_k{trial}];        
%         end       
%     end
%     disp(length(eeg) + " samples");
    % eeg: (channels x samples) double matrix
    % env: (samples x 1) double matrix
    % pc: (samples x 1) double matrix

    % [env pc]: (samples x 2) double matrix containing in the first column the
    % envelope and the PC in the second column
    % Before the computation of the GCMI the continous data must be Copula
    % Normalized. The copula normalization is performed separately on each
    % column of the matrix. 
    %env_pc=copnorm([envelope pca_k]);
    if(do_env)
        env = copnorm(envelope);
    end
    if(do_pc)
        pc = copnorm(pca_k);
    end
    % Now env_pc, env and pc are copula normalized
    % At the end of the computation, each channel will have its own value of MI
    % in the three cases
    MI_env=zeros(64,1);
    MI_pc=zeros(64,1);
    MI_env_pc=zeros(64,1);
    
    % For each channel the values (in the three ways) of the MI using the GCMI
    % estimator is computed
    for ch=1:size(eeg,1)    
        if do_env MI_env(ch)=gcmi_cc(eeg(ch,:),env); end
        if do_pc MI_pc(ch)=gcmi_cc(eeg(ch,:),pc); end
        if do_coupled MI_env_pc(ch)=gcmi_cc(eeg(ch,:),[env pc]);end
    end
end

