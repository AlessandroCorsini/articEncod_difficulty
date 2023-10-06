function [Rdn,Unq1,Unq2,Syn] = PID_copnorm(eeg,envelope,pca_k)
% This function returns the value of the PID computed using the Iccs redundanct function
% for gaussian variables, for this reason all the input signals are copula normalized before
% the computation of the PID.
%% Parameters
% 1. data_eeg: {(channels x 1) x samples} cells array containing all the channels of the
%    EEG in the rows and all the samples measures in the columns
% 2. envelope: (samples x 1) matrix containing the value of the speech
%    envelope at each sample point
% 3. pca_k: (samples x 1) matrix containing the values of the kth PC of
%    the EMA at each sample point
%     rng(1000);

    lat=lattice2d();
 
    env=copnorm(envelope);
    pc=copnorm(pca_k);
    env_pc=[env pc];

    Rdn=zeros(64,1);
    Unq1=zeros(64,1);
    Unq2=zeros(64,1);
    Syn=zeros(64,1);

    for ch=1:size(eeg,1)
        C=cov([ env pc copnorm(eeg(ch,:)')]);
        lat_PI=calc_pi_mvn(lat,C, [1 1 1] ,@Iccs_mvn,true);
        pid=lat_PI.PI;
        Rdn(ch)=pid(1);
        Unq1(ch)=pid(2);
        Unq2(ch)=pid(3);
        Syn(ch)=pid(4);
    end
end

