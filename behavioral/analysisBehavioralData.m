clear
clc
cd '/Volumes/My Passport/Entrainment__'

%% Find the subjects with complete and successful preprocessing
trigger_name = "startSound";
d=dir('analisi');
subs_files ={};
subs={};
subs_n={};

for i=3:length(d)
    if endsWith(d(i).name,"_filtered_epoched" + trigger_name + "_ica_interp-epo.fif")
        subs_files=[subs_files; d(i).name];
        tmp=split(d(i).name,'_f');
        tmp2=split(tmp(1),'_');
        subs=[subs; tmp(1)];
        subs_n=[subs_n; tmp2(2)];
    end
end

% No rhyme in one block of 34, corrupted pca in 29
INVALID_PHRASE = [34];
%% Single Subject
% SingleSubject (dim: #trials) --> answers, goodTrials, phrases, words, RT, performance
for sub=1:length(subs_n)
    goodTrials_old = load("analisi\" + subs{sub} + "_goodTrials.mat").good_trials + 1;
    meta = load("metadata\Subj"+subs_n{sub});
    numTrials = length(meta.got_ans);
    answers = meta.got_ans;
    goodTrials = zeros(1, numTrials);
    goodTrials(goodTrials_old) = 1;
    goodTrials(meta.all_sess(1:numTrials) == INVALID_PHRASE) = 0;
    phrases = (meta.all_sess).';
    words = zeros(1, numTrials);
    for phrase=1:size(meta.rand_rime_non_rime, 1)
        words(phrases == phrase) = meta.rand_rime_non_rime(phrase, :);
    end
    sum(goodTrials == 1 & phrases(1:numTrials) == INVALID_PHRASE)
    RTs = (meta.RTr).';
    behavior.answers = answers;
    behavior.goodTrials = goodTrials;
    behavior.phrases = phrases(1:numTrials);
    behavior.words = words(1:numTrials);
    behavior.RTs = RTs(1:numTrials);
    accuracy = sum(behavior.answers(behavior.goodTrials == 1) == 1) / (sum(behavior.answers(behavior.goodTrials == 1) == 1) + sum(behavior.answers(behavior.goodTrials == 1) == 0));
    alphas = [1, 1, 0.5, 0.75, 0.8];
    betas = [1, 2 * accuracy, 2, 1/0.75, 1/0.8];
    % Definition of different measures of performance
    % OLD:
    % performance = zeros(1, numTrials);
    % performance(behavior.answers == 1) = behavior.RTs(behavior.answers == 1);
    % performance(behavior.answers == 0) = behavior.RTs(behavior.answers == 0) * 2 * accuracy;
    performance{length(alphas)} = [];
    for alpha=1:length(alphas)
        performance{alpha} = zeros(1, numTrials);
        performance{alpha}(behavior.answers == 1) = behavior.RTs(behavior.answers == 1) * alphas(alpha);
        performance{alpha}(behavior.answers == 0) = behavior.RTs(behavior.answers == 0) * betas(alpha);
    end
    behavior.performance = performance;
    save("behavioral_\results\subs\" + subs{sub} + "_behavior.mat", "behavior");
end

%% All Subjects
% AllSubjects (dim:#subs) --> subs, accuracy, meanRT, stdRT
accuracy = zeros(length(subs), 1);
meanRT = zeros(length(subs), 1);
stdRT = zeros(length(subs), 1);
meanRTCorrect =  zeros(length(subs), 1);
stdRTCorrect = zeros(length(subs), 1);
meanRTWrong = zeros(length(subs), 1);
stdRTWrong = zeros(length(subs), 1);
for sub=1:length(subs)
    load("behavioral_\results\subs\" + subs{sub} + "_behavior.mat");
    accuracy(sub) = sum(behavior.answers(behavior.goodTrials == 1) == 1) / (sum(behavior.answers(behavior.goodTrials == 1) == 1) + sum(behavior.answers(behavior.goodTrials == 1) == 0));
    RTs = behavior.RTs(behavior.goodTrials == 1 & (behavior.answers == 1 | behavior.answers == 0));
    RTsCorrect = behavior.RTs(behavior.goodTrials == 1 & behavior.answers == 1);
    RTsWrong = behavior.RTs(behavior.goodTrials == 1 & behavior.answers == 0);
    meanRT(sub) = mean(RTs);
    stdRT(sub) = std(RTs);
    meanRTCorrect(sub) = mean(RTsCorrect);
    stdRTCorrect(sub) = std(RTsCorrect);
    meanRTWrong(sub) = mean(RTsWrong);
    stdRTWrong(sub) = std(RTsWrong);
end
clear behavior
behaviorAll.subs = subs;
behaviorAll.accuracy = accuracy;
behaviorAll.meanRT = meanRT;
behaviorAll.stdRT = stdRT;
behaviorAll.meanRTCorrect = meanRTCorrect;
behaviorAll.stdRTCorrect = stdRTCorrect;
behaviorAll.meanRTWrong = meanRTWrong;
behaviorAll.stdRTWrong = stdRTWrong;
save("behavioral_\results\all\behaviorAll.mat", "behaviorAll");


%% Phrases
% SinglePhrase (dim: #phrases) --> length, timesFirstWord, timesSecondWord, 
% accuracy, meanNormalizedRTCorrect, meanNormalizedRTWrong,
% meanNormalizedRTTotal, meanNoramlizedRT4WordsCorrect,
% meanNoramlizedRT4WordsWrong, meanNormalizedRT4WordsAll 

NUM_PHRASES = 50;
NUM_WORDS = 4;
lengths = [];
timesFirstWord = zeros(NUM_PHRASES, 2);
timesSecondWord = zeros(NUM_PHRASES, 2);
meanAccuracy = zeros(NUM_PHRASES, 1);
meanNormalizedRTCorrect = zeros(NUM_PHRASES, 1);
meanNormalizedRTWrong = zeros(NUM_PHRASES, 1);
meanNormalizedRTTotal = zeros(NUM_PHRASES, 1);
meanNormalizedRT2Words = zeros(NUM_PHRASES, 2);
stdNormalizedRT2Words = zeros(NUM_PHRASES, 2);
meanAccuracy2Words = zeros(NUM_PHRASES, 2);
stdAccuracy2Words = zeros(NUM_PHRASES, 2);

% Computation of the lengths of the phrases
cd wav_cut/
d = dir(pwd);
numsFiles = [];
for file=3:length(d)
    tmp = split(d(file).name, ".");
    tmp2 = split(tmp(1), "_");
    num = cell2mat(tmp2(2));
    num = str2double(num);
    numsFiles = [numsFiles; num];
end
[numsFiles_ordered, idxs] = sort(numsFiles);
% exclusion of the first 2 elements of the directory
idxs = idxs + 2;

data_allPhrases = {};
fs_allPhrases = [];
for i =1:length(idxs)
    [data, fs] = audioread(d(idxs(i)).name);
    data_allPhrases = [data_allPhrases data];
    fs_allPhrases = [fs_allPhrases fs];
    length_sample = length(data);
    % computing the length of the single phrase in ms
    lengths = [lengths; (length_sample/fs) * 1000];
end
lengths(INVALID_PHRASE) = nan;
cd ..\

% Definition of the times of the words
times_start = [3.075, 5.200;
    0.490, 1.130;
    2.965, 5.225;
    1.925, 3.775;
    0.685, 2.170;
    1.040, 5.135;
    0.100, 3.970;
    0.475, 4.355;
    0.680, 5.250;
    0.225, 3.635;
    1.000, 4.585;
    1.330, 5.350;
    1.700, 5.950;
    0.900, 5.210;
    1.300, 4.950;
    1.640, 2.630;
    1.820, 3.300;
    1.250, 3.140;
    0.800, 3.535;
    1.890, 3.790;
    1.215, 4.645;
    2.645, 4.655;
    1.070, 4.035;
    3.225, 4.535;
    3.865, 3.060;
    1.540, 4.830;
    1.230, 3.225;
    0.330, 4.850;
    2.285, 5.115;
    2.530, 5.130;
    0.640, 3.155;
    0.885, 3.910;
    1.050, 3.475;
    2.190, nan;
    0.065, 3.520;
    0.310, 7.980;
    0.675, 5.130;
    1.855, 6.030;
    0.420, 4.445;
    1.680, 4.815;
    1.070, 3.990;
    1.215, 5.825;
    0.825, 5.850;
    0.300, 2.180;
    0.835, 6.805;
    0.625, 3.735;
    2.140, 5.230;
    0.745, 3.500;
    0.630, 2.890;
    0.390, 3.990;];
    

times_end = [3.585, 5.680;
    0.980, 1.630;
    3.495, 5.800;
    2.395, 4.630;
    1.240, 2.520;
    1.740, 5.490;
    0.750, 4.460;
    1.365, 4.730;
    1.100, 5.600;
    0.535, 4.545;
    1.780, 5.160;
    1.860, 6.090;
    2.430, 6.240;
    1.670, 5.630;
    1.580, 5.275;
    2.270, 3.150;
    2.670, 3.620;
    1.580, 3.740;
    1.500, 3.945;
    2.340, 4.230;
    1.915, 5.300;
    3.015, 5.130;
    1.785, 4.340;
    3.990, 5.110;
    4.520, 3.865;
    1.860, 5.210;
    1.825, 3.640;
    0.655, 5.125;
    2.580, 5.565;
    3.020, 5.720;
    1.250, 3.695;
    1.445, 4.695;
    1.650, 4.100;
    2.550, nan;
    0.850, 4.340;
    0.660, 8.540;
    1.070, 5.600;
    2.300, 6.385;
    1.165, 5.085;
    2.305, 5.355;
    1.600, 4.315;
    1.750, 6.305;
    1.075, 6.280;
    0.630, 2.615;
    1.380, 7.520;
    1.095, 4.210;
    2.655, 5.725;
    1.120, 4.190;
    1.220, 3.460;
    0.990, 4.680;];
timesFirstWord(:, 1) = times_start(:, 1);
timesFirstWord(:, 2) = times_end(:, 1);
timesSecondWord(:, 1) = times_start(:, 2);
timesSecondWord(:, 2) = times_end(:, 2);

% Computation of the rest of the information
correctAns = zeros(NUM_PHRASES, 4);
wrongAns = zeros(NUM_PHRASES, 4);
correctNormalizedRT = zeros(NUM_PHRASES, 4);
wrongNormalizedRT = zeros(NUM_PHRASES, 4);

for sub=1:length(subs)
    load("behavioral_\results\subs\" + subs{sub} + "_behavior.mat");
    load("metadata\Subj" + subs_n{sub});
    validTrials = behavior.goodTrials == 1 & (behavior.answers == 1 | behavior.answers == 0);
    
    for phrase=1:NUM_PHRASES        
        for word=1:NUM_WORDS
            normalizedRT = behavior.RTs(validTrials & behavior.phrases == phrase & behavior.words == word)./mean(behavior.RTs(validTrials), 'omitnan'); % o max?
            if(behavior.answers(validTrials & behavior.phrases == phrase & behavior.words == word) == 1)
                correctAns(phrase, word) = correctAns(phrase, word) + 1;
                correctNormalizedRT(phrase, word) = correctNormalizedRT(phrase, word) + normalizedRT;
            elseif(behavior.answers(validTrials & behavior.phrases == phrase & behavior.words == word) == 0)
                wrongAns(phrase, word) = wrongAns(phrase, word) + 1;
                wrongNormalizedRT(phrase, word) = wrongNormalizedRT(phrase, word) + normalizedRT;
            end
        end
    end
end

clear behavior
behaviorPhrases.lengths = lengths;
behaviorPhrases.timesFirstWord = timesFirstWord;
behaviorPhrases.timesSecondWord = timesSecondWord;
behaviorPhrases.accuracy = sum(correctAns, 2)./sum(correctAns + wrongAns, 2);
behaviorPhrases.accuracy4Words = correctAns./(correctAns + wrongAns);
behaviorPhrases.meanNormalizedRTCorrect = sum(correctNormalizedRT, 2)./sum(correctAns, 2);
behaviorPhrases.meanNormalizedRTWrong = sum(wrongNormalizedRT, 2)./sum(wrongAns, 2);
behaviorPhrases.meanNormalizedRTAll = sum(correctNormalizedRT + wrongNormalizedRT, 2)./sum(correctAns + wrongAns, 2);
behaviorPhrases.meanNormalizedRT4WordsCorrect = correctNormalizedRT./correctAns;
behaviorPhrases.meanNormalizedRT4WordsWrong = wrongNormalizedRT./wrongAns;
behaviorPhrases.meanNormalizedRT4WordsAll = (correctNormalizedRT + wrongNormalizedRT)./(correctAns + wrongAns);
save("behavioral_\results\phrases\behaviorPhrases.mat", "behaviorPhrases");
