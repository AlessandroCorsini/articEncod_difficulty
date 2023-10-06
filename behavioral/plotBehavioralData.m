clear
clc
%cd 'D:\Entrainment__\behavioral_\'

%% Single Subject
% distribution of trials based on RT, distribution of trials based on
% performance
%cd 'D:\Entrainment__\behavioral_\'
subs = load("results/all/behaviorAll.mat").behaviorAll.subs;
accuracy = load("results/all/behaviorAll.mat").behaviorAll.accuracy;
load("results/phrases/behaviorPhrases.mat");

for sub=1:length(subs)
    load("results\subs\" + subs{sub} + "_behavior.mat");
    validTrials = behavior.goodTrials == 1 & (behavior.answers == 1 | behavior.answers == 0);

    histogram(behavior.RTs(validTrials));
    title('Distribution of the trials based on RTs')
    xlabel('RT')
    ylabel('Quantity of Trials')
    %saveas(gcf, "figures\subs\" + subs{sub} + "_distributionRT.fig");

    for perf=1:size(behavior.performance, 1)
        figure
        histogram(behavior.performance{1,perf}(validTrials));
        hold off
        title('Distribution of trials based on performance')
        xlabel('Performance')
        ylabel('Quantity of Trials')
        %saveas(gcf, "figures\subs\" + subs{sub} + "_distributionPerformance");
    end

%     plot(behavior.RTs(behavior.goodTrials == 1 & (behavior.answers == 0 | behavior.answers == 1)), 'o-');
%     title('RTs')
%     xlabel('Trial');
%     ylabel('RT');
%     %saveas(gcf, "figures\subs\" + subs{sub} + "_RTs.fig");
    
    lengths = zeros(1, length(behaviorPhrases.lengths));
    for phrase=1:length(behaviorPhrases.lengths)
        lengths(behavior.phrases == phrase) = behaviorPhrases.lengths(phrase);
    end
    [phrases_ordered, idxs] = sort(lengths/1000);
    yyaxis left
    plot(phrases_ordered(behavior.goodTrials(idxs) == 1 & (behavior.answers(idxs) == 1 | behavior.answers(idxs) == 0)), 'o-')
    title("Phrase Length vs Trial Performance")
    hold on
    yyaxis right
    performance_ordered = behavior.performance(idxs);
    plot(performance_ordered(behavior.goodTrials(idxs) == 1 & (behavior.answers(idxs) == 1 | behavior.answers(idxs) == 0)), 'o-')
    hold off
    %saveas(gcf, "figures\subs\" + subs{sub} + "_length-performance");
    close(gcf)
end

%% All Subjects
% accuracy vs RT, scatter accuracy-RT, RT-answers
cd 'D:\Entrainment__\behavioral_\'
load("results\all\behaviorAll.mat");

[accuracy_ordered, idxs] = sort(behaviorAll.accuracy);
yyaxis left
p1 = plot(accuracy_ordered, 'ko-');
ylabel('Accuracy');
set(gca, 'YColor', 'k');
yyaxis right
p2 = errorbar(behaviorAll.meanRT(idxs), behaviorAll.stdRT(idxs), 'bo-');
ylabel('RTs');
title("Accuracy vs RTs");
xlabel('Subjects');
ylabel('RTs');
set(gca, 'YColor', 'b')
legend([p1(1), p2(1)],'Accuracy', 'RTs');
%hold off
%saveas(gcf, "figures\all\Accuracy-RT");
close(gcf);

Cov_meanRTAccuracy = cov([behaviorAll.meanRT, behaviorAll.accuracy]);
[V_meanRTAccuracy, D_meanRTAccuracy] = eig(Cov_meanRTAccuracy);
avg = [mean(behaviorAll.meanRT), mean(behaviorAll.accuracy)];
scatter(behaviorAll.meanRT, behaviorAll.accuracy);
hold on
line([avg(1) - V_meanRTAccuracy(1, 2)*2; avg(1) + V_meanRTAccuracy(1, 2)*2], [avg(2) - V_meanRTAccuracy(2, 2)*2; avg(2) + V_meanRTAccuracy(2, 2)*2]);
axis equal
title("Scatter RT-Accuracy")
xlabel('RT')
ylabel('Accuracy')
hold off
%saveas(gcf, "figures\all\Accuracy-RT_scatter");
close(gcf);


[meanRTCorrect_ordered, idxs] = sort(behaviorAll.meanRTCorrect);
errorbar(meanRTCorrect_ordered, behaviorAll.stdRTCorrect(idxs), 'x-');
hold on
errorbar(behaviorAll.meanRTWrong(idxs), behaviorAll.stdRTWrong(idxs), 'x-');
delta = (behaviorAll.meanRTWrong - behaviorAll.meanRTCorrect)./behaviorAll.meanRT;
bar(delta(idxs));
title("RTs for correct and wrong answers");
xlabel("Subject")
ylabel("RT(seconds)")
legend("RTs Correct", "RTs Wrong")
%saveas(gcf, "figures\all\RT-answer.fig")
close(gcf);

[delta_ordered, idxs] = sort(delta);
Cov_deltaRTAccuracy = cov([delta, behaviorAll.accuracy]);
[V_deltaRTAccuracy, D_deltaRTAccuracy] = eig(Cov_deltaRTAccuracy);
avg = [mean(delta), mean(behaviorAll.accuracy)];
scatter(delta_ordered, behaviorAll.accuracy(idxs))
hold on
axis equal
line([avg(1) - V_deltaRTAccuracy(1, 2)/5; avg(1) + V_deltaRTAccuracy(1, 2)/5], [avg(2) - V_deltaRTAccuracy(2, 2)/5; avg(2) + V_deltaRTAccuracy(2, 2)/5]);
title("Scatter deltaRT Accuracy")
xlabel("DeltaRT")
ylabel("Accuracy")
hold off
%saveas(gcf, "figures\all\deltaRT-Accuracy_scatter.fig");


%% Phrases
% orderedLength vs accuracy, orderedLength vs meanNormalizedRTWrong vs
% meanNormalizedRTCorrect vs meanNormalizedRTAll
cd 'D:\Entrainment__\behavioral_\'
load("results\phrases\behaviorPhrases.mat");

[lengths_ordered, idxs] = sort(behaviorPhrases.lengths);
yyaxis left
plot((lengths_ordered/1000), 'o-');
ylabel('Length')
hold on
yyaxis right
plot(behaviorPhrases.accuracy(idxs), 'o-');
ylabel('Accuracy')
hold off
title("Length vs Accuracy")
xlabel('Phrase')
%saveas(gcf, "figures\phrases\Length-Accuracy.fig")
close(gcf);

Cov_LengthAccuracy = cov([behaviorPhrases.lengths/1000, behaviorPhrases.accuracy], 'omitrows');
[V_LengthAccuracy, D_LengthAccuracy] = eig(Cov_LengthAccuracy);
avg = [mean(behaviorPhrases.lengths/1000, 'omitnan'), mean(behaviorPhrases.accuracy, 'omitnan')];
scatter(behaviorPhrases.lengths/1000, behaviorPhrases.accuracy);
hold on
line([avg(1) - V_LengthAccuracy(1, 2); avg(1) + V_LengthAccuracy(1, 2)], [avg(2) - V_LengthAccuracy(2, 2); avg(2) + V_LengthAccuracy(2, 2)]);
axis equal
title(" Scatter Length Accuracy")
xlabel('Lengths')
ylabel('Accuracy')
hold off
%saveas(gcf, "figures\phrases\Length-Accuracy_scatter");
close(gcf);

yyaxis left
plot((lengths_ordered/1000), 'o-');
ylabel('Length')
hold on
yyaxis right
plot(behaviorPhrases.meanNormalizedRTAll(idxs), 'o-');
ylabel('Mean Normalized RT')
hold off
title("Length vs NormalizedRT")
xlabel('Phrase')
%saveas(gcf, "figures\phrases\Length-RTNormalizedAll.fig")
close(gcf);

yyaxis left
plot((lengths_ordered/1000), 'o-');
ylabel('Length')
hold on
yyaxis right
plot(behaviorPhrases.meanNormalizedRTCorrect(idxs), 'o-');
ylabel('Mean Normalized RT')
hold off
title("Length vs NormalizedRT")
xlabel('Phrase')
legend("Length", "NormalizedRTCorrect")
%saveas(gcf, "figures\phrases\Length-RTNormalizedCorrect.fig")
close(gcf);

yyaxis left
plot((lengths_ordered/1000), 'o-');
ylabel('Length')
hold on
yyaxis right
plot(behaviorPhrases.meanNormalizedRTWrong(idxs), 'o-');
ylabel('Mean Normalized RT')
hold off
title("Length vs NormalizedRT")
xlabel('Phrase')
legend("Length", "NormalizedRTWrong")
%saveas(gcf, "figures\phrases\Length-RTNormalizedWrong.fig")
close(gcf);

[meanNormalizedRTCorrect_ordered, idxs] = sort(behaviorPhrases.meanNormalizedRTCorrect);
yyaxis left
plot(meanNormalizedRTCorrect_ordered, 'o-')
hold on
plot(behaviorPhrases.meanNormalizedRTWrong(idxs), 'o-');
delta = (behaviorPhrases.meanNormalizedRTWrong - behaviorPhrases.meanNormalizedRTCorrect)./behaviorPhrases.meanNormalizedRTAll;
plot(behaviorPhrases.accuracy(idxs), 'o-');
yyaxis right
ylim([0 4])
bar(behaviorPhrases.meanNormalizedRTWrong(idxs)./meanNormalizedRTCorrect_ordered);
legend("RT Correct", "RT Wrong", "Accuracy", "Delta");
hold off
close(gcf)


Cov_LengthRTCorrect = cov([behaviorPhrases.lengths/1000, behaviorPhrases.meanNormalizedRTCorrect], 'omitrows');
Cov_LengthRTWrong = cov([behaviorPhrases.lengths/1000, behaviorPhrases.meanNormalizedRTWrong], 'omitrows');
[V_LengthRTCorrect, D_LengthRTCorrect] = eig(Cov_LengthRTCorrect);
[V_LengthRTWrong, D_LengthRTWrong] = eig(Cov_LengthRTWrong);
avgCorrect = [mean(behaviorPhrases.lengths/1000, 'omitnan'), mean(behaviorPhrases.meanNormalizedRTCorrect, 'omitnan')];
avgWrong = [mean(behaviorPhrases.lengths/1000, 'omitnan'), mean(behaviorPhrases.meanNormalizedRTWrong, 'omitnan')];
scatter(behaviorPhrases.lengths/1000, behaviorPhrases.meanNormalizedRTCorrect, 'blue');
hold on
scatter(behaviorPhrases.lengths/1000, behaviorPhrases.meanNormalizedRTWrong, 'red');
line([avgCorrect(1) - V_LengthRTCorrect(1, 2); avgCorrect(1) + V_LengthRTCorrect(1, 2)], [avgCorrect(2) - V_LengthRTCorrect(2, 2); avgCorrect(2) + V_LengthRTCorrect(2, 2)], 'color', 'blue');
line([avgWrong(1) - V_LengthRTWrong(1, 2); avgWrong(1) + V_LengthRTWrong(1, 2)], [avgWrong(2) - V_LengthRTWrong(2, 2); avgWrong(2) + V_LengthRTWrong(2, 2)], 'color','red');
hold off
title('Scatter Length-RTCorrect&Wrong')
xlabel('Length (seconds)')
ylabel('Mean Normalized RT')
legend("Correct", "Wrong")
%saveas(gcf, "figures\phrases\Length-RTNormalizedCorrect&Wrong_scatter.fig");
close(gcf);

[accuracy_ordered, idxs] = sort(behaviorPhrases.accuracy);
yyaxis left
plot(accuracy_ordered, 'o-')
hold on
ylabel("Accuracy")
yyaxis right 
plot(behaviorPhrases.meanNormalizedRTAll(idxs), 'o-')
hold off
title("Accuracy vs RTNormalizedAll")
xlabel("Phrase")
ylabel("RTNormalizedAll")
%saveas(gcf, "figures\phrases\Accuracy-RTNormalizedAll.fig")
close(gcf)

Cov_AccuracyRTAll = cov([behaviorPhrases.accuracy, behaviorPhrases.meanNormalizedRTAll], 'omitrows');
[V_AccuracyRTAll, D_AccuracyRTAll] = eig(Cov_AccuracyRTAll);
avg = [mean(behaviorPhrases.accuracy, 'omitnan'), mean(behaviorPhrases.meanNormalizedRTAll, 'omitnan')];
scatter(behaviorPhrases.accuracy, behaviorPhrases.meanNormalizedRTAll);
hold on
line([avg(1) - V_AccuracyRTAll(1, 2)/5; avg(1) + V_AccuracyRTAll(1, 2)/5], [avg(2) - V_AccuracyRTAll(2, 2)/5; avg(2) + V_AccuracyRTAll(2, 2)/5]);
axis equal
title("Scatter Accuracy-NormalizedRT")
xlabel('Accuracy')
ylabel('RT All')
hold off
%saveas(gcf, "figures\phrases\Accuracy-RTNormalizedAll_scatter");
close(gcf);

yyaxis left
plot(accuracy_ordered, 'o-')
hold on
ylabel("Accuracy")
yyaxis right 
plot(behaviorPhrases.meanNormalizedRTCorrect(idxs), 'o-')
hold off
title("Accuracy vs RTNormalizedCorrect")
xlabel("Phrase")
ylabel("RTNormalizedCorrect")
%saveas(gcf, "figures\phrases\Accuracy-RTNormalizedCorrect.fig")
close(gcf)

Cov_AccuracyRTCorrect = cov([behaviorPhrases.accuracy, behaviorPhrases.meanNormalizedRTCorrect], 'omitrows');
[V_AccuracyRTCorrect, D_AccuracyRTCorrect] = eig(Cov_AccuracyRTCorrect);
avg = [mean(behaviorPhrases.accuracy, 'omitnan'), mean(behaviorPhrases.meanNormalizedRTCorrect, 'omitnan')];
scatter(behaviorPhrases.accuracy, behaviorPhrases.meanNormalizedRTCorrect);
hold on
line([avg(1) - V_AccuracyRTCorrect(1, 2)/5; avg(1) + V_AccuracyRTCorrect(1, 2)/5], [avg(2) - V_AccuracyRTCorrect(2, 2)/5; avg(2) + V_AccuracyRTCorrect(2, 2)/5]);
axis equal
title("Scatter Accuracy-NormalizedRTCorrect")
xlabel('Accuracy')
ylabel('RT Correct')
hold off
%saveas(gcf, "figures\phrases\Accuracy-RTNormalizedCorrect_scatter");
close(gcf);

yyaxis left
plot(accuracy_ordered, 'o-')
hold on
ylabel("Accuracy")
yyaxis right 
plot(behaviorPhrases.meanNormalizedRTWrong(idxs), 'o-')
hold off
title("Accuracy vs RTNormalizedWrong")
xlabel("Phrase")
ylabel("RTNormalizedWrong")
%saveas(gcf, "figures\phrases\Accuracy-RTNormalizedWrong.fig")
close(gcf)

Cov_AccuracyRTWrong = cov([behaviorPhrases.accuracy, behaviorPhrases.meanNormalizedRTWrong], 'omitrows');
[V_AccuracyRTWrong, D_AccuracyRTWrong] = eig(Cov_AccuracyRTWrong);
avg = [mean(behaviorPhrases.accuracy, 'omitnan'), mean(behaviorPhrases.meanNormalizedRTWrong, 'omitnan')];
scatter(behaviorPhrases.accuracy, behaviorPhrases.meanNormalizedRTWrong);
hold on
line([avg(1) - V_AccuracyRTWrong(1, 2)/5; avg(1) + V_AccuracyRTWrong(1, 2)/5], [avg(2) - V_AccuracyRTWrong(2, 2)/5; avg(2) + V_AccuracyRTWrong(2, 2)/5]);
axis equal
title("Scatter Accuracy-NormalizedRTWrong")
xlabel('Accuracy')
ylabel('RT Wrong')
hold off
%saveas(gcf, "figures\phrases\Accuracy-RTNormalizedWrong_scatter");
close(gcf);

delta = (behaviorPhrases.meanNormalizedRTWrong - behaviorPhrases.meanNormalizedRTCorrect)./behaviorPhrases.meanNormalizedRTAll;
yyaxis left
plot(accuracy_ordered, 'o-')
hold on
ylabel("Accuracy")
yyaxis right 
plot(delta(idxs), 'o-')
hold off
title("Accuracy vs RTNormalized Ratio Correct-Wrong")
xlabel("Phrase")
ylabel("RT Normalized Ratio Correct Wrong")
%saveas(gcf, "figures\phrases\Accuracy-RTNormalizedRatioCorrectWrong.fig")
close(gcf)

Cov_AccuracyRTRatioCorrectWrong = cov([behaviorPhrases.accuracy, delta], 'omitrows');
[V_AccuracyRTRatioCorrectWrong, D_AccuracyRTRatioCorrectWrong] = eig(Cov_AccuracyRTRatioCorrectWrong);
avg = [mean(behaviorPhrases.accuracy, 'omitnan'), mean(delta, 'omitnan')];
scatter(behaviorPhrases.accuracy, delta);
hold on
line([avg(1) - V_AccuracyRTRatioCorrectWrong(1, 2)/5; avg(1) + V_AccuracyRTRatioCorrectWrong(1, 2)/5], [avg(2) - V_AccuracyRTRatioCorrectWrong(2, 2)/5; avg(2) + V_AccuracyRTRatioCorrectWrong(2, 2)/5]);
title("Scatter Accuracy-RTNormalized (Ratio Correct/Wrong)")
xlabel('Accuracy')
ylabel('RT Ratio Correct Wrong')
hold off
%saveas(gcf, "figures\phrases\Accuracy-RTNormalizedRatioCorrectWrong_scatter");
close(gcf);

[accuracy_ordered, idxs] = sort(reshape(behaviorPhrases.accuracy4Words, [], 1));
yyaxis left
plot(accuracy_ordered, 'o-')
hold on
yyaxis right
meanNormalizedRT4WordsAll_reshaped = reshape(behaviorPhrases.meanNormalizedRT4WordsAll, [], 1);
plot(meanNormalizedRT4WordsAll_reshaped(idxs), 'o-')
hold off
%saveas(gcf,"figures\phrases\Accuracy-RTNormalizedAll4Words")
close(gcf);

Cov_AccuracyRT4Words = cov([reshape(behaviorPhrases.accuracy4Words, [], 1), meanNormalizedRT4WordsAll_reshaped], 'omitrows');
[V_AccuracyRT4Words, D_AccuracyRT4Words] = eig(Cov_AccuracyRT4Words);
avg = [mean(reshape(behaviorPhrases.accuracy4Words, [], 1), 'omitnan'), mean(meanNormalizedRT4WordsAll_reshaped, 'omitnan')];
scatter(reshape(behaviorPhrases.accuracy4Words, [], 1), meanNormalizedRT4WordsAll_reshaped);
hold on
line([avg(1) - V_AccuracyRT4Words(1, 1)/5; avg(1) + V_AccuracyRT4Words(1, 1)/5], [avg(2) - V_AccuracyRT4Words(2, 1)/5; avg(2) + V_AccuracyRT4Words(2, 1)/5]);
title("Scatter Accuracy-RTNormalized (Ratio Correct/Wrong)")
xlabel('Accuracy')
ylabel('RT Ratio Correct Wrong')
hold off
%saveas(gcf, "figures\phrases\Accuracy-RTNormalizedAll4Words_scatter");
close(gcf);

