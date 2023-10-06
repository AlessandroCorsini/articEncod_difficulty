% Verificare ulteriormente la correttezza
clear
clc
% cd 'D:\Entrainment__\behavioral_\'

load("results/all/behaviorAll.mat");
load("results/phrases/behaviorPhrases.mat");
[accuracy_ordered, idxs_phrase] = sort(behaviorPhrases.accuracy);
%disp(idxs_phrase);
for sub=1:length(behaviorAll.subs)
    load("results\subs\" + behaviorAll.subs{sub} + "_behavior.mat");
    validTrials = behavior.goodTrials == 1 & (behavior.answers == 1 | behavior.answers == 0);
    % Split based oon objective difficulty of the phrase
    easyTrials = zeros(length(behavior.goodTrials), 1);
    hardTrials = zeros(length(behavior.goodTrials), 1);
    hardTrials(validTrials & ismember(behavior.phrases, idxs_phrase(1:25))) = 1;
    easyTrials(validTrials & ismember(behavior.phrases, idxs_phrase(26:end))) = 1;
    %disp([sum(easyTrials), sum(hardTrials)]);
    % Split based on Performance of the subject indepemdent from the
    % difficulty of the phrase
    for split=1:length(behavior.performance)
        goods = zeros(sum(validTrials), 1);
        bads = zeros(sum(validTrials), 1);
        [performance_ordered, idxs] = sort(behavior.performance{split}(validTrials));    
        goods(idxs(1:floor(length(idxs)/2))) = 1;
        bads(idxs(floor((length(idxs)/2) + 1):end)) = 1;
        goodPerformanceTrials = zeros(length(behavior.goodTrials), 1);
        badPerformanceTrials = zeros(length(behavior.goodTrials), 1);
        goodPerformanceTrials(validTrials) = goods;
        badPerformanceTrials(validTrials) = bads;
        %disp([sum(goodPerformanceTrials), sum(badPerformanceTrials)]);
        %disp([mean(behavior.performance(goodPerformanceTrials == 1), 'omitnan'), mean(behavior.performance(badPerformanceTrials == 1), 'omitnan')]);
        dataSplit.goodPerformanceTrials{split} = goodPerformanceTrials;
        dataSplit.badPerformanceTrials{split} = badPerformanceTrials;
    end
    dataSplit.easyTrials = easyTrials;
    dataSplit.hardTrials = hardTrials;
    %save("results\subs\" + behaviorAll.subs{sub} + "_dataSplit.mat", "dataSplit");
end

lengths_easyTrials_subjects = [];
lengths_hardTrials_subjects = [];
lengths_goodPerformanceTrials_subjects = {};
lengths_badPerformanceTrials_subjects = {};
subs = behaviorAll.subs;
load("results\phrases\behaviorPhrases.mat");
for sub=1:length(subs)
    load("results\subs\" + subs{sub} + "_behavior.mat");
    load("results\subs\" + subs{sub} + "_dataSplit.mat");
    lengths_easyTrials_subjects = [lengths_easyTrials_subjects, sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.easyTrials == 1)))];
    lengths_hardTrials_subjects = [lengths_hardTrials_subjects, sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.hardTrials == 1)))];
    for split=1:length(behavior.performance)
        if(sub == 1)
            lengths_goodPerformanceTrials_subjects{split} = sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.goodPerformanceTrials{split} == 1)));
            lengths_badPerformanceTrials_subjects{split} = sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.badPerformanceTrials{split} == 1)));
        else
            lengths_goodPerformanceTrials_subjects{split} = [lengths_goodPerformanceTrials_subjects{split},  sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.goodPerformanceTrials{split} == 1)))];
            lengths_badPerformanceTrials_subjects{split} = [lengths_badPerformanceTrials_subjects{split},  sum(behaviorPhrases.lengths(behavior.phrases(dataSplit.badPerformanceTrials{split} == 1)))];
        end
    end
end

figure
plot(lengths_easyTrials_subjects/1000, '-o')
hold on
plot(lengths_hardTrials_subjects/1000, '-o')
legend("Easy", "Hard")
hold off
%savefig("figures\all\dataSplit\dataSplitDifficulty.fig");

for split=1:length(lengths_goodPerformanceTrials_subjects)
    figure
    plot(lengths_goodPerformanceTrials_subjects{split}/1000, '-o')
    hold on
    plot(lengths_badPerformanceTrials_subjects{split}/1000, '-o')
    legend("Good", "Bad")
    hold off
    %savefig("figures\all\dataSplit\dataSplitPerformance_split" + split + ".fig");

    figure
    commonTrials_goodEasy = zeros(1, length(subs));
    commonTrials_goodHard = zeros(1, length(subs));
    for sub=1:length(subs)
        load("results\subs\" + subs{sub} + "_dataSplit.mat");
        commonTrials_goodEasy(sub) = sum(dataSplit.goodPerformanceTrials{split} & dataSplit.easyTrials);
        commonTrials_goodHard(sub) = sum(dataSplit.goodPerformanceTrials{split} & dataSplit.hardTrials);
    end
    [~, idxs] = sort(commonTrials_goodEasy./commonTrials_goodHard);
    plot(commonTrials_goodEasy(idxs),'o-');
    hold on
    plot(commonTrials_goodHard(idxs), 'o-');
    hold off
    savefig("figures\all\dataSplit\commonTrials_good_split" + split + ".fig");
    
    figure
    goodTrials_answereRight = zeros(1, length(subs));
    goodTrials_answereWrong = zeros(1, length(subs));
    wrongAnswers = zeros(1, length(subs));
    for sub=1:length(subs)
        load("results\subs\" + subs{sub} + "_dataSplit.mat");
        load("results\subs\" + subs{sub} + "_behavior.mat");
        goodTrials_answereRight(sub) = sum(behavior.answers == 1 & (dataSplit.goodPerformanceTrials{split})');
        goodTrials_answereWrong(sub) = sum(behavior.answers == 0 & (dataSplit.goodPerformanceTrials{split})');
        wrongAnswers(sub) = sum(behavior.answers == 0);
    end
    [goodTrials_answereRight_ordered, idxs] = sort(goodTrials_answereRight);
    plot(goodTrials_answereRight_ordered, 'o-');
    hold on
    plot(goodTrials_answereWrong(idxs), 'o-');
    plot(wrongAnswers(idxs));
    hold off
   %savefig("figures\all\dataSplit\goodTrials_right_wrong_split" + split + ".fig");

end


