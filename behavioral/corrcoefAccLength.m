clear
clc

lengths = load("results/phrases/behaviorPhrases.mat").behaviorPhrases.lengths;
accuracies = load("results/PhrasesPooledSubjectsAccuracies.mat").accuracy;


InvalidPhs = [29,34];

lengths(InvalidPhs) = [];
accuracies(InvalidPhs) = [];

[r,p] = corrcoef(lengths, accuracies);