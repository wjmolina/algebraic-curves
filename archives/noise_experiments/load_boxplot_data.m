clear;
clc;

% Historical inputs and outputs are relative to the original experiment root.
experiment_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
load(fullfile(experiment_root, "statistics", "error_inf.mat"));
data = [];
data = [data; error_inf];
for j = 70 : - 5 : 30
    scores = load(fullfile(experiment_root, "statistics", "error_" + j + ".mat"));
    if isfield(scores, 'error_n') && isfield(scores, 'error_d')
        data = [data; scores.error_n; scores.error_d];
    else
        field = sprintf('error_%d', j);
        data = [data; scores.(field); scores.([field '_d'])];
    end
end
boxplot(data', "Labels", ["∞ dB", repelem((70 : - 5 : 30) + " dB", 2) + repmat(["", " d"], 1, 9)]);
ylim([0, 1]);
xtickangle(90);
h = gcf;
set(h, "PaperOrientation", "landscape");
set(h, "PaperUnits", "normalized");
set(h, "PaperPosition", [0, 0, 1, 1]);
print(gcf, "-dpdf", "test3.pdf");