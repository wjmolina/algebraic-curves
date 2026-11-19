% Show the triangular Bernstein scores from the four-method comparison.
if ~exist('d_array', 'var')
    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    load(fullfile(repo_root, 'data', 'statistics', 'comparison_scores.mat'), 'd_array');
end
figure;
boxplot(d_array(:, 2), 'Labels', {'Triangular Bernstein'}, 'Colors', 'k', 'Whisker', 1000000000);
ylabel('Sørensen–Dice Coefficient');
