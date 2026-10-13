function replay_matlab(dataDir, outputDir)
%REPLAY_MATLAB Plot six vector PDFs from packaged CSVs; no experiments.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin < 1 || strlength(string(dataDir)) == 0
    dataDir = fullfile(root, 'data', 'publication');
end
if nargin < 2 || strlength(string(outputDir)) == 0
    outputDir = fullfile(root, 'outputs', ...
        "matlab_replay_" + string(datetime('now',Format='yyyyMMdd_HHmmss_SSS')));
end
assert(~isfolder(outputDir) && ~isfile(outputDir), ...
    'Output must be a NEW directory; existing files will never be overwritten.');
coverage = readtable(fullfile(dataDir,'coverage.csv'), TextType='string',Delimiter=',');
assert(height(coverage)==450, 'Expected 450 current coverage rows.');
for model = ["linear","plate","nonlinear"]
    c = coverage(coverage.model==model,:);
    assert(height(c)==150);
    if model=="nonlinear"
        assert(all(c.M==50000 & c.P==50000 & c.reference_N==500000000));
    else
        assert(all(c.M==30000));
    end
end
mkdir(outputDir);
plot_nested_M_paper_style(dataDir,["linear","plate","nonlinear"],outputDir);
plot_intervals_pdf_label_only(outputDir, ...
    fullfile(dataDir,'intervals_linear_nonlinear.csv'),fullfile(dataDir,'intervals_plate.csv'));
fprintf('Six publication PDFs replayed from saved CSVs; no new trials.\n');
end
