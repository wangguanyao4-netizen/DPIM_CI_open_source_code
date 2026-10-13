function root = setup_dpim()
% Path setup only. Does not start Python, tests, or an experiment.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'matlab'));
end
