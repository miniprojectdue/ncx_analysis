restoredefaultpath;

projectRoot = fullfile(getenv('HOME'), 'Desktop', 'CM_Mini_Project');
repoRoot    = fullfile(projectRoot, 'Sex_Specific_Human_Electromechanics');
funcRoot    = fullfile(projectRoot, 'functions');
scriptRoot  = fullfile(projectRoot, 'scripts');
resultRoot  = fullfile(projectRoot, 'results');

cd(projectRoot);
addpath(genpath(repoRoot));
addpath(funcRoot);
addpath(scriptRoot);

disp(projectRoot)
disp(repoRoot)

which modelRunner
which getStartingState
which getCurrentsStructure
which model_ToRORd_Land_Female
which model_ToRORd_Land_Male