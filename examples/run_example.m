function run_example()
%RUN_EXAMPLE Run a short smoke-test simulation from a clean checkout.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);

global Externalforcerate Gamma SwingRateFore Eta KbT
global Time2Stop MaxSwingDis RateReleaseADP StiffThin ThinSpacing IsometricForce
global StiffThick MotorStiff filename addForce kthick
global F01 F10 Rate01 Rate10 RPi SPi ATPhydrolysis CPi frequency

MotorStiff = 2.0;
StiffThick = 45813;
KbT = 4.14;
IsometricForce = 5.0;
ThinSpacing = 5.5;
StiffThin = 22725;
MaxSwingDis = 6.0;
Time2Stop = 0.02;
SwingRateFore = 1700;
RateReleaseADP = 300;
Eta = 4.15;
Externalforcerate = 5000;

Gamma = 200;
frequency = 40;
addForce = 100;
CPi = 5;

kthick = 45813;
F01 = 40;
F10 = F01;
Rate01 = 200;
Rate10 = 500;
RPi = 250;
SPi = 2555;
ATPhydrolysis = 150;

outputRoot = fullfile(repoRoot, 'examples', 'output');
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

filename = fullfile(outputRoot, 'Data4readRigid.txt');
datacreat(filename);

caseDir = fullfile(outputRoot, 'example_short_run');
if ~isfolder(caseDir)
    mkdir(caseDir);
end

rng(1);
main(caseDir);

fprintf('Example simulation finished. Output folder:\n%s\n', caseDir);
end

