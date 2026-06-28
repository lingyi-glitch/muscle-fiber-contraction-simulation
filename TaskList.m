function TaskList()
global Externalforcerate Gamma SwingRateFore Eta KbT ...
    Time2Stop MaxSwingDis RateReleaseADP StiffThin ThinSpacing IsometricForce ...
    StiffThick MotorStiff filename addForce kthick ...
    F01 F10 Rate01 Rate10 RPi SPi ATPhydrolysis CPi frequency

warning('off');

rootDir = fileparts(mfilename('fullpath'));

MotorStiff = 2.0;
StiffThick = 45813;
KbT = 4.14;
IsometricForce = 5.0;
ThinSpacing = 5.5;
StiffThin = 22725;
MaxSwingDis = 6.0;
Time2Stop = 0.5;
SwingRateFore = 1700;
RateReleaseADP = 300;
Eta = 4.15;
Externalforcerate = 5000;

filename = fullfile(rootDir, 'Data4readRigid.txt');
datacreat(filename);

kthick = 45813;
F01 = 40;
F10 = F01;
Rate01 = 200;
Rate10 = 500;
RPi = 250;
SPi = 2555;
ATPhydrolysis = 150;
CPi = 5;

logFile = fullfile(rootDir, 'simulation_progress.log');
Processfid = fopen(logFile, 'a');
if Processfid == -1
    error('TaskList:LogOpenFailed', 'Cannot open %s for writing.', logFile);
end

for kk = 1:5
    for repeatID = 1:8
        Gamma = 200;
        frequency = 20 * kk;

        for cycleID = 8:8
            addForce = 12.5 * cycleID;
            DataDir = makeOutputDir(rootDir, repeatID, CPi, Gamma, frequency);
            runOneCase(DataDir, cycleID, addForce, Processfid);
        end

        for cycleID = 1:6
            addForce = 25 * cycleID + 100;
            DataDir = makeOutputDir(rootDir, repeatID, CPi, Gamma, frequency);
            runOneCase(DataDir, cycleID, addForce, Processfid);
        end
    end
end

fclose(Processfid);
 
 end

function DataDir = makeOutputDir(rootDir, repeatID, CPi, Gamma, frequency)
dirName = ['0', int2str(repeatID), '-CPi=', int2str(CPi), ...
    '-Gamma=', int2str(Gamma), 'frequency=', int2str(frequency)];
DataDir = fullfile(rootDir, dirName);
if ~isfolder(DataDir)
    mkdir(DataDir);
end
end

function runOneCase(DataDir, cycleID, addForce, Processfid)
global cuowu

main(DataDir);
trynumber = 1;

while cuowu == 2
    fprintf('Retry case with early excessive displacement.\n');
    main(DataDir);
    trynumber = trynumber + 1;
    if trynumber > 20
        break;
    end
end

if cuowu == 0
    message = sprintf('Cycle %2d, Force %g has been done!\n', cycleID, addForce);
else
    message = sprintf('Cycle %2d, Force %g has been done with excessive displacement!\n', cycleID, addForce);
end

fprintf('%s', message);
fprintf(Processfid, '%s', message);
end


