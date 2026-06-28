function main(DataDir)
global NumThinNode NumThickNode Time2Stop NumThinEle NumThickEle Externalforcerate filename ...
    cuowu addForce Emergency

[ThinElement, InitPosThin, ThickElement, InitPosThick] = Dataread(filename);
NumThinEle = size(ThinElement, 1);
NumThinNode = numel(InitPosThin);
NumThickEle = size(ThickElement, 1);
NumThickNode = numel(InitPosThick);
cuowu = 0;

RyR = zeros(1, 10);
RyR(1) = 1;

Calast = zeros(1, 26);
Calast(1) = 35;
Calast(2) = 1480;
Calast(3) = 35;
Calast(4) = 1480;
Calast(6) = 10;
Calast(7) = 10;
Calast(8) = 50;
Calast(9) = 50;
Calast(10) = 14300;
Calast(11) = 14300;
Calast(12) = 50;
Calast(13) = 50;
Calast(14:19) = 1000;
Calast(23) = 0;
Calast(24) = 5;
Calast(26) = 135;

Emergency = 1;
InitSK = InitSKDD(ThinElement, ThickElement);

MotorWorkState = zeros(NumThickNode, 1);
MotorBindSite = zeros(NumThickNode, 1);
BoundMotor4Thin = zeros(NumThinNode, 1);
ThinNodeState = zeros(NumThinNode, 1);

InitialBoundMotorNum = 5;
MotorWorkState(1:InitialBoundMotorNum) = 2;
for i = 1:InitialBoundMotorNum
    freeThinNode = find(ThinNodeState == 0);
    [~, localID] = min(abs(InitPosThin(freeThinNode) - InitPosThick(i)));
    bindSite = freeThinNode(localID);
    MotorBindSite(i) = bindSite;
    BoundMotor4Thin(bindSite) = i;
    ThinNodeState(bindSite) = 1;
end

SwingDis = zeros(NumThickNode, 1);
memory = zeros(NumThickNode, 1);

initPosition = [InitPosThin; InitPosThick];
Position = initPosition;
ForceThickFilament = 0;

if ~isfolder(DataDir)
    mkdir(DataDir);
end
FileName4work = fullfile(DataDir, ['Force=', num2str(addForce)]);
WriteFile4Work = fopen(FileName4work, 'w');
if WriteFile4Work == -1
    error('main:OutputOpenFailed', 'Cannot open %s for writing.', FileName4work);
end

for jjj = 1:50
    MinTime = 0.01;
    TimeNow = jjj * MinTime;
    [Activationrate, RyR, Calast] = Ca2CaT(RyR, Calast, MinTime, TimeNow);
end

Calast(24) = (140 / 223) * nnz(MotorWorkState == 2);
Calast(25) = (140 / 223) * nnz(MotorWorkState == -2);
Calast(23) = 140 * Activationrate - Calast(24) - Calast(25);
GammaRate =  openratio(RyR, Calast, TimeNow);

TimeNow = 0;
Counter = 0;
while TimeNow < Time2Stop
    Counter = Counter + 1;

    [MinTime, NewWorkState, NewBindSite, NewBoundMotor4Thin, NewThinState, ...
        NewSwingDis, Newmemory, RandomVar, MotorID, motorForce, ...
        motorForceOut, deltaSwingDisOut] = MonteCarlo(Position, MotorWorkState, ...
        MotorBindSite, BoundMotor4Thin, ThinNodeState, SwingDis, memory, ...
        ForceThickFilament, GammaRate);

    TimeNow = TimeNow + MinTime;

    % Calcium kinetics intentionally use the motor state from the beginning
    % of this reaction interval. The updated motor state is applied below.
    [Activationrate, RyR, Calast] = Ca2CaT(RyR, Calast, MinTime, TimeNow);
    Calast(24) = (140 / 223) * nnz(MotorWorkState == 2);
    Calast(25) = (140 / 223) * nnz(MotorWorkState == -2);
    Calast(23) = 140 * Activationrate - Calast(24) - Calast(25);
    GammaRate = openratio(RyR, Calast, TimeNow);

    MotorWorkState = NewWorkState;
    MotorBindSite = NewBindSite;
    BoundMotor4Thin = NewBoundMotor4Thin;
    ThinNodeState = NewThinState;
    SwingDis = NewSwingDis;
    memory = Newmemory;

    ExternalForce = min(Externalforcerate * TimeNow, addForce);
    MotorOffNum = nnz(MotorWorkState == 0) + nnz(MotorWorkState == 10);
    MotorOnNum = nnz(abs(MotorWorkState) == 1);
    MotorWorkNum = nnz(abs(MotorWorkState) == 2);

    if MotorWorkNum < 3
        Emergency = 10;
    else
        Emergency = 1;
    end

    [MotorSK, ForceMotorVector] = MotorSKDD(MotorBindSite, MotorWorkState, SwingDis, initPosition);
    SKDD = InitSK + MotorSK;
    outForce = -ForceMotorVector;
    outForce(NumThinNode) = outForce(NumThinNode) + ExternalForce;

    Displace = SKDD \ outForce;
    Position = initPosition + Displace;

    ForceThickFilament = sum(motorForce);

    fprintf(WriteFile4Work, '%g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t\n ', ...
        Counter, TimeNow, RandomVar, MotorID, motorForceOut, deltaSwingDisOut, ...
        Displace(NumThinNode), MotorOffNum, MotorOnNum, MotorWorkNum, ...
        ForceThickFilament, Activationrate);

    if Displace(NumThinNode) >= 1550
        cuowu = 1;
        if TimeNow < 0.1
            cuowu = 2;
        end
        break;
    end
end

fclose(WriteFile4Work);
end
