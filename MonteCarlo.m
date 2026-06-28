function [MinTime, NewWorkState, NewBindSite, NewBoundMotor4Thin, NewThinState, ...
    NewSwingDis, Newmemory, RandomVar, MotorID, motorForce, NewmotorForceOut, ...
    NewdeltaSwingDisOut] = MonteCarlo(Position, MotorWorkState, MotorBindSite, ...
    BoundMotor4Thin, ThinNodeState, SwingDis, memory, ForceThickFilament, GammaRate)

global Gamma SwingRateFore Eta NumThickNode MaxSwingDis RateReleaseADP ...
    ThinSpacing IsometricForce KbT NumThinNode MotorStiff F01 F10 Rate01 Rate10 ...
    RPi SPi ATPhydrolysis CPi Emergency

Beta = 3 * ThinSpacing^2 / (2 * KbT);
TempRandomVar = zeros(NumThickNode, 1);

NewWorkState = MotorWorkState;
NewBindSite = MotorBindSite;
NewBoundMotor4Thin = BoundMotor4Thin;
NewThinState = ThinNodeState;
NewSwingDis = SwingDis;
Newmemory = memory;

tempmemory = memory;
TempWorkState = MotorWorkState;
TempBindSite = MotorBindSite;
TempBoundMotor4Thin = BoundMotor4Thin;
TempThinState = ThinNodeState;
TempSwingDis = SwingDis;

PosThin = Position(1:NumThinNode);
PosThick = Position(NumThinNode + 1:NumThinNode + NumThickNode);
TimeInterval = zeros(NumThickNode, 1);
motorForce = zeros(NumThickNode, 1);
motorForceOut = zeros(NumThickNode, 1);
deltaSwingDisOut = zeros(NumThickNode, 1);

for i = 1:NumThickNode
    if MotorWorkState(i) == 0
        RateOn = Rate01 * exp(ForceThickFilament / F01);
        TimeInterval(i) = -log(rand) / RateOn;
        TempRandomVar(i) = 0;
        TempWorkState(i) = 1;

    elseif MotorWorkState(i) == 10
        RateOn = Rate01 * exp(ForceThickFilament / F01);
        TimeInterval(i) = -log(rand) / RateOn;
        TempRandomVar(i) = 0;
        TempWorkState(i) = -1;

    elseif MotorWorkState(i) == -1
        FreeThinNode = find(ThinNodeState == 0);
        if isempty(FreeThinNode)
            error('MonteCarlo:NoFreeThinNode', 'No free thin-filament binding site is available.');
        end

        NormalDis = (PosThin(FreeThinNode(1)) - PosThick(i)) / ThinSpacing;
        SelectThin = FreeThinNode(1);
        Rate = bindingRate(Beta, Gamma, Emergency, NormalDis);
        selectnum = -log(rand) / Rate;

        for j = 1:numel(FreeThinNode)
            normalizedDistance = (PosThin(FreeThinNode(j)) - PosThick(i)) / ThinSpacing;
            Rate = bindingRate(Beta, Gamma, Emergency, normalizedDistance);
            tempnum = -log(rand) / Rate;
            if tempnum < selectnum
                SelectThin = FreeThinNode(j);
                selectnum = tempnum;
            end
        end

        selecttime = findPredictedBindingTime(selectnum, GammaRate);

        RateOff = 0.01 * Rate10 * exp(-ForceThickFilament / F10);
        offTime = -log(rand) / RateOff;
        if offTime < selecttime && i > 5
            selecttime = offTime;
            TempRandomVar(i) = 5;
            TimeInterval(i) = selecttime;
            TempWorkState(i) = 10;
        else
            TempBindSite(i) = SelectThin;
            TempThinState(SelectThin) = 1;
            TempSwingDis(i) = memory(i);
            TempWorkState(i) = 2;
            TempRandomVar(i) = -1;
            TimeInterval(i) = selecttime;
        end

    elseif MotorWorkState(i) == 1
        selecttime = -log(rand) / (ATPhydrolysis * Emergency);

        RateOff = Rate10 * exp(-ForceThickFilament / F10);
        offTime = -log(rand) / RateOff;
        if offTime < selecttime && i > 5
            selecttime = offTime;
            TempRandomVar(i) = 5;
            TimeInterval(i) = selecttime;
            TempWorkState(i) = 0;
        else
            TempWorkState(i) = -1;
            TempRandomVar(i) = -6;
            TimeInterval(i) = selecttime;
        end

    else
        Distance = PosThin(MotorBindSite(i)) - PosThick(i);
        motorForce(i) = MotorStiff * (SwingDis(i) + Distance);

        if SwingDis(i) < MaxSwingDis
            if motorForce(i) < IsometricForce
                SwingRate = SwingRateFore * exp((IsometricForce - motorForce(i)) / IsometricForce);
            else
                SwingRate = 0.01;
            end
            TimeSwing = -log(rand) / SwingRate;

            RateBondBreak = bondBreakRate(MotorWorkState(i), motorForce(i), Eta, Emergency);
            TimeBondBreak = -log(rand) / RateBondBreak;

            if MotorWorkState(i) == 2
                RatePi = RPi * exp(-CPi / SPi);
            else
                RatePi = 1000 * CPi / (CPi + 1900);
            end
            TimePi = -log(rand) / RatePi;

            if TimeSwing < TimeBondBreak && TimeSwing < TimePi
                DeltaSwingDis = (IsometricForce - motorForce(i)) / MotorStiff;
                if DeltaSwingDis > 0
                    motorForceOut(i) = motorForce(i);
                    deltaSwingDisOut(i) = DeltaSwingDis;
                    TempSwingDis(i) = min(MaxSwingDis, SwingDis(i) + DeltaSwingDis);
                else
                    TempSwingDis(i) = max(0, SwingDis(i));
                end
                TempRandomVar(i) = 2;
                TimeInterval(i) = TimeSwing;

            elseif TimePi < TimeBondBreak && TimePi < TimeSwing
                TempRandomVar(i) = 6;
                TimeInterval(i) = TimePi;
                TempWorkState(i) = -MotorWorkState(i);
                motorForceOut(i) = motorForce(i);

            else
                if MotorWorkState(i) == 2
                    sign = -1;
                    tempmemory(i) = SwingDis(i);
                else
                    sign = 1;
                    tempmemory(i) = 0;
                end
                TempWorkState(i) = sign;
                TempBindSite(i) = 0;
                TempBoundMotor4Thin(MotorBindSite(i)) = 0;
                TempThinState(MotorBindSite(i)) = 0;
                TempSwingDis(i) = 0;
                TempRandomVar(i) = 3 * sign;
                TimeInterval(i) = TimeBondBreak;
                motorForce(i) = 0;
            end

        else
            if MotorWorkState(i) == -2
                TimeReleaseADP = -log(rand) / RateReleaseADP;
            else
                TimeReleaseADP = 1000;
            end

            RateBondBreak = bondBreakRate(MotorWorkState(i), motorForce(i), Eta, Emergency);
            TimeBondBreak = -log(rand) / RateBondBreak;

            TimeInterval(i) = min(TimeReleaseADP, TimeBondBreak);
            TempWorkState(i) = 1;
            TempBindSite(i) = 0;
            TempBoundMotor4Thin(MotorBindSite(i)) = 0;
            TempThinState(MotorBindSite(i)) = 0;
            TempSwingDis(i) = 0;
            TempRandomVar(i) = 4;
            motorForce(i) = 0;
        end
    end
end

[MinTime, MotorID] = min(TimeInterval);
if MinTime < 0
    error('MonteCarlo:NegativeTime', 'Negative event time encountered.');
end

switch abs(TempRandomVar(MotorID))
    case 0
        NewWorkState(MotorID) = TempWorkState(MotorID);
        NewmotorForceOut = 0;
        NewdeltaSwingDisOut = 0;
    case 1
        NewThinState(TempBindSite(MotorID)) = TempThinState(TempBindSite(MotorID));
        NewBoundMotor4Thin(TempBindSite(MotorID)) = MotorID;
        NewWorkState(MotorID) = TempWorkState(MotorID);
        NewSwingDis(MotorID) = TempSwingDis(MotorID);
        NewBindSite(MotorID) = TempBindSite(MotorID);
        NewmotorForceOut = 0;
        NewdeltaSwingDisOut = 0;
    case 2
        NewSwingDis(MotorID) = TempSwingDis(MotorID);
        NewmotorForceOut = motorForceOut(MotorID);
        NewdeltaSwingDisOut = deltaSwingDisOut(MotorID);
    case 3
        NewThinState(MotorBindSite(MotorID)) = TempThinState(MotorBindSite(MotorID));
        NewBoundMotor4Thin(MotorBindSite(MotorID)) = TempBoundMotor4Thin(MotorBindSite(MotorID));
        NewWorkState(MotorID) = TempWorkState(MotorID);
        NewSwingDis(MotorID) = TempSwingDis(MotorID);
        NewBindSite(MotorID) = TempBindSite(MotorID);
        Newmemory(MotorID) = tempmemory(MotorID);
        NewmotorForceOut = 0;
        NewdeltaSwingDisOut = 0;
    case 4
        NewThinState(MotorBindSite(MotorID)) = TempThinState(MotorBindSite(MotorID));
        NewBoundMotor4Thin(MotorBindSite(MotorID)) = TempBoundMotor4Thin(MotorBindSite(MotorID));
        NewWorkState(MotorID) = TempWorkState(MotorID);
        NewSwingDis(MotorID) = TempSwingDis(MotorID);
        NewBindSite(MotorID) = TempBindSite(MotorID);
        NewmotorForceOut = 0;
        NewdeltaSwingDisOut = 0;
    case 5
        NewWorkState(MotorID) = TempWorkState(MotorID);
        NewmotorForceOut = 0;
        NewdeltaSwingDisOut = 0;
    case 6
        NewWorkState(MotorID) = TempWorkState(MotorID);
        NewSwingDis(MotorID) = TempSwingDis(MotorID);
        NewmotorForceOut = motorForceOut(MotorID);
        NewdeltaSwingDisOut = 0;
end

RandomVar = TempRandomVar(MotorID);
end

function Rate = bindingRate(Beta, Gamma, Emergency, normalizedDistance)
Rate = sqrt(Beta / pi) * 2 * Gamma * Emergency * ...
    exp(-Beta * normalizedDistance^2) / (1.0 + erf(0.3 * sqrt(Beta)));
end

function selecttime = findPredictedBindingTime(selectnum, GammaRate)
jjj = 0;
integral = 0;
while integral < selectnum
    jjj = jjj + 1;
    integral = GammaRate(jjj);
    if jjj > 10000
        break;
    end
end

if jjj >= 2
    selecttime = 1e-6 * (jjj - 1 + ...
        (selectnum - GammaRate(jjj - 1)) / (GammaRate(jjj) - GammaRate(jjj - 1)));
else
    selecttime = 1e-6 * (jjj - 1 + selectnum / GammaRate(jjj));
end
end

function RateBondBreak = bondBreakRate(MotorState, motorForce, Eta, Emergency)
if MotorState == 2
    RateBondBreak = (Eta / Emergency) * (40 * exp(-motorForce / 3) + 0.8 * exp(motorForce / 1.2));
    if motorForce > 5
        RateBondBreak = (Eta / Emergency) * (43 + 40 * exp(-motorForce / 3) + 3 * exp(motorForce / 5));
    end
else
    RateBondBreak = (Eta / Emergency) * (48 * exp(-motorForce / 1.5) + 2 * exp(motorForce / 5));
end
end
