function [Activationrate, RyRnow, Canow] = Ca2CaT(RyR, Calast, MinTime, TimeNow)
timestart = 1000 * (TimeNow - MinTime);
[Activationrate, RyRnow, Canow] = CalciumDynamics(RyR, Calast, MinTime, timestart, false);
end
