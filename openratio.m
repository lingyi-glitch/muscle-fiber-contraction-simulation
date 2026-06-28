function integration = openratio(RyR, Calast, TimeNow)
MinTime = 0.01;
timestart = 1000 * TimeNow;
[~, ~, ~, integration] = CalciumDynamics(RyR, Calast, MinTime, timestart, true);
end
