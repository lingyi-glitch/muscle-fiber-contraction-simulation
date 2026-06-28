function [Activationrate, RyRnow, Canow, integration] = CalciumDynamics(RyR, Calast, MinTime, timestart, withIntegration)
global frequency

if nargin < 5
    withIntegration = false;
end
if isempty(frequency) || frequency <= 0
    error('CalciumDynamics:BadFrequency', 'frequency must be a positive scalar.');
end

timestep1 = 10000;
timestep = 10001;
dt = 1000 * MinTime / timestep1;
t = timestart + dt * (0:timestep - 1);

v = -80 * ones(1, timestep);
tv = t - (1000 / frequency) * floor(frequency * t / 1000);
activeVoltage = tv <= 4;
v(activeVoltage) = 115 * exp(-1.5 * tv(activeVoltage)) - 80;

% RyR state variables.
aerfa1 = 0.2;
K = 4.5;
Vb = -20;
kc = 0.5 * aerfa1 * exp((v - Vb) / (8 * K));
kcm = 0.5 * aerfa1 * exp(-(v - Vb) / (8 * K));
kl = 0.002;
klm = 900;
f = 0.175;

C0 = zeros(1, timestep); C0(1) = RyR(1);
C1 = zeros(1, timestep); C1(1) = RyR(2);
C2 = zeros(1, timestep); C2(1) = RyR(3);
C3 = zeros(1, timestep); C3(1) = RyR(4);
C4 = zeros(1, timestep); C4(1) = RyR(5);
O0 = zeros(1, timestep); O0(1) = RyR(6);
O1 = zeros(1, timestep); O1(1) = RyR(7);
O2 = zeros(1, timestep); O2(1) = RyR(8);
O3 = zeros(1, timestep); O3(1) = RyR(9);
O4 = zeros(1, timestep); O4(1) = RyR(10);

for i = 1:timestep - 1
    dC0 = -kl * C0(i) + klm * O0(i) - 4 * kc(i) * C0(i) + kcm(i) * C1(i);
    dO0 = kl * C0(i) - klm * O0(i) - 4 * kc(i) * O0(i) / f + kcm(i) * O1(i) * f;
    dC1 = 4 * kc(i) * C0(i) - kcm(i) * C1(i) - kl * C1(i) / f + klm * f * O1(i) - 3 * kc(i) * C1(i) + 2 * kcm(i) * C2(i);
    dO1 = 4 * kc(i) * O0(i) / f - kcm(i) * O1(i) * f + kl * C1(i) / f - klm * f * O1(i) - 3 * kc(i) * O1(i) / f + 2 * kcm(i) * O2(i) * f;
    dC2 = 3 * kc(i) * C1(i) - 2 * kcm(i) * C2(i) - kl * C2(i) / (f * f) + klm * (f * f) * O2(i) - 2 * kc(i) * C2(i) + 3 * kcm(i) * C3(i);
    dO2 = 3 * kc(i) * O1(i) / f - 2 * kcm(i) * O2(i) * f + kl * C2(i) / (f * f) - klm * (f * f) * O2(i) - 2 * kc(i) * O2(i) / f + 3 * kcm(i) * O3(i) * f;
    dC3 = 2 * kc(i) * C2(i) - 3 * kcm(i) * C3(i) - kl * C3(i) / (f * f * f) + klm * (f * f * f) * O3(i) - kc(i) * C3(i) + 4 * kcm(i) * C4(i);
    dO3 = 2 * kc(i) * O2(i) / f - 3 * kcm(i) * O3(i) * f + kl * C3(i) / (f * f * f) - klm * (f * f * f) * O3(i) - kc(i) * O3(i) / f + 4 * kcm(i) * O4(i) * f;
    dC4 = kc(i) * C3(i) - 4 * kcm(i) * C4(i) - kl * C4(i) / (f * f * f * f) + klm * (f * f * f * f) * O4(i);
    dO4 = kc(i) * O3(i) / f - 4 * kcm(i) * O4(i) * f + kl * C4(i) / (f * f * f * f) - klm * (f * f * f * f) * O4(i);

    C0(i + 1) = C0(i) + dt * dC0;
    C1(i + 1) = C1(i) + dt * dC1;
    C2(i + 1) = C2(i) + dt * dC2;
    C3(i + 1) = C3(i) + dt * dC3;
    C4(i + 1) = C4(i) + dt * dC4;
    O0(i + 1) = O0(i) + dt * dO0;
    O1(i + 1) = O1(i) + dt * dO1;
    O2(i + 1) = O2(i) + dt * dO2;
    O3(i + 1) = O3(i) + dt * dO3;
    O4(i + 1) = O4(i) + dt * dO4;
end

D0 = zeros(1, timestep);
D1 = zeros(1, timestep);
D2 = zeros(1, timestep);
A1 = zeros(1, timestep);
A2 = zeros(1, timestep);
Ttot = 140 * ones(1, timestep);
T0 = zeros(1, timestep);
Mg1P = zeros(1, timestep);
ATP1 = zeros(1, timestep);
Mg2P = zeros(1, timestep);
ATP2 = zeros(1, timestep);
Mg1 = zeros(1, timestep);
Mg2 = zeros(1, timestep);
Mg1ATP = zeros(1, timestep);
Mg2ATP = zeros(1, timestep);
Ca1SR = zeros(1, timestep);
Ca2SR = zeros(1, timestep);
Ca1 = zeros(1, timestep);
Ca2 = zeros(1, timestep);
Ca1P = zeros(1, timestep);
Ca1ATP = zeros(1, timestep);
Cstot = 31000;
Ca1Cs = zeros(1, timestep);
Ca2Cs = zeros(1, timestep);
Ca2T = zeros(1, timestep);
Ca2CaT = zeros(1, timestep);
Ca2P = zeros(1, timestep);
Ca2ATP = zeros(1, timestep);

Ca1(1) = Calast(1); Ca1SR(1) = Calast(2); Ca2(1) = Calast(3); Ca2SR(1) = Calast(4);
Ca2T(1) = Calast(5); Ca1P(1) = Calast(6); Ca2P(1) = Calast(7);
Mg1P(1) = Calast(8); Mg2P(1) = Calast(9); Ca1Cs(1) = Calast(10); Ca2Cs(1) = Calast(11);
Ca1ATP(1) = Calast(12); Ca2ATP(1) = Calast(13); Mg1ATP(1) = Calast(14); Mg2ATP(1) = Calast(15);
ATP1(1) = Calast(16); ATP2(1) = Calast(17); Mg1(1) = Calast(18); Mg2(1) = Calast(19);
Ca2CaT(1) = Calast(20); D0(1) = Calast(21); D1(1) = Calast(22); D2(1) = Calast(23);
A1(1) = Calast(24); A2(1) = Calast(25); T0(1) = Calast(26);

if withIntegration
    integration = zeros(1, timestep);
    integration(1) = (dt / 1000) * (D2(1) + A2(1) + A1(1)) / 140;
else
    integration = [];
end

i2 = 60;
taoR = 0.75;
taoRSR = 0.75;
taoATP = 0.375;
taoMg = 1.5;
V0 = 0.864;
V1 = 0.01 * 0.95 * V0;
V2 = 0.99 * 0.95 * V0;
V1SR = 0.01 * 0.05 * V0;
V2SR = 0.99 * 0.05 * V0;
vSR = 2.4375;
KSR = 1;
Le = 0.00004;
Ptot = 1500;
KPon = 0;
KPoff = 0;
KCATPon = 0.15;
KCATPoff = 30;
KCsoff = 0.01;
KCson = KCsoff * 1.25 * 0.001;
KTon = 0.0885;
KToff = 0.115;
K0on = 0;
K0off = 0.15;
KMgon = 0;
KMgoff = 0;
KMATPon = 0.0015;
KMATPoff = 0.15;
KCaon = 0.15;
KCaoff = 0.05;

for i = 1:timestep - 1
    dCa1 = i2 * (O0(i) + O1(i) + O2(i) + O3(i) + O4(i)) * (Ca1SR(i) - Ca1(i)) / V1 - vSR * Ca1(i) / (V1 * (Ca1(i) + KSR)) + Le * (Ca1SR(i) - Ca1(i)) / V1 - taoR * (Ca1(i) - Ca2(i)) / V1 - (KPon * Ca1(i) * (Ptot - Ca1P(i) - Mg1P(i)) - KPoff * Ca1P(i)) - (KCATPon * Ca1(i) * ATP1(i) - KCATPoff * Ca1ATP(i));
    dCa1SR = -i2 * (O0(i) + O1(i) + O2(i) + O3(i) + O4(i)) * (Ca1SR(i) - Ca1(i)) / V1SR + vSR * Ca1(i) / (V1SR * (Ca1(i) + KSR)) - Le * (Ca1SR(i) - Ca1(i)) / V1SR - taoRSR * (Ca1SR(i) - Ca2SR(i)) / V1SR - (KCson * Ca1SR(i) * (Cstot - Ca1Cs(i)) - KCsoff * Ca1Cs(i));
    dCa2 = -vSR * Ca2(i) / (V2 * (Ca2(i) + KSR)) + Le * (Ca2SR(i) - Ca2(i)) / V2 + taoR * (Ca1(i) - Ca2(i)) / V2 - (KTon * Ca2(i) * T0(i) - KToff * Ca2T(i) + KTon * Ca2(i) * Ca2T(i) - KToff * Ca2CaT(i) + KTon * Ca2(i) * D0(i) - KToff * D1(i) + KTon * Ca2(i) * D1(i) - KToff * D2(i)) - (KPon * Ca2(i) * (Ptot - Ca2P(i) - Mg2P(i)) - KPoff * Ca2P(i)) - (KCATPon * Ca2(i) * ATP2(i) - KCATPoff * Ca2ATP(i));
    dCa2SR = vSR * Ca2(i) / (V2SR * (Ca2(i) + KSR)) - Le * (Ca2SR(i) - Ca2(i)) / V2SR + taoRSR * (Ca1SR(i) - Ca2SR(i)) / V2SR - (KCson * Ca2SR(i) * (Cstot - Ca2Cs(i)) - KCsoff * Ca2Cs(i));
    dCa2T = KTon * Ca2(i) * T0(i) - KToff * Ca2T(i) - KTon * Ca2(i) * Ca2T(i) + KToff * Ca2CaT(i) - K0on * Ca2T(i) + K0off * D1(i);
    dCa1P = KPon * Ca1(i) * (Ptot - Ca1P(i) - Mg1P(i)) - KPoff * Ca1P(i);
    dCa2P = KPon * Ca2(i) * (Ptot - Ca2P(i) - Mg2P(i)) - KPoff * Ca2P(i);
    dMg1P = KMgon * Mg1(i) * (Ptot - Ca1P(i) - Mg1P(i)) - KMgoff * Mg1P(i);
    dMg2P = KMgon * Mg2(i) * (Ptot - Ca2P(i) - Mg2P(i)) - KMgoff * Mg2P(i);
    dCa1Cs = KCson * Ca1SR(i) * (Cstot - Ca1Cs(i)) - KCsoff * Ca1Cs(i);
    dCa2Cs = KCson * Ca2SR(i) * (Cstot - Ca2Cs(i)) - KCsoff * Ca2Cs(i);
    dCa1ATP = KCATPon * Ca1(i) * ATP1(i) - KCATPoff * Ca1ATP(i) - taoATP * (Ca1ATP(i) - Ca2ATP(i)) / V1;
    dCa2ATP = KCATPon * Ca2(i) * ATP2(i) - KCATPoff * Ca2ATP(i) + taoATP * (Ca1ATP(i) - Ca2ATP(i)) / V2;
    dMg1ATP = KMATPon * Mg1(i) * ATP1(i) - KMATPoff * Mg1ATP(i) - taoATP * (Mg1ATP(i) - Mg2ATP(i)) / V1;
    dMg2ATP = KMATPon * Mg2(i) * ATP2(i) - KMATPoff * Mg2ATP(i) + taoATP * (Mg1ATP(i) - Mg2ATP(i)) / V2;
    dATP1 = -(KCATPon * Ca1(i) * ATP1(i) - KCATPoff * Ca1ATP(i)) - (KMATPon * Mg1(i) * ATP1(i) - KMATPoff * Mg1ATP(i)) - taoATP * (ATP1(i) - ATP2(i)) / V1;
    dATP2 = -(KCATPon * Ca2(i) * ATP2(i) - KCATPoff * Ca2ATP(i)) - (KMATPon * Mg2(i) * ATP2(i) - KMATPoff * Mg2ATP(i)) + taoATP * (ATP1(i) - ATP2(i)) / V2;
    dMg1 = -(KMgon * Mg1(i) * (Ptot - Ca1P(i) - Mg1P(i)) - KMgoff * Mg1P(i)) - (KMATPon * Mg1(i) * ATP1(i) - KMATPoff * Mg1ATP(i)) - taoMg * (Mg1(i) - Mg2(i)) / V1;
    dMg2 = -(KMgon * Mg2(i) * (Ptot - Ca2P(i) - Mg2P(i)) - KMgoff * Mg2P(i)) - (KMATPon * Mg2(i) * ATP2(i) - KMATPoff * Mg2ATP(i)) + taoMg * (Mg1(i) - Mg2(i)) / V2;
    dCa2CaT = KTon * Ca2(i) * Ca2T(i) - KToff * Ca2CaT(i) - KCaon * Ca2CaT(i) + KCaoff * D2(i);
    dD0 = -KTon * Ca2(i) * D0(i) + KToff * D1(i) + K0on * T0(i) - K0off * D0(i);
    dD1 = KTon * Ca2(i) * D0(i) - KToff * D1(i) + K0on * Ca2T(i) - K0off * D1(i) - KTon * Ca2(i) * D1(i) + KToff * D2(i);
    dD2 = KTon * Ca2(i) * D1(i) - KToff * D2(i) + KCaon * Ca2CaT(i) - KCaoff * D2(i);

    Ca1(i + 1) = Ca1(i) + dt * dCa1;
    Ca1SR(i + 1) = Ca1SR(i) + dt * dCa1SR;
    Ca2(i + 1) = Ca2(i) + dt * dCa2;
    Ca2SR(i + 1) = Ca2SR(i) + dt * dCa2SR;
    Ca2T(i + 1) = Ca2T(i) + dt * dCa2T;
    Ca1P(i + 1) = Ca1P(i) + dt * dCa1P;
    Ca2P(i + 1) = Ca2P(i) + dt * dCa2P;
    Mg1P(i + 1) = Mg1P(i) + dt * dMg1P;
    Mg2P(i + 1) = Mg2P(i) + dt * dMg2P;
    Ca1Cs(i + 1) = Ca1Cs(i) + dt * dCa1Cs;
    Ca2Cs(i + 1) = Ca2Cs(i) + dt * dCa2Cs;
    Mg1(i + 1) = Mg1(i) + dt * dMg1;
    Mg2(i + 1) = Mg2(i) + dt * dMg2;
    Ca2CaT(i + 1) = Ca2CaT(i) + dt * dCa2CaT;
    Ca1ATP(i + 1) = Ca1ATP(i) + dt * dCa1ATP;
    Ca2ATP(i + 1) = Ca2ATP(i) + dt * dCa2ATP;
    Mg1ATP(i + 1) = Mg1ATP(i) + dt * dMg1ATP;
    Mg2ATP(i + 1) = Mg2ATP(i) + dt * dMg2ATP;
    ATP1(i + 1) = ATP1(i) + dt * dATP1;
    ATP2(i + 1) = ATP2(i) + dt * dATP2;
    D0(i + 1) = D0(i) + dt * dD0;
    D1(i + 1) = D1(i) + dt * dD1;
    D2(i + 1) = D2(i) + dt * dD2;
    A1(i + 1) = A1(i);
    A2(i + 1) = A2(i);
    T0(i + 1) = Ttot(i + 1) - Ca2T(i + 1) - Ca2CaT(i + 1) - D0(i + 1) - D1(i + 1) - D2(i + 1) - A1(i + 1) - A2(i + 1);

    if withIntegration
        integration(i + 1) = integration(i) + (dt / 1000) * (D2(i + 1) + A2(i + 1) + A1(i + 1)) / 140;
    end
end

RyRnow = [C0(timestep), C1(timestep), C2(timestep), C3(timestep), C4(timestep), ...
    O0(timestep), O1(timestep), O2(timestep), O3(timestep), O4(timestep)];
Canow = [Ca1(timestep), Ca1SR(timestep), Ca2(timestep), Ca2SR(timestep), ...
    Ca2T(timestep), Ca1P(timestep), Ca2P(timestep), Mg1P(timestep), ...
    Mg2P(timestep), Ca1Cs(timestep), Ca2Cs(timestep), Ca1ATP(timestep), ...
    Ca2ATP(timestep), Mg1ATP(timestep), Mg2ATP(timestep), ATP1(timestep), ...
    ATP2(timestep), Mg1(timestep), Mg2(timestep), Ca2CaT(timestep), ...
    D0(timestep), D1(timestep), D2(timestep), A1(timestep), A2(timestep), T0(timestep)];
Activationrate = (D2(timestep) + A2(timestep) + A1(timestep)) / 140;
end
