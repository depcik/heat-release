function EXRGYout = exergyHR(inp, calc, pdata, XCV, YCV, TCV, ECV, HR1out, HR2out, mu, mfa, mfl, mCV, dEHTfldt, dEindt, fuel, fuelchem)
% --- Routine specific variables
NM = inp(13);                               % Total number of chemical species
To = inp(84);                   % Exergy reference temperature [K]
po = inp(85);                   % Exergy reference pressure [Pa]
rps = calc(24);                 % Radians per second [rad/s] = [deg/s] * 3.14156 rad/180 deg
nps = calc(2);                              % Total number of data points
Twall = inp(82);                            % Wall temperature [K]
dt = calc(12);                              % Time-step of data [s] 
IVC = calc(3);                             % Index on IVC [-]
EVO = calc(4);                              % Index on EVO [-]
dps = calc(11);                             % Degrees per second [deg/s] = [rev/min]*[360 deg/rev]*[min/60 s]
MO2 = 2*inp(102);                           % Molecular mass of O2 [kg/mol]
MN2 = 2*inp(103);                           % Molecular mass of N2 [kg/mol]
MAr = inp(105);                             % Molecular mass of Ar [kg/mol]
MCO2 = inp(101) + MO2;                      % Molecular mass of CO2 [kg/mol]
MH2O = 2*inp(104) + inp(102);               % Molecular mass of H2O [kg/mol]
MCO = inp(101) + inp(102);      % Molecular mass of CO [kg/mol]
MOH = inp(102) + inp(104);      % Molecular mass of OH [kg/mol]
MH2 = 2*inp(104);               % Molecular mass of H2 [kg/mol]
MH = inp(104);                  % Molecular mass of H [kg/mol]
MO = inp(102);                  % Molecular mass of O [kg/mol]
% index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O, NM+1-Mixture
MM = [MO2 MN2 MAr MCO2 MH2O MCO MH2 MOH MH MO]; % Molecular mass array [kg/mol]
Mfa = fuel(24);                             % Molecular mass of added gas [kg/mol]
Mfl = fuel(5);                              % Molecular mass of liquid fuel [kg/mol]
Runiv = inp(100);                           % Universal gas constant [J/molK]
pref = inp(83);                             % CHEMKIN reference pressure [Pa]

% Orhter variables
% data(:,7) = dV/d(theta) [m3/rad] 
% CHEMKIN Curve fits
cd 'CHEMKIN III Data';
fidCO2=fopen('CO2.txt','r');
fidN2=fopen('N2.txt','r');
fidO2=fopen('O2.txt','r');
fidH2O=fopen('H2O.txt','r');
fidAr=fopen('Ar.txt','r');
fidCO=fopen('CO.txt','r');
fidH2=fopen('H2.txt','r');
fidOH=fopen('OH.txt','r');
fidH=fopen('H.txt','r');
fidO=fopen('O.txt','r');
fCO2=textscan(fidCO2,'%f');
fN2=textscan(fidN2,'%f');
fO2=textscan(fidO2,'%f');
fH2O=textscan(fidH2O,'%f');
fAr=textscan(fidAr,'%f');
fCO=textscan(fidCO,'%f');
fH2=textscan(fidH2,'%f');
fOH=textscan(fidOH,'%f');
fH=textscan(fidH,'%f');
fO=textscan(fidO,'%f');
CO2c=fCO2{1};
N2c=fN2{1};
O2c=fO2{1};
H2Oc=fH2O{1};
Arc=fAr{1};
COc=fCO{1};
H2c=fH2{1};
OHc=fOH{1};
Hc=fH{1};
Oc=fO{1};
fclose(fidCO2);
fclose(fidN2);
fclose(fidO2);
fclose(fidH2O);
fclose(fidAr);
fclose(fidCO);
fclose(fidH2);
fclose(fidOH);
fclose(fidH);
fclose(fidO);
cd ..

% --- Zero needed arrays for specific subroutines
Xu = zeros(nps,NM+1);
Yu = zeros(nps,NM+1);
cpmu = zeros(nps,NM+1);
cvmu = zeros(nps,NM+1);
hmu = zeros(nps,NM+1);
umu = zeros(nps,NM+1);
somu = zeros(nps,NM+1);
cpu = zeros(nps,NM+1);
cvu = zeros(nps,NM+1);
hu = zeros(nps,NM+1);
uu = zeros(nps,NM+1);
sou = zeros(nps,NM+1);
% --- zero arrays for output
EXRGYout = zeros(nps,40);                        % Output of heat release
% 1 - dEw/dt - exergy due to work [W]
% 2 - dEw/deg - exergy due to work [J/deg]
% 3 - dEHTc/dt - exergy due to convective heat transfer [W]
% 4 - dEHTc/deg - exergy due to convective heat transfer [J/deg]
% 5 - dEHTr/dt - exergy due to radiative heat transfer [W]
% 6 - dEHTr/deg - exergy due to radiative heat transfer [J/deg]
% 7 - dEHTfl/dt - exergy due to the heating of the fuel [W]
% 8 - dEHTfl/deg - exergy due to the heating of the fuel [J/deg]
% 9 - dEHT/dt - total exergy due to heat transfer [W]
% 10 - dEHT/deg - total exergy due to heat transfer [J/deg]
% 11 - dEin/dt - exergy due to the flow in [W]
% 12 - dEin/deg - exergy due to the flow in [J/deg]
% 13 - dEd/dt - exergy destruction [W]
% 14 - dEd/deg - exergy destruction [J/deg]
% 15 - dECV/dt - exergy change in the control volume [W]
% 16 - dECV/deg - exergy change in the control volume [J/deg]
% 17 - Exergy Heat Release [W]
% dEHR/dt = dECV/dt - dEH/dt + dEw/dt - dEin/dt + dEd/dt
% 18 - Exergy Heat Release [J/deg]
% 19 - dEHT/dt corrected - [W]. This is the corrected exergy heat transfer
% that balances the Conservation of Energy equation
% 20 - dEHT/deg corrected [J/deg]
% 21 - dEd/dt corrected - [W]. This is the corrected exergy destruction
% using the corrected heat transfer values that balances the Conservation
% of Energy equation
% 22 - dEd/deg corrected [J/deg]
% 23 - Corrected Exergy Heat Release [W]
% dEHR/dt,corrected = dECV/dt - dEH/dt,corrected + dEw/dt - dEin/dt + dEd/dt,corrected
% 24 - dEHR/deg corrected [J/deg]
% 25 - dEHR/dt [W] - theoretical exergy heat release
% 26 - dEHR/deg [J/deg] - theoretical exergy heat release
% 27 - EHR [J] - total theoretical exergy heat release
% Cumulative values
% 28 - Ew [J] - total work exergy
% 29 - EHTrc [J] - total heat transfer exergy due to convection and radiation
% 30 - EHTfl [J] - total heat transfer exergy to the fuel
% 31 - Ein [J] - total exergy flow in
% 32 - Ed [J] - total exergy destruction
% 33 - ECV [J] - total exergy change in CV
% 34 - EHR [J] - total exergy heat release 
% HR1out
% 1 - dWCV/dt [W]
% 4 - dQHTc/dt [W]
% 7 - dQHTr/dt [W]
% HR2out
% 9 - sigmaCV [W/K]
for j=IVC+1:EVO     % IVC to EVO
    % Calculate the work exergy [W]
    % dEw/dt = dWCV/dt - po*dV/dt
    EXRGYout(j,1) = HR1out(j,1) - po*pdata(j,7)*rps;            % [W]
    EXRGYout(j,2) = EXRGYout(j,1)/dps;                          % [J/deg]
    EXRGYout(j,28) = EXRGYout(j-1,28) + EXRGYout(j,1)*dt;       % [J]
    % Heat transfer exergy [W]
    % Convection
    % dEHTc/dt = dQHTc/dt*(1-To/Twall)
    EXRGYout(j,3) = HR1out(j,4)*(1 - To/Twall);
    EXRGYout(j,4) = EXRGYout(j,3)/dps;
    % Radiation
    % dEHTr/dt = dQHTr/dt*(1-To/Twall)
    EXRGYout(j,5) = HR1out(j,7)*(1 - To/Twall);
    EXRGYout(j,6) = EXRGYout(j,5)/dps;
    % Heat transfer to fuel - calculated in prior routine
    % dEHTfl/dt = dQHTfl/dt*(1-To/Tb)
    EXRGYout(j,7) = dEHTfldt(j);        % [W]
    EXRGYout(j,8) = EXRGYout(j,7)/dps;  % [J/deg]
    % Total heat transfer exergy
    % dEHT/dt
    EXRGYout(j,9) = EXRGYout(j,3) + EXRGYout(j,5) + EXRGYout(j,7);
    EXRGYout(j,10) = EXRGYout(j,9)/dps;
    EXRGYout(j,29) = EXRGYout(j-1,29) + (EXRGYout(j,3)+EXRGYout(j,5))*dt;
    EXRGYout(j,30) = EXRGYout(j-1,30) + EXRGYout(j,7)*dt;
    % Specific flow exergy in - calculated in prior routine
    % dEin/dt = dmfl/dt*efl
    EXRGYout(j,11) = dEindt(j);             % [W]
    EXRGYout(j,12) = EXRGYout(j,11)/dps;    % [J/deg]
    EXRGYout(j,31) = EXRGYout(j-1,31) + EXRGYout(j,11)*dt;      % [J]
    % Exergy destruction
    % dEd/dt = To*sigmadot
    EXRGYout(j,13) = To*HR2out(j,9); % [W]
    EXRGYout(j,14) = EXRGYout(j,13)/dps;    % [J/deg]
    EXRGYout(j,32) = EXRGYout(j-1,32) + EXRGYout(j,13)*dt;  % [J]
    % Exergy of the control volume
    % dECV/dt
    EXRGYout(j,15) = (ECV(j)-ECV(j-1))/dt;      % [W]
    EXRGYout(j,16) = EXRGYout(j,15)/dps;        % [J/deg]
    EXRGYout(j,33) = EXRGYout(j-1,33) + EXRGYout(j,15)*dt;  % [J]
    % This equation typically equals zero when the species are not changing
    EXRGYout(j,17) = EXRGYout(j,15) - EXRGYout(j,9) + EXRGYout(j,1) - EXRGYout(j,11) + EXRGYout(j,13);  % [W]
    % However, reviewing the Conservation of Energy equation, I do believe
    % the resulting value here is the exergetic heat release. 
    % dEHR/dt = dECV/dt - dEH/dt + dEw/dt - dEin/dt + dEd/dt
    EXRGYout(j,18) = EXRGYout(j,17)/dps;        % [J/deg]
    EXRGYout(j,34) = EXRGYout(j-1,34) + EXRGYout(j,17)*dt;      % [J]
    % What if we use the theoretical values of convective heat transfer (no
    % radiation) and the corresponding entropy generation?
    EXRGYout(j,19) = HR1out(j,27)*(1-To/Twall);         % Corrected exergy heat transfer [W]
    EXRGYout(j,20) = EXRGYout(j,19)/dps;            % Corrected exergy heat transfer [J/deg]
    % Corrected exergy destruction that comes from the corrected heat
    % transfer values in the Conservation of Entropy equation
    EXRGYout(j,21) = To*HR2out(j,22);                    % Corrected exergy destruction [W]
    EXRGYout(j,22) = EXRGYout(j,21)/dps;        % Corrected exergy destruction [J/deg]
    % Corrected exergy heat release
    % dEHR/dt,corrected = dECV/dt - dEHT/dt,corrected - dEHTfl/dt + dEw/dt
    % - dEin/dt + dEd/dt,corrected
    EXRGYout(j,23) = EXRGYout(j,15) - EXRGYout(j,19) - EXRGYout(j,7) + EXRGYout(j,1) - EXRGYout(j,11) + EXRGYout(j,21); % [W]
    EXRGYout(j,24) = EXRGYout(j,23)/dps;        % Corrected exergy heat release [J/deg]
end

% --- Theoretical exergy heat release
for j=IVC+1:EVO    % IVC to EVO
    % Find the internal energy of the species at the control volume
    % temperature
    % Unburned & burned - works for both. We just need the internal energy
    % of the different species at the temperture of the control volume
    mou_O2 = mu(j,1)/MO2;
    mou_N2 = mu(j,2)/MN2;
    mou_Ar = mu(j,3)/MAr;
    mou_CO2 = mu(j,4)/MCO2;
    mou_H2O = mu(j,5)/MH2O;
    mou_CO = mu(j,6)/MCO;
    mou_H2 = mu(j,7)/MH2;
    mou_OH = mu(j,8)/MOH;
    mou_H = mu(j,9)/MH;
    mou_O = mu(j,10)/MO;
    [~, ~, ~, ~, ~, ~, ~, ~, somu, ~, ~, ~, uu, ~, ~, ~] = ...
        unburnprop(j, Runiv, pref, pdata(j,13), TCV(j), MM, ...
        mou_O2, mou_N2, mou_Ar, mou_CO2, mou_H2O, mou_CO, mou_H2, mou_OH, mou_H, mou_O, ...
        O2c, N2c, Arc, CO2c, H2Oc, COc, H2c, OHc, Hc, Oc, ...
        mu, Xu, Yu, cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou);
    % Actual entropy of the species
    s(1) = (somu(j,1) - Runiv*log(XCV(j,1)*pdata(j,13)/pref))/MO2;
    s(2) = (somu(j,2) - Runiv*log(XCV(j,2)*pdata(j,13)/pref))/MN2; 
    s(3) = (somu(j,3) - Runiv*log(XCV(j,3)*pdata(j,13)/pref))/MAr; 
    s(4) = (somu(j,4) - Runiv*log(XCV(j,4)*pdata(j,13)/pref))/MCO2;
    s(5) = (somu(j,5) - Runiv*log(XCV(j,5)*pdata(j,13)/pref))/MH2O; 
    s(6) = (somu(j,6) - Runiv*log(XCV(j,6)*pdata(j,13)/pref))/MCO; 
    s(7) = (somu(j,7) - Runiv*log(XCV(j,7)*pdata(j,13)/pref))/MH2; 
    s(8) = (somu(j,8) - Runiv*log(XCV(j,8)*pdata(j,13)/pref))/MOH; 
    s(9) = (somu(j,9) - Runiv*log(XCV(j,9)*pdata(j,13)/pref))/MH; 
    s(10) = (somu(j,10) - Runiv*log(XCV(j,10)*pdata(j,13)/pref))/MO; 
    u(1) = uu(j,1);
    u(2) = uu(j,2);
    u(3) = uu(j,3);
    u(4) = uu(j,4);
    u(5) = uu(j,5);
    u(6) = uu(j,6);
    u(7) = uu(j,7);
    u(8) = uu(j,8);
    u(9) = uu(j,9);
    u(10) = uu(j,10);
    % Added gas
    if (mfa(j) > 0)
        mofa = mfa(j)/Mfa;
        [~, ~, ~, ~, ~, ~, ~, somfa, ~, ~, ~, u(NM+1), ~] = faprop(Runiv, pref, pdata(j,13), TCV(j), Mfa, mofa, fuelchem);
        s(NM+1) = (somfa - Runiv*log(XCV(j,NM+1)*pdata(j,13)/pref))/Mfa; 
    else
        u(NM+1) = 0;
        s(NM+1) = 0;
    end
    % Liquid fuel
    if (mfl(j) > 0) 
        mofl = mfl(j)/Mfl;
        [~, ~, ~, ~, ~, ~, ~, somfl, ~, ~, ~, u(NM+2), ~] = flprop(Runiv, pref, pdata(j,13), TCV(j), Mfl, mofl, fuelchem);
        s(NM+2) = (somfl - Runiv*log(XCV(j,NM+2)*pdata(j,13)/pref))/Mfl; 
    else
        u(NM+2) = 0;
        s(NM+2) = 0;
    end
    % Let's compute the change in mass fractions
    % index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
    % NM+1-added gas, NM+2-liquid direct injected fuel
    btm = YCV(j,1)/MO2 + YCV(j,2)/MN2 + YCV(j,3)/MAr + YCV(j,4)/MCO2 + YCV(j,5)/MH2O + ...
        YCV(j,6)/MCO + YCV(j,7)/MH2 + YCV(j,8)/MOH + YCV(j,9)/MH + YCV(j,10)/MO + ... 
        YCV(j,NM+1)/Mfa + YCV(j,NM+2)/Mfl;
    % Calculate dYdt
    for k=1:NM+2 
        dYdt(k) = (YCV(j,k)-YCV(j-1,k))/dt;
    end
    top = dYdt(1)/MO2 + dYdt(2)/MN2 + dYdt(3)/MAr + dYdt(4)/MCO2 + dYdt(5)/MH2O + ...
        dYdt(6)/MCO + dYdt(7)/MH2 + dYdt(8)/MOH + dYdt(9)/MH + dYdt(10)/MO + ...
        dYdt(NM+1)/Mfa + dYdt(NM+2)/Mfl;
    % We need to find dR/dt. This will require dX/dt and dMmix/dt
    dXdt(1) = (dYdt(1)/MO2*btm - YCV(j,1)/MO2*top)/(btm*btm);
    dXdt(2) = (dYdt(2)/MN2*btm - YCV(j,2)/MN2*top)/(btm*btm);
    dXdt(3) = (dYdt(3)/MAr*btm - YCV(j,3)/MAr*top)/(btm*btm);
    dXdt(4) = (dYdt(4)/MCO2*btm - YCV(j,4)/MCO2*top)/(btm*btm);
    dXdt(5) = (dYdt(5)/MH2O*btm - YCV(j,5)/MH2O*top)/(btm*btm);
    dXdt(6) = (dYdt(6)/MCO*btm - YCV(j,6)/MCO*top)/(btm*btm);
    dXdt(7) = (dYdt(7)/MH2*btm - YCV(j,7)/MH2*top)/(btm*btm);
    dXdt(8) = (dYdt(8)/MOH*btm - YCV(j,8)/MOH*top)/(btm*btm);
    dXdt(9) = (dYdt(9)/MH*btm - YCV(j,9)/MH*top)/(btm*btm);
    dXdt(10) = (dYdt(10)/MO*btm - YCV(j,10)/MO*top)/(btm*btm);
    dXdt(NM+1) = (dYdt(NM+1)/Mfa*btm - YCV(j,NM+1)/Mfa*top)/(btm*btm);
    dXdt(NM+2) = (dYdt(NM+2)/Mfl*btm - YCV(j,NM+2)/Mfl*top)/(btm*btm);
    % dMmix/dt
    dMmixdt = MO2*dXdt(1) + MN2*dXdt(2) + MAr*dXdt(3) + MCO2*dXdt(4) + MH2O*dXdt(5) + ...
        MCO*dXdt(6) + MH2*dXdt(7) + MOH*dXdt(8) + MH*dXdt(9) + MO*dXdt(10) + ...
        Mfa*dXdt(NM+1) + Mfl*dXdt(NM+2);
    % dR/dt
    Mmix = XCV(j,1)*MO2 + XCV(j,2)*MN2 + XCV(j,3)*MAr + XCV(j,4)*MCO2 + XCV(j,5)*MH2O + ...
        XCV(j,6)*MCO + XCV(j,7)*MH2 + XCV(j,8)*MOH + XCV(j,9)*MH + XCV(j,10)*MO + ...
        XCV(j,NM+1)*Mfa + XCV(j,NM+2)*Mfl;
    dRdt = -Runiv/(Mmix*Mmix)*dMmixdt;
    % dEHR/dt (from theory)
    for k=1:NM+2
        if (dYdt(k) ~= 0)
            EXRGYout(j,25) = EXRGYout(j,25) - mCV(j)*(u(k) + po*TCV(j)/pdata(j,13)*dRdt/dYdt(k) - To*s(k))*dYdt(k);
        end
    end
    EXRGYout(j,26) = EXRGYout(j,25)/dps;
    EXRGYout(j,27) = EXRGYout(j-1,27) + EXRGYout(j,25)*dt;
end