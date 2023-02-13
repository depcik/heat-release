function HR2out = secondlawHR(pdata, TCV, SCV, XCV, YCV, cpCV, RCV, dSHTfldt, dSindt, HR1out, mu, mfa, mfl, mCV, inp, calc, fuel, fuelchem)
% --- Routine specific variables
NM = inp(13);                               % Total number of chemical species
nps = calc(2);                              % Total number of data points
IVC = calc(3);                             % Index on IVC [-]
EVO = calc(4);                              % Index on EVO [-]
dt = calc(12);                              % Time-step of data [s] 
Twall = inp(82);                            % Wall temperature [K]
dps = calc(11);                             % Degrees per second [deg/s] = [rev/min]*[360 deg/rev]*[min/60 s]
MO2 = 2*inp(102);                           % Molecular mass of O2 [kg/mol]
MN2 = 2*inp(103);                           % Molecular mass of N2 [kg/mol]
MAr = inp(105);                             % Molecular mass of Ar [kg/mol]
MCO2 = inp(101) + MO2;                      % Molecular mass of CO2 [kg/mol]
MH2O = 2*inp(104) + inp(102);               % Molecular mass of H2O [kg/mol]
MCO = inp(101) + inp(102);                  % Molecular mass of CO [kg/mol]
MOH = inp(102) + inp(104);                  % Molecular mass of OH [kg/mol]
MH2 = 2*inp(104);                           % Molecular mass of H2 [kg/mol]
MH = inp(104);                              % Molecular mass of H [kg/mol]
MO = inp(102);                              % Molecular mass of O [kg/mol]
% index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O, NM+1-Mixture
MM = [MO2 MN2 MAr MCO2 MH2O MCO MH2 MOH MH MO]; % Molecular mass array [kg/mol]
Mfa = fuel(24);                             % Molecular mass of added gas [kg/mol]
Mfl = fuel(5);                              % Molecular mass of liquid fuel [kg/mol]
Runiv = inp(100);                           % Universal gas constant [J/molK]
pref = inp(83);                             % CHEMKIN reference pressure [Pa]

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

% --- zero arrays for output
HR2out = zeros(nps,30);                        % Output of heat release
% 1 - dSCV/dt [W/K]
% 2 - dSCV/deg [J/(deg*K)]
% 3 - dSHT/dt [W/K]
% 4 - dSHT/deg [J/(deg*K)]
% 5 - dSHTfl/dt [W/K]
% 6 - dSHTfl/deg [J/(deg*K)]
% 7 - dSin/dt [W/K]
% 8 - dSin/deg [J/(deg*K)]
% 9 - sigmaCV [W/K]
% 10 - sigmaCVdeg [J/degK]
% 11 - sigma,total [J/K] - total entropy generated
% 12 - dSHR/dt [W/K]
% 13 - dSHR/deg [J/degK]
% 14 - SHR [J/K]
% 15 - Entropy generated through conv. & radiative heat transfer [J/K]
% 16 - Entropy generated through the flow entering [J/K]
% 17 - Entropy generated through the change in control volume [J/K]
% 18 - CV entropy based on temperature change [W/K]
% 19 - Entropy generated through CV temperature change [J/K]
% 20 - CV entropy based on pressure change [W/K]
% 21 - Entropy generated through CV pressure change [J/K]
% 22 - sigmaCV, theory [W/K] - using the theoretical values of heat
% transfer that causes the Conservation of Energy to balance
% 23 - sigmaCVdeg, theory [J/degK]
% 24 - sigma,theory total [J/K] - total entropy generated using the
% balanced version of convective heat transfer
% 25 - Entropy generated through heat transfer to liquid fuel [J/K] 

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

% HR1out(j,4) = dQHTc/dt
% HR1out(j,7) = dQHTr/dt
% Compute the entropy generation at each stage
for j=IVC+1:EVO     % IVC to EVO
    % Change in control volume entropy
    % dSCV/dt
    HR2out(j,1) = (SCV(j)-SCV(j-1))/dt;       % [J/s] = [W]
    HR2out(j,2) = HR2out(j,1)/dps;
    % Heat transfer to and from the wall
    % (dQHTc/dt + dQHTr/dt)/Twall
    HR2out(j,3) = (HR1out(j,4)+HR1out(j,7))/Twall;        % [W]/[K]
    HR2out(j,4) = HR2out(j,3)/dps;
    % We already have the entropy changes due to heat transfer to the
    % liquid fuel and the entropy coming in from the liquid fuel being
    % added
    % (dQHTfl/dt)/Tb
    % dSin/dt
    HR2out(j,5) = dSHTfldt(j);
    HR2out(j,6) = HR2out(j,5)/dps;
    HR2out(j,7) = dSindt(j);
    HR2out(j,8) = HR2out(j,7)/dps;
    % Entropy generation - using the calculated values of heat transfer
    % Sigma,dot
    HR2out(j,9) = HR2out(j,1) - HR2out(j,3) - HR2out(j,5) - HR2out(j,7);    % [W/K]
    HR2out(j,10) = HR2out(j,9)/dps;
    % Total entropy generated
    HR2out(j,11) = HR2out(j-1,11) + HR2out(j,9)*dt;       % [J/K]
    % Entropy generated through conv. & radiative heat transfer
    HR2out(j,15) = HR2out(j-1,15) - HR2out(j,3)*dt;
    % Entropy generated through heat transfer to the liquid fuel
    HR2out(j,25) = HR2out(j-1,25) - HR2out(j,5)*dt;
    % Entropy generated through the liquid fuel entering
    HR2out(j,16) = HR2out(j-1,16) - HR2out(j,7)*dt;
    % Entropy generated through the change in the control volume
    HR2out(j,17) = HR2out(j-1,17) + HR2out(j,1)*dt;
    % As stated prior, the heat transfer used is a correlation that does
    % not cause the Conservation of Energy equation to balance at each
    % crank angle. We computed what that heat transfer should be to
    % balance: HR1out(j,27) 
    % Let us compute entropy generation using this result. Keep in mind
    % that we are ignoring radiation (it is lumped in the convective
    % result)
    HR2out(j,22) = HR2out(j,1) - HR1out(j,27)/Twall - HR2out(j,5) - HR2out(j,7);
    HR2out(j,23) = HR2out(j,22)/dps;
    HR2out(j,24) = HR2out(j-1,24) + HR2out(j,22)*dt;    % [J/K]
    % Can we determine an entropy heat release
    % No, entropy generation includes this term - it is the catch all in
    % the balance
end

% Theory of entropy change due to heat release
for j=IVC+1:EVO    % IVC to EVO
    % Find the entropy of the species at the control volume temperature
    % Unburned & burned - works for both. We just need the entropy
    % of the different species at the temperture of the control volume
    % index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
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
    [~, ~, ~, ~, ~, ~, ~, ~, somu, ~, ~, ~, ~, ~, ~, ~] = ...
        unburnprop(j, Runiv, pref, pdata(j,13), TCV(j), MM, ...
        mou_O2, mou_N2, mou_Ar, mou_CO2, mou_H2O, mou_CO, mou_H2, mou_OH, mou_H, mou_O, ...
        O2c, N2c, Arc, CO2c, H2Oc, COc, H2c, OHc, Hc, Oc, ...
        mu, Xu, Yu, cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou);
    for kk=1:NM
        if (XCV(j,kk) ~= 0) 
            s(kk) = (somu(j,kk) - Runiv*log(XCV(j,kk)*pdata(j,13)/pref))/MM(kk);
        else
            s(kk) = 0;
        end
    end
    % Added gas
    if (mfa(j) > 0)
        mofa = mfa(j)/Mfa;
        [~, ~, ~, ~, ~, ~, ~, somfa, ~, ~, ~, ~, ~] = faprop(Runiv, pref, pdata(j,13), TCV(j), Mfa, mofa, fuelchem);
        s(NM+1) = (somfa - Runiv*log(XCV(j,NM+1)*pdata(j,13)/pref))/Mfa; 
    else
        s(NM+1) = 0;
    end
    % Liquid fuel
    if (mfl(j) > 0)
        mofl = mfl(j)/Mfl;
        [~, ~, ~, ~, ~, ~, ~, somfl, ~, ~, ~, ~, ~] = flprop(Runiv, pref, pdata(j,13), TCV(j), Mfl, mofl, fuelchem);
        s(NM+2) = (somfl - Runiv*log(XCV(j,NM+2)*pdata(j,13)/pref))/Mfl; 
    else
        s(NM+2) = 0;
    end
    % Let's compute the change in mass fractions
    % index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
    % NM+1-added gas, NM+2-liquid direct injected fuel
    HR2out(j,12) = 0;
    for k=1:NM+2 
        dYdt = (YCV(j,k)-YCV(j-1,k))/dt;
        HR2out(j,12) = HR2out(j,12) + mCV(j)*s(k)*dYdt;
    end
    HR2out(j,13) = HR2out(j,12)/dps;
    HR2out(j,14) = HR2out(j-1,14) + HR2out(j,12)*dt;
    % Other terms in the control volume
    HR2out(j,18) = mCV(j)*cpCV(j)/TCV(j)*((TCV(j)-TCV(j-1))/dt);
    HR2out(j,19) = HR2out(j-1,19) + HR2out(j,18)*dt;
    vCV = RCV(j)*TCV(j)/pdata(j,13);    % specific volume in control volume
    HR2out(j,20) = -mCV(j)*vCV/TCV(j)*((pdata(j,13)-pdata(j-1,13))/dt);
    HR2out(j,21) = HR2out(j-1,21) + HR2out(j,20)*dt;
end
% Unlike our first law heat release where as we add fuel, they are bringing
% energy with them. Thus, we see an increase in the heat release. Fuels
% have a large negative value of internal energy. Since we are subtracting
% off this large negative, we see a corresponding energy change. Think of
% it like changing the datum states as we bring in species with their own
% heats of formation. However, the entropy datum state is the same.
% In addition, what "is" entropy heat release? It is the entropy caused by the
% conversion of the fuel to products. When we add liquid fuel into the
% cylinder, we are not creating any products. Instead, we are adding fuel
% which is the opposite. 