function [mCV, YCV, XCV, RCV, TCV, cpCV, cvCV, UCV, HCV, SCV, ECV] = ...
    cvprop(i, inp, Runiv, pref, pCV, VCV, YCV, XCV, MM, ...
    Mfa, Mfl, mu, mfa, mfl, mb, ...
    O2c, N2c, Arc, CO2c, H2Oc, COc, H2c, OHc, Hc, Oc, fuel, fuelchem)
% Set the molecular masses
% index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O, NM+1-Mixture
NM = length(MM);
MO2 = MM(1);
MN2 = MM(2);
MAr = MM(3);
MCO2 = MM(4);
MH2O = MM(5);
MCO = MM(6);
MH2 = MM(7);
MOH = MM(8);
MH = MM(9);
MO = MM(10);
% Routine values
To = inp(84);                   % Exergy reference temperature [K]
po = inp(85);                   % Exergy reference pressure [Pa]
% Molar chemical exergy of the species [J/mol]
% http://web.mit.edu/2.813/www/readings/APPENDIX.pdf
% https://www.sciencedirect.com/science/article/pii/S0360544206000569
emchO2 = 3.97*1000;
emchN2 = 0.72*1000;
emchAr = 12*1000;
emchCO2 = 19.48*1000;
emchH2O = 9.5*1000;
emchCO = 274.71*1000;   
emchH2 = 236.09*1000;
emchOH = 153940;    % http://publicationslist.org.s3.amazonaws.com/data/miguel.a.lozano/ref-93/ar15_a15(cp0027).pdf
emchH = 331.3*1000;
emchO = 233.7*1000;
emchfa = fuel(30);    
emchfl = fuel(11);
% cv, cp, u, and h for ideal gases are only a function of temperature
% s for ideal gases are a function of temperature and partial pressure

%--- Zero Arrays --%
Xu = zeros(i,NM+1);
Yu = zeros(i,NM+1);
% At TCV and PCV
cpmu = zeros(i,NM+1);
cvmu = zeros(i,NM+1);
hmu = zeros(i,NM+1);
umu = zeros(i,NM+1);
somu = zeros(i,NM+1);
cpu = zeros(i,NM+1);
cvu = zeros(i,NM+1);
hu = zeros(i,NM+1);
uu = zeros(i,NM+1);
sou = zeros(i,NM+1);
% At To and po
cpmuo = zeros(i,NM+1);
cvmuo = zeros(i,NM+1);
hmuo = zeros(i,NM+1);
umuo = zeros(i,NM+1);
somuo = zeros(i,NM+1);
cpuo = zeros(i,NM+1);
cvuo = zeros(i,NM+1);
huo = zeros(i,NM+1);
uuo = zeros(i,NM+1);
souo = zeros(i,NM+1);

% We do not need burned gas arrays since they will be evaluated at the same
% control volume temperature

% Calculate the mass of the total control volume
mCV = mu(i,NM+1) + mfa + mfl + mb(i,NM+1);    % Total mass [kg] in control volume
% Calculate the mass fractions in the control volume
% The burned mass values will also be included in this once we add them
% into the routine. Everything is at "one" temperature when solving for the
% control volume
for j=1:NM
    YCV(i,j) = (mu(i,j)+mb(i,j))/mCV;
end
YCV(i,NM+1) = mfa/mCV;
YCV(i,NM+2) = mfl/mCV;
% Mole fractions within the control volume
% 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
btm = YCV(i,1)/MO2 + YCV(i,2)/MN2 + YCV(i,3)/MAr + YCV(i,4)/MCO2 + YCV(i,5)/MH2O + ...
    YCV(i,6)/MCO + YCV(i,7)/MH2 + YCV(i,8)/MOH + YCV(i,9)/MH + YCV(i,10)/MO + ...
    YCV(i,NM+1)/Mfa + YCV(i,NM+2)/Mfl;
MCV = 1/btm;   % Mixture molecular mass [kg/mol]
XCV(i,1) = YCV(i,1)/MO2/btm;
XCV(i,2) = YCV(i,2)/MN2/btm;
XCV(i,3) = YCV(i,3)/MAr/btm;
XCV(i,4) = YCV(i,4)/MCO2/btm;
XCV(i,5) = YCV(i,5)/MH2O/btm;
XCV(i,6) = YCV(i,6)/MCO/btm;
XCV(i,7) = YCV(i,7)/MH2/btm;
XCV(i,8) = YCV(i,8)/MOH/btm;
XCV(i,9) = YCV(i,9)/MH/btm;
XCV(i,10) = YCV(i,10)/MO/btm;
XCV(i,NM+1) = YCV(i,NM+1)/Mfa/btm;
XCV(i,NM+2) = YCV(i,NM+2)/Mfl/btm;
% Gas constant
RCV = Runiv/MCV;            % J/kgK
% Temperature of the control volume
TCV = pCV*VCV/(mCV*RCV);  % [K]
% Find the properties of the control volume
% Get the unburned gas internal energy at the temperature of the control volume
% 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
mou_O2 = mu(i,1)/MO2;           % [mol] = [kg] * [mol/kg]
mou_N2 = mu(i,2)/MN2;
mou_Ar = mu(i,3)/MAr;
mou_CO2 = mu(i,4)/MCO2;
mou_H2O = mu(i,5)/MH2O;
mou_CO = mu(i,6)/MCO;
mou_H2 = mu(i,7)/MH2;
mou_OH = mu(i,8)/MOH;
mou_H = mu(i,9)/MH;
mou_O = mu(i,10)/MO;
[~, ~, ~, ~, ~, ~, ~, ~, somu, cpu, cvu, hu, uu, ~, ~, ~] = ...
    unburnprop(i, Runiv, pref, pCV, TCV, MM, ...
        mou_O2, mou_N2, mou_Ar, mou_CO2, mou_H2O, mou_CO, mou_H2, mou_OH, mou_H, mou_O, ...
        O2c, N2c, Arc, CO2c, H2Oc, COc, H2c, OHc, Hc, Oc, ...
        mu, Xu, Yu, cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou);
% so entropy values are returned (standard-state); thus, we will need to find the partial
% Reference condition values for exergy: po & To
% uuo - returns the reference values of internal energy
% somuo - returns the standard state values of entropy
[~, ~, ~, ~, ~, ~, ~, ~, somuo, ~, ~, ~, uuo, ~, ~, ~] = ...
    unburnprop(i, Runiv, pref, po, To, MM, ...
        mou_O2, mou_N2, mou_Ar, mou_CO2, mou_H2O, mou_CO, mou_H2, mou_OH, mou_H, mou_O, ...
        O2c, N2c, Arc, CO2c, H2Oc, COc, H2c, OHc, Hc, Oc, ...
        mu, Xu, Yu, cpmuo, cvmuo, hmuo, umuo, somuo, cpuo, cvuo, huo, uuo, souo);
% pressure values in [Pa]
pO2 = pCV*XCV(i,1);
pN2 = pCV*XCV(i,2);
pAr = pCV*XCV(i,3);
pCO2 = pCV*XCV(i,4);
pH2O = pCV*XCV(i,5);
pCO = pCV*XCV(i,6);
pH2 = pCV*XCV(i,7);
pOH = pCV*XCV(i,8);
pH = pCV*XCV(i,9);
pO = pCV*XCV(i,10);
pmfa = pCV*XCV(i,NM+1);
pmfl = pCV*XCV(i,NM+2);
% Find the true molar entropy values
% True molar entropy values for exergy at dead state requires a pressure
% correction
if (pO2 ~= 0)
    smCV(1) = somu(i,1) - Runiv*log(pO2/pref);   % [J/molK]
    smoCV(1) = somuo(i,1) - Runiv*log(XCV(i,1)*po/pref);
else
    smCV(1) = 0;
    smoCV(1) = 0;
end
if (pN2 ~= 0)
    smCV(2) = somu(i,2) - Runiv*log(pN2/pref);   % [J/molK]
    smoCV(2) = somuo(i,2) - Runiv*log(XCV(i,2)*po/pref);
else
    smCV(2) = 0;
    smoCV(2) = 0;
end
if (pAr ~= 0)
    smCV(3) = somu(i,3) - Runiv*log(pAr/pref);
    smoCV(3) = somuo(i,3) - Runiv*log(XCV(i,3)*po/pref);
else
    smCV(3) = 0;
    smoCV(3) = 0;
end
if (pCO2 ~= 0)
    smCV(4) = somu(i,4) - Runiv*log(pCO2/pref);
    smoCV(4) = somuo(i,4) - Runiv*log(XCV(i,4)*po/pref);
else
    smCV(4) = 0;
    smoCV(4) = 0;
end
if (pH2O ~= 0)
    smCV(5) = somu(i,5) - Runiv*log(pH2O/pref);
    smoCV(5) = somuo(i,5) - Runiv*log(XCV(i,5)*po/pref);
else
    smCV(5) = 0;
    smoCV(5) = 0;
end
if (pCO ~= 0)
    smCV(6) = somu(i,6) - Runiv*log(pCO/pref);
    smoCV(6) = somuo(i,6) - Runiv*log(XCV(i,6)*po/pref);
else
    smCV(6) = 0;
    smoCV(6) = 0;
end
if (pH2 ~= 0)
    smCV(7) = somu(i,7) - Runiv*log(pH2/pref);
    smoCV(7) = somuo(i,7) - Runiv*log(XCV(i,7)*po/pref);
else
    smCV(7) = 0;
    smoCV(7) = 0;
end
if (pOH ~= 0)
    smCV(8) = somu(i,8) - Runiv*log(pOH/pref);
    smoCV(8) = somuo(i,8) - Runiv*log(XCV(i,8)*po/pref);
else
    smCV(8) = 0;
    smoCV(8) = 0;
end
if (pH ~= 0)
    smCV(9) = somu(i,9) - Runiv*log(pH/pref);
    smoCV(9) = somuo(i,9) - Runiv*log(XCV(i,9)*po/pref);
else
    smCV(9) = 0;
    smoCV(9) = 0;
end
if (pO ~= 0)
    smCV(10) = somu(i,10) - Runiv*log(pO/pref); 
    smoCV(10) = somuo(i,10) - Runiv*log(XCV(i,10)*po/pref);
else
    smCV(10) = 0;
    smoCV(10) = 0;
end
% Get the added fuel properties
mofa = mfa/Mfa;
if (mfa ~= 0)
    [~, ~, ~, ~, ~, ~, ~, somfa, cpfa, cvfa, hfa, ufa, ~] = faprop(Runiv, pref, pCV, TCV, Mfa, mofa, fuelchem);
    smCV(NM+1) = somfa - Runiv*log(pmfa/pref);     % [J/molK]
    % Reference conditions for exergy
    [~, ~, ~, ~, ~, ~, ~, somfao, ~, ~, ~, ufao, ~] = faprop(Runiv, pref, po, To, Mfa, mofa, fuelchem);
    smoCV(NM+1) = somfao - Runiv*log(XCV(i,NM+1)*po/pref);  
else
    cpfa = 0; cvfa = 0; hfa = 0; ufa = 0; ufao = 0;
    smCV(NM+1) = 0;
    smoCV(NM+1) = 0;
end
% Get the liquid fuel properties
mofl = mfl/Mfl;
if (mfl ~= 0)
    [~, ~, ~, ~, ~, ~, ~, somfl, cpfl, cvfl, hfl, ufl, ~] = flprop(Runiv, pref, pCV, TCV, Mfl, mofl, fuelchem);
    smCV(NM+2) = somfl - Runiv*log(pmfl/pref);     % [J/molK]
    % Reference conditions for exergy
    [~, ~, ~, ~, ~, ~, ~, somflo, ~, ~, ~, uflo, ~] = flprop(Runiv, pref, po, To, Mfl, mofl, fuelchem);
    smoCV(NM+2) = somflo - Runiv*log(XCV(i,NM+2)*po/pref);
else
    cpfl = 0; cvfl = 0; hfl = 0; ufl = 0; uflo = 0; 
    smCV(NM+2) = 0;
    smoCV(NM+2) = 0;
end

%% Get the burned fuel properties
% This is NOT needed since we would evaluate this at the same pressure and
% temperature of the unburned properties

%% Compute the mixture properties
cpCV = 0;
cvCV = 0;
uCV = 0;
hCV = 0;
for j=1:NM
    cpCV = cpCV + YCV(i,j)*cpu(i,j);
    cvCV = cvCV + YCV(i,j)*cvu(i,j);
    uCV = uCV + YCV(i,j)*uu(i,j);
    hCV = hCV + YCV(i,j)*hu(i,j);
end
cpCV = cpCV + YCV(i,NM+1)*cpfa + YCV(i,NM+2)*cpfl;
cvCV = cvCV + YCV(i,NM+1)*cvfa + YCV(i,NM+2)*cvfl;
uCV = uCV + YCV(i,NM+1)*ufa + YCV(i,NM+2)*ufl;
hCV = hCV + YCV(i,NM+1)*hfa + YCV(i,NM+2)*hfl;
% Entropy
smCVt = 0;
for j=1:NM+2
    smCVt = smCVt + XCV(i,j)*smCV(j);
end
sCV = smCVt/MCV;     % [J/kgK] = [J/molK] * [mol/kg]

%% Exergy
% Standard state given the species in the control volume
uo = YCV(i,1)*uuo(i,1) + YCV(i,2)*uuo(i,2) + YCV(i,3)*uuo(i,3) + YCV(i,4)*uuo(i,4) + YCV(i,5)*uuo(i,5) + ...
    YCV(i,6)*uuo(i,6) + YCV(i,7)*uuo(i,7) + YCV(i,8)*uuo(i,8) + YCV(i,9)*uuo(i,9) + YCV(i,10)*uuo(i,10) + ...
    YCV(i,NM+1)*ufao + YCV(i,NM+2)*uflo;
% Actual state given the species in the control volume at the temperature
u1 = YCV(i,1)*uu(i,1) + YCV(i,2)*uu(i,2) + YCV(i,3)*uu(i,3) + YCV(i,4)*uu(i,4) + YCV(i,5)*uu(i,5) + ...
    YCV(i,6)*uu(i,6) + YCV(i,7)*uu(i,7) + YCV(i,8)*uu(i,8) + YCV(i,9)*uu(i,9) + YCV(i,10)*uu(i,10) + ...
    YCV(i,NM+1)*ufa + YCV(i,NM+2)*ufl;
% Standard state specific volume
vo = RCV*To/po;
% Current specific volume
v1 = VCV/mCV;
% Find the standard state mixture
smo = 0;
for j=1:NM+2
    smo = smo + XCV(i,j)*smoCV(j);
end
% Conversion to mass-based units
so = smo/MCV;
% The actual entropy with all the species at their partial pressures
sm1 = XCV(i,1)*smCV(1) + XCV(i,2)*smCV(2) + XCV(i,3)*smCV(3) + XCV(i,4)*smCV(4) + XCV(i,5)*smCV(5) + ...
    XCV(i,6)*smCV(6) + XCV(i,7)*smCV(7) + XCV(i,8)*smCV(8) + XCV(i,9)*smCV(9) + XCV(i,10)*smCV(10) + ...
    XCV(i,NM+1)*smCV(NM+1) + XCV(i,NM+2)*smCV(NM+2);
% Conversion to mass based units
s1 = sm1/MCV;
% Chemical exergy [J/mol]
% Standard state values
emch = XCV(i,1)*emchO2 + XCV(i,2)*emchN2 + XCV(i,3)*emchAr + XCV(i,4)*emchCO2 + XCV(i,5)*emchH2O + ...
    XCV(i,6)*emchCO + XCV(i,7)*emchH2 + XCV(i,8)*emchOH + XCV(i,9)*emchH + XCV(i,10)*emchO + ...
    XCV(i,NM+1)*emchfa + XCV(i,NM+2)*emchfl;
% Conversion to mass based units
ech = emch/MCV; % [J/mol] * [mol/kg]
eCV = (u1-uo) + po*(v1-vo) - To*(s1-so) + ech;
% Extensive properties
UCV = mCV*uCV;
HCV = mCV*hCV;
SCV = mCV*sCV;
ECV = mCV*eCV;
