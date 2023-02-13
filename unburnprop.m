% We only return the so values of entropy since we need to know partial
% pressures to get the true entropy
function [mu, Xu, Yu, Mu, cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou, Ru, Vu] = ...
    unburnprop(i, Runiv, pref, pu, Tu, MM, ...
    mou_O2, mou_N2, mou_Ar, mou_CO2, mou_H2O, mou_CO, mou_H2, mou_OH, mou_H, mou_O, ...
    O2c, N2c, Arc, CO2c, H2Oc, COc, H2c, OHc, Hc, Oc, ...
    mu, Xu, Yu, cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou)
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
% Find unburned masses
mu(i,1) = mou_O2*MO2;           
mu(i,2) = mou_N2*MN2;
mu(i,3) = mou_Ar*MAr;
mu(i,4) = mou_CO2*MCO2;
mu(i,5) = mou_H2O*MH2O;
mu(i,6) = mou_CO*MCO;
mu(i,7) = mou_H2*MH2;
mu(i,8) = mou_OH*MOH;
mu(i,9) = mou_H*MH;
mu(i,10) = mou_O*MO;
%% --- Determine the mole fractions of the unburned gas
% In addition, we have Amagat's model version of the ideal gas law that 
% will tell us the respective volumes at IVC. To do this, we need to figure
% out the respective gas constants of each zone.
% unburned
btm = mou_O2+mou_N2+mou_Ar+mou_CO2+mou_H2O+mou_CO+mou_H2+mou_OH+mou_H+mou_O;
Xu(i,1) = mou_O2/btm;
Xu(i,2) = mou_N2/btm;
Xu(i,3) = mou_Ar/btm;
Xu(i,4) = mou_CO2/btm;
Xu(i,5) = mou_H2O/btm;
Xu(i,6) = mou_CO/btm;
Xu(i,7) = mou_H2/btm;
Xu(i,8) = mou_OH/btm;
Xu(i,9) = mou_H/btm;
Xu(i,10) = mou_O/btm;
Xu(i,NM+1) = 0;
% Mixture molecular mass
Mu = Xu(i,1)*MO2 + Xu(i,2)*MN2 + Xu(i,3)*MAr + Xu(i,4)*MCO2 + Xu(i,5)*MH2O + ...
    Xu(i,6)*MCO + Xu(i,7)*MH2 + Xu(i,8)*MOH + Xu(i,9)*MH + Xu(i,10)*MO; % Molecular mass [kg/mol]
% Determine the mass fractions
Yu(i,1) = Xu(i,1)*MO2/Mu;
Yu(i,2) = Xu(i,2)*MN2/Mu;
Yu(i,3) = Xu(i,3)*MAr/Mu;
Yu(i,4) = Xu(i,4)*MCO2/Mu;
Yu(i,5) = Xu(i,5)*MH2O/Mu;
Yu(i,6) = Xu(i,6)*MCO/Mu;
Yu(i,7) = Xu(i,7)*MH2/Mu;
Yu(i,8) = Xu(i,8)*MOH/Mu;
Yu(i,9) = Xu(i,9)*MH/Mu;
Yu(i,10) = Xu(i,10)*MO/Mu;
Yu(i,NM+1) = 0;
%% --- Determine the constant pressure, internal energy, enthalpy, and entropy of the unburned components
% Let us fill our arrays of constant pressure specific heat, internal
% energy, enthalpy, and entropy of the mixtures
% CHEMKIN curve-fits: first seven coefficients are those for 1000-5000K
% Second seven coefficients are those for 200-1000K
% Oxygen at 298.15 K: s(1bar) = 205.152 J/molK; cp ~ 29.4 J/molK; h = 0
% Nitrogen at 298.15 K: s(1bar) = 191.609 J/molK; cp ~ 29.1 J/molK; h = 0
% Argon at 298.15 K: s(1bar) = 151.84 J/molK; cp ~ 20.78 J/molK; h = 0
% Carbon Dioxide at 298.15: s(1bar) = 213.785 J/molK; cp ~ 37.1 J/molK; h = -393.51 kJ/mol
% Water at 298.15: s(1bar) = 188.835 J/molK; cp ~ 33.5 J/molK; h = -241.826 kJ/mol
% https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1#Thermo-Gas
T = Tu;
if (T < 1000)
    ic = 8;
else
    ic = 1;
end
% Oxygen molar values
cpmu(i,1) = Runiv*(O2c(ic) + O2c(ic+1)*T + O2c(ic+2)*T^2 + O2c(ic+3)*T^3 + O2c(ic+4)*T^4);                            % [J/molK]
cvmu(i,1) = cpmu(i,1) - Runiv;
hmu(i,1) = Runiv*T*(O2c(ic) + O2c(ic+1)/2*T + O2c(ic+2)/3*T^2 + O2c(ic+3)/4*T^3 + O2c(ic+4)/5*T^4 + O2c(ic+5)/T);     % [J/mol]
umu(i,1) = hmu(i,1) - Runiv*T;  % [J/mol]
somu(i,1) = Runiv*(O2c(ic)*log(T) + O2c(ic+1)*T + O2c(ic+2)/2*T^2 + O2c(ic+3)/3*T^3 + O2c(ic+4)/4*T^4 + O2c(ic+6));   % [J/molK]
% Oxygen mass values
cpu(i,1) = cpmu(i,1)/MO2;
cvu(i,1) = cvmu(i,1)/MO2;
hu(i,1) = hmu(i,1)/MO2;
uu(i,1) = umu(i,1)/MO2;
sou(i,1) = somu(i,1)/MO2;
% Nitrogen molar values
cpmu(i,2) = Runiv*(N2c(ic) + N2c(ic+1)*T + N2c(ic+2)*T^2 + N2c(ic+3)*T^3 + N2c(ic+4)*T^4);                          % [J/molK]
cvmu(i,2) = cpmu(i,2) - Runiv;
hmu(i,2) = Runiv*T*(N2c(ic) + N2c(ic+1)/2*T + N2c(ic+2)/3*T^2 + N2c(ic+3)/4*T^3 + N2c(ic+4)/5*T^4 + N2c(ic+5)/T);     % [J/mol]
umu(i,2) = hmu(i,2) - Runiv*T;  % [J/mol]
somu(i,2) = Runiv*(N2c(ic)*log(T) + N2c(ic+1)*T + N2c(ic+2)/2*T^2 + N2c(ic+3)/3*T^3 + N2c(ic+4)/4*T^4 + N2c(ic+6));   % [J/molK]
% Nitrogen mass values
cpu(i,2) = cpmu(i,2)/MN2;
cvu(i,2) = cvmu(i,2)/MN2;
hu(i,2) = hmu(i,2)/MN2;
uu(i,2) = umu(i,2)/MN2;
sou(i,2) = somu(i,2)/MN2;
% Argon molar values
cpmu(i,3) = Runiv*(Arc(ic) + Arc(ic+1)*T + Arc(ic+2)*T^2 + Arc(ic+3)*T^3 + Arc(ic+4)*T^4);                          % [J/molK]
cvmu(i,3) = cpmu(i,3) - Runiv;
hmu(i,3) = Runiv*T*(Arc(ic) + Arc(ic+1)/2*T + Arc(ic+2)/3*T^2 + Arc(ic+3)/4*T^3 + Arc(ic+4)/5*T^4 + Arc(ic+5)/T);     % [J/mol]
umu(i,3) = hmu(i,3) - Runiv*T;  % [J/mol]
somu(i,3) = Runiv*(Arc(ic)*log(T) + Arc(ic+1)*T + Arc(ic+2)/2*T^2 + Arc(ic+3)/3*T^3 + Arc(ic+4)/4*T^4 + Arc(ic+6));   % [J/molK]
% Argon mass values
cpu(i,3) = cpmu(i,3)/MAr;
cvu(i,3) = cvmu(i,3)/MAr;
hu(i,3) = hmu(i,3)/MAr;
uu(i,3) = umu(i,3)/MAr;
sou(i,3) = somu(i,3)/MAr;
% Carbon dioxide molar values
cpmu(i,4) = Runiv*(CO2c(ic) + CO2c(ic+1)*T + CO2c(ic+2)*T^2 + CO2c(ic+3)*T^3 + CO2c(ic+4)*T^4);                          % [J/molK]
cvmu(i,4) = cpmu(i,4) - Runiv;
hmu(i,4) = Runiv*T*(CO2c(ic) + CO2c(ic+1)/2*T + CO2c(ic+2)/3*T^2 + CO2c(ic+3)/4*T^3 + CO2c(ic+4)/5*T^4 + CO2c(ic+5)/T);     % [J/mol]
umu(i,4) = hmu(i,4) - Runiv*T;  % [J/mol]
somu(i,4) = Runiv*(CO2c(ic)*log(T) + CO2c(ic+1)*T + CO2c(ic+2)/2*T^2 + CO2c(ic+3)/3*T^3 + CO2c(ic+4)/4*T^4 + CO2c(ic+6));   % [J/molK]
% Carbon dioxide mass values
cpu(i,4) = cpmu(i,4)/MCO2;
cvu(i,4) = cvmu(i,4)/MCO2;
hu(i,4) = hmu(i,4)/MCO2;
uu(i,4) = umu(i,4)/MCO2;
sou(i,4) = somu(i,4)/MCO2;
% Water molar values
cpmu(i,5) = Runiv*(H2Oc(ic) + H2Oc(ic+1)*T + H2Oc(ic+2)*T^2 + H2Oc(ic+3)*T^3 + H2Oc(ic+4)*T^4);                          % [J/molK]
cvmu(i,5) = cpmu(i,5) - Runiv;
hmu(i,5) = Runiv*T*(H2Oc(ic) + H2Oc(ic+1)/2*T + H2Oc(ic+2)/3*T^2 + H2Oc(ic+3)/4*T^3 + H2Oc(ic+4)/5*T^4 + H2Oc(ic+5)/T);     % [J/mol]
umu(i,5) = hmu(i,5) - Runiv*T;  % [J/mol]
somu(i,5) = Runiv*(H2Oc(ic)*log(T) + H2Oc(ic+1)*T + H2Oc(ic+2)/2*T^2 + H2Oc(ic+3)/3*T^3 + H2Oc(ic+4)/4*T^4 + H2Oc(ic+6));   % [J/molK]
% Water mass values
cpu(i,5) = cpmu(i,5)/MH2O;
cvu(i,5) = cvmu(i,5)/MH2O;
hu(i,5) = hmu(i,5)/MH2O;
uu(i,5) = umu(i,5)/MH2O;
sou(i,5) = somu(i,5)/MH2O;
% Carbon monoxide molar values
cpmu(i,6) = Runiv*(COc(ic) + COc(ic+1)*T + COc(ic+2)*T^2 + COc(ic+3)*T^3 + COc(ic+4)*T^4);                          % [J/molK]
cvmu(i,6) = cpmu(i,6) - Runiv;
hmu(i,6) = Runiv*T*(COc(ic) + COc(ic+1)/2*T + COc(ic+2)/3*T^2 + COc(ic+3)/4*T^3 + COc(ic+4)/5*T^4 + COc(ic+5)/T);     % [J/mol]
umu(i,6) = hmu(i,6) - Runiv*T;  % [J/mol]
somu(i,6) = Runiv*(COc(ic)*log(T) + COc(ic+1)*T + COc(ic+2)/2*T^2 + COc(ic+3)/3*T^3 + COc(ic+4)/4*T^4 + COc(ic+6));   % [J/molK]
% Carbon monoxide mass values
cpu(i,6) = cpmu(i,6)/MCO;
cvu(i,6) = cvmu(i,6)/MCO;
hu(i,6) = hmu(i,6)/MCO;
uu(i,6) = umu(i,6)/MCO;
sou(i,6) = somu(i,6)/MCO;
% Hydrogen molar values
cpmu(i,7) = Runiv*(H2c(ic) + H2c(ic+1)*T + H2c(ic+2)*T^2 + H2c(ic+3)*T^3 + H2c(ic+4)*T^4);                          % [J/molK]
cvmu(i,7) = cpmu(i,7) - Runiv;
hmu(i,7) = Runiv*T*(H2c(ic) + H2c(ic+1)/2*T + H2c(ic+2)/3*T^2 + H2c(ic+3)/4*T^3 + H2c(ic+4)/5*T^4 + H2c(ic+5)/T);     % [J/mol]
umu(i,7) = hmu(i,7) - Runiv*T;  % [J/mol]
somu(i,7) = Runiv*(H2c(ic)*log(T) + H2c(ic+1)*T + H2c(ic+2)/2*T^2 + H2c(ic+3)/3*T^3 + H2c(ic+4)/4*T^4 + H2c(ic+6));   % [J/molK]
% Hydrogen mass values
cpu(i,7) = cpmu(i,7)/MH2;
cvu(i,7) = cvmu(i,7)/MH2;
hu(i,7) = hmu(i,7)/MH2;
uu(i,7) = umu(i,7)/MH2;
sou(i,7) = somu(i,7)/MH2;
% OH molar values
cpmu(i,8) = Runiv*(OHc(ic) + OHc(ic+1)*T + OHc(ic+2)*T^2 + OHc(ic+3)*T^3 + OHc(ic+4)*T^4);                          % [J/molK]
cvmu(i,8) = cpmu(i,8) - Runiv;
hmu(i,8) = Runiv*T*(OHc(ic) + OHc(ic+1)/2*T + OHc(ic+2)/3*T^2 + OHc(ic+3)/4*T^3 + OHc(ic+4)/5*T^4 + OHc(ic+5)/T);     % [J/mol]
umu(i,8) = hmu(i,8) - Runiv*T;  % [J/mol]
somu(i,8) = Runiv*(OHc(ic)*log(T) + OHc(ic+1)*T + OHc(ic+2)/2*T^2 + OHc(ic+3)/3*T^3 + OHc(ic+4)/4*T^4 + OHc(ic+6));   % [J/molK]
% OH mass values
cpu(i,8) = cpmu(i,8)/MOH;
cvu(i,8) = cvmu(i,8)/MOH;
hu(i,8) = hmu(i,8)/MOH;
uu(i,8) = umu(i,8)/MOH;
sou(i,8) = somu(i,8)/MOH;
% H molar values
cpmu(i,9) = Runiv*(Hc(ic) + Hc(ic+1)*T + Hc(ic+2)*T^2 + Hc(ic+3)*T^3 + Hc(ic+4)*T^4);                          % [J/molK]
cvmu(i,9) = cpmu(i,9) - Runiv;
hmu(i,9) = Runiv*T*(Hc(ic) + Hc(ic+1)/2*T + Hc(ic+2)/3*T^2 + Hc(ic+3)/4*T^3 + Hc(ic+4)/5*T^4 + Hc(ic+5)/T);     % [J/mol]
umu(i,9) = hmu(i,9) - Runiv*T;  % [J/mol]
somu(i,9) = Runiv*(Hc(ic)*log(T) + Hc(ic+1)*T + Hc(ic+2)/2*T^2 + Hc(ic+3)/3*T^3 + Hc(ic+4)/4*T^4 + Hc(ic+6));   % [J/molK]
% H mass values
cpu(i,9) = cpmu(i,9)/MH;
cvu(i,9) = cvmu(i,9)/MH;
hu(i,9) = hmu(i,9)/MH;
uu(i,9) = umu(i,9)/MH;
sou(i,9) = somu(i,9)/MH;
% O molar values
cpmu(i,10) = Runiv*(Oc(ic) + Oc(ic+1)*T + Oc(ic+2)*T^2 + Oc(ic+3)*T^3 + Oc(ic+4)*T^4);                          % [J/molK]
cvmu(i,10) = cpmu(i,10) - Runiv;
hmu(i,10) = Runiv*T*(Oc(ic) + Oc(ic+1)/2*T + Oc(ic+2)/3*T^2 + Oc(ic+3)/4*T^3 + Oc(ic+4)/5*T^4 + Oc(ic+5)/T);     % [J/mol]
umu(i,10) = hmu(i,10) - Runiv*T;  % [J/mol]
somu(i,10) = Runiv*(Oc(ic)*log(T) + Oc(ic+1)*T + Oc(ic+2)/2*T^2 + Oc(ic+3)/3*T^3 + Oc(ic+4)/4*T^4 + Oc(ic+6));   % [J/molK]
% O mass values
cpu(i,10) = cpmu(i,10)/MO;
cvu(i,10) = cvmu(i,10)/MO;
hu(i,10) = hmu(i,10)/MO;
uu(i,10) = umu(i,10)/MO;
sou(i,10) = somu(i,10)/MO;
%% Mixture properties
% Calculate the mixture molecular mass, total mass, gas constant and volume
% of the unburned gas
mu(i,NM+1) = 0;
cpmu(i,NM+1) = 0;
cvmu(i,NM+1) = 0;
hmu(i,NM+1) = 0;
umu(i,NM+1) = 0;
somu(i,NM+1) = 0;
for j=1:NM
    mu(i,NM+1) = mu(i,NM+1) + mu(i,j);
    cpmu(i,NM+1) = cpmu(i,NM+1) + Xu(i,j)*cpmu(i,j);        % [J/molK]
    cvmu(i,NM+1) = cvmu(i,NM+1) + Xu(i,j)*cvmu(i,j);        % [J/molK]
    hmu(i,NM+1) = hmu(i,NM+1) + Xu(i,j)*hmu(i,j);           % [J/mol]
    umu(i,NM+1) = umu(i,NM+1) + Xu(i,j)*umu(i,j);           % [J/mol]
    somu(i,NM+1) = somu(i,NM+1) + Xu(i,j)*somu(i,j);        % [J/molK]
end
Ru = Runiv/Mu;                      % Gas constant of unburned zone [J/kgK]
Vu = mu(i,NM+1)*Ru*Tu/pu;           % Volume[m]
cpu(i,NM+1) = cpmu(i,NM+1)/Mu;      % [J/kgK]
cvu(i,NM+1) = cvmu(i,NM+1)/Mu;      % [J/kgK]
hu(i,NM+1) = hmu(i,NM+1)/Mu;        % [J/kg]
uu(i,NM+1) = umu(i,NM+1)/Mu;        % [J/kg]
sou(i,NM+1) = somu(i,NM+1)/Mu;      % [J/kgK]




