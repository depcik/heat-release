function [mfl, Rfl, Vfl, cpmfl, cvmfl, hmfl, umfl, somfl, cpfl, cvfl, hfl, ufl, sofl] = ...
    flprop(Runiv, pref, pfl, Tfl, Mfl, mofl, fuelchem)
% Return only the so values of entropy; we will handle partial pressures
% separately in the different routines
% Determine the properties
mfl = mofl*Mfl;               % mass of added gas at IVC - includes residual & EGR
Rfl = Runiv/Mfl;                        % Gas constant of added fuel [J/kgK]
Vfl = mfl*Rfl*Tfl/pfl; % Volume of added gas [m3]
% CHEMKIN fits - starts at fuelchem(20) with the higher temperature
% values
T = Tfl;
if (T < 1000)
    ic = 8;
else
    ic = 1;
end
% Biodiesel heat of formation - https://www.sciencedirect.com/science/article/pii/S0009308409003521?via%3Dihub
cpmfl = Runiv*(fuelchem(ic) + fuelchem(ic+1)*T + fuelchem(ic+2)*T^2 + fuelchem(ic+3)*T^3 + fuelchem(ic+4)*T^4);
cvmfl = cpmfl - Runiv;
hmfl = Runiv*T*(fuelchem(ic) + fuelchem(ic+1)/2*T + fuelchem(ic+2)/3*T^2 + fuelchem(ic+3)/4*T^3 + fuelchem(ic+4)/5*T^4 + fuelchem(ic+5)/T);     % [J/mol]
umfl = hmfl - Runiv*T;  % [J/mol]
somfl = Runiv*(fuelchem(ic)*log(T) + fuelchem(ic+1)*T + fuelchem(ic+2)/2*T^2 + fuelchem(ic+3)/3*T^3 + fuelchem(ic+4)/4*T^4 + fuelchem(ic+6));   % [J/molK]
% Here, we would use Dalton's law. But, for the assisted gas we are
% assuming there is only the blended fuel filling the particular
% volume. Thus, the mole fraction would be one and it is at the
% pressure of the chamber
%smfl = somfl - Runiv*log(pfl/pref); % [J/molK]
cpfl = cpmfl/Mfl;    % J/kgK
cvfl = cvmfl/Mfl;   % J/kgK
hfl = hmfl/Mfl;      % J/kg
ufl = umfl/Mfl;      % J/kg
sofl = somfl/Mfl;       % J/kgK
