function [mfa, Rfa, Vfa, cpmfa, cvmfa, hmfa, umfa, somfa, cpfa, cvfa, hfa, ufa, sofa] = ...
    faprop(Runiv, pref, pfa, Tfa, Mfa, mofa, fuelchem)
% Return only the so values of entropy, we will handle partial pressures
% later
% Determine the properties
mfa = mofa*Mfa;               % mass of added gas at IVC - includes residual & EGR
Rfa = Runiv/Mfa;                        % Gas constant of added fuel [J/kgK]
Vfa = mfa*Rfa*Tfa/pfa; % Volume of added gas [m3]
% CHEMKIN fits - starts at fuelchem(20) with the higher temperature
% values
T = Tfa;
if (T < 1000)
    ic = 27;
else
    ic = 20;
end
% DME @ 298.15K: h=-184.1 kJ/mol; cp=65.57 J/molK; 266.8 J/molK - slightly different
% than the CHEMKIN curve-fits
cpmfa = Runiv*(fuelchem(ic) + fuelchem(ic+1)*T + fuelchem(ic+2)*T^2 + fuelchem(ic+3)*T^3 + fuelchem(ic+4)*T^4);
cvmfa = cpmfa - Runiv;
hmfa = Runiv*T*(fuelchem(ic) + fuelchem(ic+1)/2*T + fuelchem(ic+2)/3*T^2 + fuelchem(ic+3)/4*T^3 + fuelchem(ic+4)/5*T^4 + fuelchem(ic+5)/T);     % [J/mol]
umfa = hmfa - Runiv*T;  % [J/mol]
somfa = Runiv*(fuelchem(ic)*log(T) + fuelchem(ic+1)*T + fuelchem(ic+2)/2*T^2 + fuelchem(ic+3)/3*T^3 + fuelchem(ic+4)/4*T^4 + fuelchem(ic+6));   % [J/molK]
% Here, we would use Dalton's law. But, for the assisted gas we are
% assuming there is only the blended fuel filling the particular
% volume. Thus, the mole fraction would be one and it is at the
% pressure of the chamber
%smfa = somfa - Runiv*log(pfa/pref); % [J/molK]
cpfa = cpmfa/Mfa;    % J/kgK
cvfa = cvmfa/Mfa;   % J/kgK
hfa = hmfa/Mfa;      % J/kg
ufa = umfa/Mfa;      % J/kg
sofa = somfa/Mfa;       % J/kgK
