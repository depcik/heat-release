function HR1out = firstlawHR(iS, pdata, TCV, cpCV, RCV, XCV, YCV, UCV, dQHTfldt, dHindt, ...
    mu, mfa, mfl, mCV, inp, calc, cr, fuel, fuelchem)
% --- Routine specific variables
NM = inp(13);                               % Total number of chemical species
nps = calc(2);                              % Total number of data points
epseng = inp(78);                           % Conversion criteria on energies [J]
IVC = calc(3);                             % Index on IVC [-]
EVO = calc(4);                              % Index on EVO [-]
rps = calc(24);                             % Radians per second [rad/s] = [deg/s] * 3.14156 rad/180 deg
dps = calc(11);                             % Degrees per second [deg/s] = [rev/min]*[360 deg/rev]*[min/60 s]
bore = inp(1);                              % Bore [m]
Twall = inp(82);                            % Wall temperature [K]
Spm = calc(23);                             % Mean piston speed [m/s]
ahtc = inp(79);                             % Coefficient in front of heat transfer expression [-]
bhtc = inp(80);                             % Coefficient on Reynolds number for heat transfer [-]
chtc = inp(81);                             % Coefficient on Prandtl number of heat transfer [-]
sigma = inp(114);                           % Stefan-Boltzman constant [W/(m2*K^4)]
aw = inp(71);                               % Wall absorptivity [-]
dt = calc(12);                              % Time-step of data [s]
mfaadd = cr(5);                             % Added gas fuel mass per cycle [kg/cycle]
mfladd = cr(3);                             % Total liquid fuel mass per cycle [kg/cycle]
Qlhvfl = fuel(7);                           % Lower Heating Value [J/kg] of liquid fuel
Qlhvfa = fuel(26);                          % Lower Heating Value [J/kg] of added gas fuel
ncfl = inp(65);                             % Combustion efficiency [-] of liquid fuel
ncfa = inp(64);                             % Combustion efficiency [-] of added gas fuel
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

% --- Emissivity coefficients
A(1)=1.13263*10^-1;                         %supplied coefficient
A(2)=1.28673;                               %supplied coefficient
A(3)=-1.492;                                %supplied coefficient
A(4)=6.0441*10^-1;                          %supplied coefficient
B(1)=1.3018*10^-2;                          %supplied coefficient
B(2)=-6.813*10^-2;                          %supplied coefficient
B(3)=7.9223*10^-1;                          %supplied coefficient
B(4)=-5.40546*10^-1;                        %supplied coefficient
C(1)=-1.0781*10^-2;                         %supplied coefficient
C(2)=-5.821*10^-1;                          %supplied coefficient
C(3)=2.4541*10^-1;                          %supplied coefficient
C(4)=1.01909*10^-1;                         %supplied coefficient
D(1)=-8.8921*10^-3;                         %supplied coefficient
D(2)=2.05064*10^-1;                         %supplied coefficient
D(3)=-1.6229*10^-1;                         %supplied coefficient
D(4)=1.47167*10^-2;                         %supplied coefficient

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
HR1out = zeros(nps,30);                        % Output of heat release
% 1 - dWCV/dt [W]
% 2 - dWCV/deg [J/deg]
% 3 - hc [W/m2K]
% 4 - dQHTc/dt [W]
% 5 - dQHTc/deg [J/deg]
% 6 - emissivity [-]
% 7 - dQHTr/dt [W]
% 8 - dQHTr/deg [J/deg]
% 9 - dQHTfl/dt [W]
% 10 - dQHTfl/deg [J/deg]
% 11 - mdot*hin [W]
% 12 - mdot*hin/deg [J/deg]
% 13 - dUCV/dt [W]
% 14 - dUCV/deg [J/deg]
% 15 - dQHR/dt [W]
% 16 - dQHR/deg [J/deg]
% 17 - QHR [J]
% 18 - QHTc [J]
% 19 - QHTr [J]
% 20 - WCV [J]
% 21 - UCV [J]
% 22 - QHTfl [J]
% 23 - HIN [J]
% 24 - QHRsp [J]
% 25 - dQHRsp/dt [W]
% 26 - dQHRsp/deg [J/deg]
% 27 - dQHTcorrected/dt [W] - this is the heat transfer needed to balance
% the Conservation of Energy equation at each step
% 28 - hcorrected [W/m2K] - this would be the resulting heat transfer
% coefficient if it was all convection
% 29 - dQHTcorrected/deg [J/deg]
% 30 - QHT,corrected [J]
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

% --- Calculated version of heat transfer
iter = 1;
epsht = 1000;
while (epsht > epseng)
    for j=IVC+1:EVO     % IVC to EVO
        % --- Calculate the work using the trapezoidal rule
        % Have the derivative of volume with respect to crank angle in
        % another routine
        % dWCV/dt & dWCV/deg
        HR1out(j,1) = 0.5*(pdata(j,13)+pdata(j-1,13))*pdata(j,7)*rps;    % [Pa]*[m3/rad]*[rad/s] = [W]
        HR1out(j,2) = HR1out(j,1)/dps;                                % [J/deg]
        % --- Calculate the heat transfer
        % For simplicity, I used Sutherland's equations for viscosity and
        % thermal conductivity
        % Instantaneous cylinder height is the characteristic length unless
        % it is greater than the bore/2
        % Mean piston speed is the characteristic velocity - for now
        mug = (1.4592*(TCV(j)^1.5))/(109.10 + TCV(j))*(10^-6);   % [kg/(m*s)]
        kg = (2.3340*(TCV(j)^1.5))/(164.54 + TCV(j))*(10^-3);    % [W/(m*K)]
        Pr = cpCV(j)*mug/kg;                     % Prandtl number [-]
        rhoCV = pdata(j,13)/(RCV(j)*TCV(j));    % Density [kg/m3]
        % Instantaneous cylinder height
        % Cylinder volume / piston surface area
        hl = pdata(j,5)/(pi()*bore^2/4);
        if (hl < (bore/2))
            Lc = hl;
        else
            Lc = bore/2;
        end
        Re = rhoCV*Spm*Lc/mug;          % Reynolds number [-]
        Nu = ahtc*(Re^bhtc)*(Pr^chtc); % Nusselt number [-]
        %hc(j) = kg*Nu/Lc;               % Convective heat transfer coefficient [W/m2K]
        HR1out(j,3) = kg*Nu/Lc;               % Convective heat transfer coefficient [W/m2K]
        % --- Heat transfer due to convection
        % dQHTc/dt
        HR1out(j,4) = HR1out(j,3)*pdata(j,8)*(Twall-TCV(j));   % [W]
        HR1out(j,5) = HR1out(j,4)/dps;                  % [J/deg]
        % --- Heat transfer due to radiation
        % dQHTr/dt
        Lb = pdata(j,5)/pdata(j,8);     % Beam length [m] - cylinder volume to surface area
        Te = TCV(j)/1000;               % Weighted temperature for coefficients
        for jj=1:4
            %calculates coefficients
            a(jj)=A(jj)+B(jj)*Te+C(jj)*Te^2+D(jj)*Te^3;
        end
        % Find the partial pressure of carbon dioxide and water vapor
        % 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
        % 11-added gas, 12-liquid fuel
        % Note: Dr. Mattson's code also has an adiabatic flame temperature
        % calculation embedded and that can be used to determine the
        % radiation based on a flme
        pgCO2H2O = (XCV(j,4) + XCV(j,5))*pdata(j,13)*10^-5;     % Partial pressure is in [bar] for the correlation
        HR1out(j,6) = a(1)+a(2)*(pgCO2H2O*Lb)+a(3)*(pgCO2H2O*Lb)^2+a(4)*(pgCO2H2O*Lb)^3;
        % Remember, heat transfer in is (+)
        % If the cylinder temperature is greater than the wall temperature,
        % heat is leaving the gas
        HR1out(j,7) = sigma*pdata(j,8)*(aw*Twall^4 - HR1out(j,6)*TCV(j)^4);        % [W/(m2*K4)]*[m2]*[K4] = [W]
        % Including adiabatic flame temperature - not included in model
        % dQHTrdt(j) = inp(114)*pdata(j,8)*(inp(71)*inp(82)^4 - 0.6*TAFT^4);
        HR1out(j,8) = HR1out(j,7)/dps;                   % [J/deg]
        % --- Heat transfer due to fuel injection
        % dQHTfl/dt
        if (iS(7) == 0)
            HR1out(j,9) = dQHTfldt(j);
            HR1out(j,10) = HR1out(j,9)/dps;                   % [J/deg]
        else
            HR1out(j,9) = 0;
            HR1out(j,10) = 0;
        end
        % --- Enthalpy addition
        % dHin/dt
        if (iS(7) == 0)
            HR1out(j,11) = dHindt(j);
            HR1out(j,12) = HR1out(j,11)/dps;                       % [J/deg]
        else
            HR1out(j,11) = 0;
            HR1out(j,12) = 0;
        end
        % --- Internal energy
        % dUCV/dt
        HR1out(j,13) = (UCV(j)-UCV(j-1))/dt;               % [J]/[s] = [W]
        HR1out(j,14) = HR1out(j,13)/dps;                     % [J/deg]
        % --- Heat Release
        HR1out(j,15) = HR1out(j,13) + HR1out(j,1) - HR1out(j,4) - HR1out(j,7) - HR1out(j,9) - HR1out(j,11);    % [W]
        HR1out(j,16) = HR1out(j,15)/dps;                % [J/deg]
        % --- Sum of the different components
        HR1out(j,17) = HR1out(j-1,17) + HR1out(j,15)*dt;                  % [J]
        HR1out(j,18) = HR1out(j-1,18) + HR1out(j,4)*dt;                 % [J]
        HR1out(j,19) = HR1out(j-1,19) + HR1out(j,7)*dt;              % [J]
        HR1out(j,20) = HR1out(j-1,20) + HR1out(j,1)*dt;                  % [J]
        HR1out(j,21) = HR1out(j-1,21) + HR1out(j,13)*dt;                 % [J]
        HR1out(j,22) = HR1out(j-1,22) + HR1out(j,9)*dt;             % [J]
        HR1out(j,23) = HR1out(j-1,23) + HR1out(j,11)*dt;                   % [J]
    end

    % --- Theoretical version of heat transfer
    % Heat release based on change in chemical species
    HR1out(:,24) = zeros;
    HR1out(:,25) = zeros;
    HR1out(:,26) = zeros;
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
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, uu, ~, ~, ~] = ...
            unburnprop(j, Runiv, pref, pdata(j,13), TCV(j), MM, ...
            mou_O2, mou_N2, mou_Ar, mou_CO2, mou_H2O, mou_CO, mou_H2, mou_OH, mou_H, mou_O, ...
            O2c, N2c, Arc, CO2c, H2Oc, COc, H2c, OHc, Hc, Oc, ...
            mu, Xu, Yu, cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou);
        % Added gas
        mofa = mfa(j)/Mfa;
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ufa, ~] = faprop(Runiv, pref, pdata(j,13), TCV(j), Mfa, mofa, fuelchem);
        % Liquid fuel
        mofl = mfl(j)/Mfl;
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ufl, ~] = flprop(Runiv, pref, pdata(j,13), TCV(j), Mfl, mofl, fuelchem);
        % Let's compute the change in mass fractions
        % 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
        % 11-added gas, 12-liquid fuel
        for k=1:NM+2
            dYdt = (YCV(j,k)-YCV(j-1,k))/dt;
            if (k < NM+1)
                HR1out(j,25) = HR1out(j,25) - mCV(j)*uu(j,k)*dYdt;
            elseif (k == NM+1)
                HR1out(j,25) = HR1out(j,25) - mCV(j)*ufa*dYdt;
            else
                HR1out(j,25) = HR1out(j,25) - mCV(j)*ufl*dYdt;
            end
        end
        HR1out(j,26) = HR1out(j,25)/dps;
        HR1out(j,24) = HR1out(j-1,24) + HR1out(j,25)*dt;
    end
    % The total from the two fuels
    % A simple check based on combustion efficiencies
    % The values should be close to this value
    QHRs = mfaadd*ncfa*Qlhvfa + mfladd*ncfl*Qlhvfl;
    % The total heat release between the three options should line up
    QHRc = HR1out(EVO,17); % calculated heat release
    QHRt = HR1out(EVO,24); % theoretical heat release from the species changing
    % --- NR iteration
    if (iter == 1)
        QHRold = QHRc;
        ahtcold = ahtc;
        ahtc = 1.1*ahtcold;
        iter = iter+1;
    else
        ahtcnew = ahtc - (QHRc-QHRt)/((QHRc-QHRold)/(ahtc-ahtcold));
        %ahtcnew = ahtc - (QHRc-QHRs)/((QHRc-QHRold)/(ahtc-ahtcold));
        QHRold = QHRc;
        ahtcold = ahtc;
        ahtc = ahtcnew;
    end
    epsht = abs(QHRt-QHRc);
    %epsht = abs(QHRs-QHRc);
end
[QHRs QHRc QHRt]
% --- Ideally, the plot of the heat release from theory and that from the
% calculated in-cylinder trace should line up. They do not. This has to do
% with inaccuracies in pressure, temperature, heat transfer, etc. Heat
% transfer is the largest one. What would heat transfer have to be to get
% them to align assuming that all the other values are correct
for j=IVC+1:EVO
    % Heat release at this point from theory - HR1out(:,25)
    % Heat release calculated - HR1out(:,15)
    % dQHT/dt = dUCV/dt - (dQHR/dt)theory - dQHT,fl/dt + dWCV/dt - sum(mdotin*hin)
    HR1out(j,27) = HR1out(j,13) - HR1out(j,25) - HR1out(j,9) + HR1out(j,1) - HR1out(j,11);  % [W]
    % If we assume all this heat transfer is from convection and that to
    % the liquid fuel is exact
    % dQHTc/dt = hc*As*(Twall-TCV)
    HR1out(j,28) = HR1out(j,27)/(Twall-TCV(j))/pdata(j,8);
    % Corrected heat transfer [J/deg]
    HR1out(j,29) = HR1out(j,27)/dps;
    % Total corrected heat transfer [J]
    HR1out(j,30) = HR1out(j-1,30) + HR1out(j,27)*dt;        % [J]
end
% The value of convective heat transfer here allows for the conservation of
% energy equation to balance at each step. This will also result in the
% same heat release between theory and the data.