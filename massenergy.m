function [inp, mu, mb, mfa, mfl, mCV, YCV, XCV, RCV, TCV, cpCV, cvCV, UCV, HCV, SCV, ECV, ...
    Tfl, cpfl, dQHTfldt, dHindt, dSHTfldt, dSindt, dEHTfldt, dEindt] = massenergy(iS, inp, calc, fuel, fuelchem, fuelinj, cr, pdata)
% After going through this multiple times, the only way the ideal gas law
% and internal energy balance is when all the zones are at the same
% temperature. This is because the internal energy and temperature of the
% control volume are based on these species.
global varode  % Why use a global variable here? GNU Octave does not have the same capabilities as the Matlab ODE call
% Thus, this is a fix to use both Matlab and GNU Octave
%% --- Zero arrays ---%
% Chemical species
% 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
NM = inp(13);   % total number of species
% Unburned species
% O2, N2, and Ar come from air (no CO2 & no H2O in the air)
% CO2 and H2O come from residual & EGR
% CO, H2, OH, H, and O come from partial combustion
mu = zeros(calc(2),NM+1);
mou_O2 = zeros(calc(2),1);
mou_N2 = zeros(calc(2),1);
mou_Ar = zeros(calc(2),1);
mou_CO2 = zeros(calc(2),1);
mou_H2O = zeros(calc(2),1);
mou_CO = zeros(calc(2),1);
mou_H2 = zeros(calc(2),1);
mou_OH = zeros(calc(2),1);
mou_H = zeros(calc(2),1);
mou_O = zeros(calc(2),1);
Xu = zeros(calc(2),NM+1);
Yu = zeros(calc(2),NM+1);
cpmu = zeros(calc(2),NM+1);
cvmu = zeros(calc(2),NM+1);
hmu = zeros(calc(2),NM+1);
umu = zeros(calc(2),NM+1);
somu = zeros(calc(2),NM+1);
cpu = zeros(calc(2),NM+1);
cvu = zeros(calc(2),NM+1);
hu = zeros(calc(2),NM+1);
uu = zeros(calc(2),NM+1);
sou = zeros(calc(2),NM+1);
% Burned species
% There is no oxygen in the burned species
% N2 and Ar come from the air - they are with oxygen and transition from
% the unburned to burned mixture at the same time; i.e., move with it
% CO2 and H2O come from the combustion reactions
% Partial products of combustion (CO, H2, OH, H, and O) will not show up as
% part of the burned mixture. We might need the arrays, however, to use the
% same properties routines
mb = zeros(calc(2),NM+1);
mob_O2 = zeros(calc(2),1);
mob_N2 = zeros(calc(2),1);
mob_Ar = zeros(calc(2),1);
mob_CO2 = zeros(calc(2),1);
mob_H2O = zeros(calc(2),1);
mob_CO = zeros(calc(2),1);
mob_H2 = zeros(calc(2),1);
mob_OH = zeros(calc(2),1);
mob_H = zeros(calc(2),1);
mob_O = zeros(calc(2),1);
Xb = zeros(calc(2),NM+1);
Yb = zeros(calc(2),NM+1);
cpmb = zeros(calc(2),NM+1);
cvmb = zeros(calc(2),NM+1);
hmb = zeros(calc(2),NM+1);
umb = zeros(calc(2),NM+1);
somb = zeros(calc(2),NM+1);
cpb = zeros(calc(2),NM+1);
cvb = zeros(calc(2),NM+1);
hb = zeros(calc(2),NM+1);
ub = zeros(calc(2),NM+1);
sob = zeros(calc(2),NM+1);
% Added gaseous fuel
mofa = zeros(calc(2),1);
mfa = zeros(calc(2),1);
Vfa = zeros(calc(2),1);
cpmfa = zeros(calc(2),1);
cvmfa = zeros(calc(2),1);
hmfa = zeros(calc(2),1);
umfa = zeros(calc(2),1);
somfa = zeros(calc(2),1);
cpfa = zeros(calc(2),1);
cvfa = zeros(calc(2),1);
hfa = zeros(calc(2),1);
ufa = zeros(calc(2),1);
sofa = zeros(calc(2),1);
% Injected fuel
mfl = zeros(calc(2),1);
mofl = zeros(calc(2),1);
Vfl = zeros(calc(2),1);
cpmfl = zeros(calc(2),1);
cvmfl = zeros(calc(2),1);
hmfl = zeros(calc(2),1);
umfl = zeros(calc(2),1);
somfl = zeros(calc(2),1);
cvfl = zeros(calc(2),1);
hfl = zeros(calc(2),1);
ufl = zeros(calc(2),1);
sofl = zeros(calc(2),1);
Tfl = zeros(calc(2),1);
cpfl = zeros(calc(2),1);
% Total control volume
% 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
% 11-fa (added gas fuel), 12-fl (injected liquid fuel)
XCV = zeros(calc(2),NM+2);
YCV = zeros(calc(2),NM+2);
dQHTfldt = zeros(calc(2),1);       % Heat transfer due to vaporization and heating of droplet [W]
dHindt = zeros(calc(2),1);         % Enthalpy addition to the control volume [W]
dSindt = zeros(calc(2),1);         % Entropy addition to the control volume [W/K] = [kg/s]*[J/kg*K]
dSHTfldt = zeros(calc(2),1);       % Entropy change due to heat transfer to liquid fuel [W/K]
dEHTfldt = zeros(calc(2),1);       % Exergy change due to heat transfer to the liquid fuel [W]
dEindt = zeros(calc(2),1);         % Specific exergy flow added to the control volume [W]

%% --- Get the CHEMKIN data curve fits for the different constituents ---%
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

%% --- Routine Variables --- %
MO2 = 2*inp(102);               % Molecular mass of O2 [kg/mol]
MN2 = 2*inp(103);               % Molecular mass of N2 [kg/mol]
MAr = inp(105);                 % Molecular mass of Ar [kg/mol]
MCO2 = inp(101) + 2*inp(102);   % Molecular mass of CO2 [kg/mol]
MH2O = 2*inp(104) + inp(102);   % Molecular mass of H2O [kg/mol]
MCO = inp(101) + inp(102);      % Molecular mass of CO [kg/mol]
MOH = inp(102) + inp(104);      % Molecular mass of OH [kg/mol]
MH2 = 2*inp(104);               % Molecular mass of H2 [kg/mol]
MH = inp(104);                  % Molecular mass of H [kg/mol]
MO = inp(102);                  % Molecular mass of O [kg/mol]
% index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O, NM+1-Mixture
MM = [MO2 MN2 MAr MCO2 MH2O MCO MH2 MOH MH MO]; % Molecular mass array [kg/mol]
Mfa = fuel(24);                 % Molecular mass of added fuel [kg/mol]
Mfl = fuel(5);                  % Molecular mass of liquid fuel
Runiv = inp(100);               % Universal gas constant [J/molK]
pref = inp(83);                 % CHEMKIN Reference Pressure [Pa]
dt = calc(12);                  % Time-step of the data [s]
Kfa = inp(73);                  % Added gas fuel burn rate parameter [m3/(kg*s)]
Efa = inp(74);                  % Added gas fuel burn rate activation energy [J/molK]
Kfl = inp(75);                  % Liquid fuel burn rate parameter [m3/(kg*s)]
Efl = inp(76);                  % Liquid fuel burn rate activation energy [J/molK]
% n2 = cr(50);                    % Added gas - local reaction: O2 [mol]
% n3 = cr(51);                    % Added gas - local reaction: CO2 [mol]
% n4 = cr(52);                    % Added gas - local reaction: H2O [mol]
% n5 = cr(53);                    % Added gas - local reaction: N2 [mol]
% l2 = cr(40);                    % Liquid fuel - local reaction: O2 [mol]
% l3 = cr(41);                    % Liquid fuel - local reaction: CO2 [mol]
% l4 = cr(42);                    % Liquid fuel - local reaction: H2O [mol]
% l5 = cr(43);                    % Liquid fuel - local reaction: N2 [mol]
delta = cr(27);                 % Air unburned - global reaction [mol]
g2 = cr(28);                    % Air unburned - global reaction: O2 [-]
g3 = cr(29);                    % Air unburned - global reaction: N2 [-]
g4 = cr(30);                    % Air unburned - global reaction: Ar [-]
g5 = cr(31);                    % Res+EGR - global reaction: O2 [-]
g6 = cr(32);                    % Res+EGR - global reaction: CO2 [-]
g7 = cr(33);                    % Res+EGR - global reaction: H2O [-]
g8 = cr(34);                    % Res+EGR - global reaction: N2 [-]
g9 = cr(35);                    % Res+EGR - global reaction: Ar [-]
epsilon = cr(26);               % Epsilon component for Res+EGR
Xig = cr(21);                   % Moles of gaseous fuel [-]
Zetag = cr(23);                 % Moles of gaseous fuel left after combustion [-]
etac_fat = inp(64);             % Target combustion efficiency of added gas [-]
etac_flt = inp(65);             % Target combustion efficiency of liquid fuel [-]
To = inp(84);                   % Exergy reference temperature [K]
po = inp(85);                   % Exergy reference pressure [Pa]

% We have four zones: unburned, burned, added gas, liquid fuel (as a gas)

%% --- Set up indexes
% Calculations begin at IVC and end at EVO
IVC = calc(3);      % index of IVC
EVO = calc(4);      % index of EVO
IGNfa = calc(19);   % index of ignition point of gaseous fuel
IGNfl = calc(21);   % index of ignition point of liquid fuel
INJ1 = calc(5);     % index of first injection
%INJ2 = calc(6);     % index of second injection
%INJ3 = calc(7);     % index of third injection
%INJ4 = calc(8);     % index of fourth injection
%INJ5 = calc(9);     % index of fifth injection

%% --- Determine the unburned massses of O2, N2, Ar, CO2, and H2O at IVC
% At IVC, we can have added gas, air, residual, and EGR. The air, residual,
% and EGR are all part of the unburned mixture. The added gas is its own
% entity (fa). Let's set up the masses of these components at IVC
% oxygen, nitrogen, argon, carbon dioxide, and water
mou_O2(IVC) = cr(27)*cr(28) + cr(26)*cr(31);     % moles of O2 (unburned) at IVC
mou_N2(IVC) = cr(27)*cr(29) + cr(26)*cr(34);     % moles of N2 (unburned at IVC)
mou_Ar(IVC) = cr(27)*cr(30) + cr(26)*cr(35);     % moles of Ar (unburned at IVC)
mou_CO2(IVC) = cr(26)*cr(32);                    % moles of CO2 (unburned at IVC)
mou_H2O(IVC) = cr(26)*cr(33);                    % moles of H2O (unburned at IVC)
mou_CO(IVC) = 0;            % moles of CO (unburned at IVC)
mou_OH(IVC) = 0;            % moles of OH (unburned at IVC)
mou_H2(IVC) = 0;            % moles of H2 (unburned at IVC)
mou_O(IVC) = 0;             % moles of O (unburned at IVC)
mou_H(IVC) = 0;             % moles of H (unburned at IVC)
Tu(IVC) = cr(10);   % Unburned gas temperature [K] - Assumed that of temperature at IVC
% Call the unburned properties routine
% index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O, NM+1-Mixture
[mu, Xu, Yu, Mu(IVC), cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou, Ru(IVC), Vu(IVC)] = ...
    unburnprop(IVC, Runiv, pref, pdata(IVC,13), Tu(IVC), MM, ...
    mou_O2(IVC), mou_N2(IVC), mou_Ar(IVC), mou_CO2(IVC), mou_H2O(IVC), mou_CO(IVC), mou_H2(IVC), mou_OH(IVC), mou_H(IVC), mou_O(IVC), ...
    O2c, N2c, Arc, CO2c, H2Oc, COc, H2c, OHc, Hc, Oc, ...
    mu, Xu, Yu, cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou);

%% --- Check for an added gas and specify its mass, temperature, gas constant, and volume
% Also get the constant pressure specific heat, internal energy, enthalpy, and entropy
if (iS(2) == 1)
    mofa(IVC) = cr(21) + cr(23)*cr(26);     % moles of added gas at IVC - includes residual and EGR
    Tfa(IVC) = cr(10);                      % Added gas temperature [K]
    % Call the added gas properties routine
    [mfa(IVC), Rfa, Vfa(IVC), cpmfa(IVC), cvmfa(IVC), hmfa(IVC), umfa(IVC), somfa(IVC), ...
        cpfa(IVC), cvfa(IVC), hfa(IVC), ufa(IVC), sofa(IVC)] = faprop(Runiv, pref, pdata(IVC,13), Tfa(IVC), Mfa, mofa(IVC), fuelchem);
else
    mofa(IVC) = 0;
    mfa(IVC) = 0;
    Tfa(IVC) = 0;
    Vfa(IVC) = 0;
    cpmfa(IVC) = 0;
    hmfa(IVC) = 0;
    umfa(IVC) = 0;
    somfa(IVC) = 0;
    somfa(IVC) = 0;
    cpfa(IVC) = 0;
    hfa(IVC) = 0;
    ufa(IVC) = 0;
    sofa(IVC) = 0;
end

%% --- Total Control Volume --- %
[mCV(IVC), YCV, XCV, RCV(IVC), TCV(IVC), cpCV(IVC), cvCV(IVC), UCV(IVC), HCV(IVC), SCV(IVC), ECV(IVC)] = ...
    cvprop(IVC, inp, Runiv, pref, pdata(IVC,13), pdata(IVC,5), YCV, XCV, MM, Mfa, Mfl, ...
    mu, mfa(IVC), mfl(IVC), mb, O2c, N2c, Arc, CO2c, H2Oc, COc, H2c, OHc, Hc, Oc, fuel, fuelchem);

%% --- Iteration on combustion efficiency --- %%
epsetac = 100;
itetac = 1;
while (epsetac > inp(77))
    %% --- Now we proceed based on what we have in the mixture at the time of IVC
    % Do we have an added gas?
    % We only need one routine and we just check when combustion happens
    if (iS(2) == 1) % Yes
        indx = IVC+1;
        while (indx <= EVO)
            itit = 0;           % Check for first time in the routine
            epsit = 1000;
            while ((epsit > 1e-6) && (itit < 10))
                if (itit == 0)      % If first time, use old values
                    rhoCV = pdata(indx-1,13)/(RCV(indx-1)*TCV(indx-1));
                    YfaCV = YCV(indx-1,NM+1);
                    YflCV = YCV(indx-1,NM+2);
                    YO2CV = YCV(indx-1,1);
                    VolCV = pdata(indx-1,5);
                    TCVi = TCV(indx-1);
                    itit = 1;
                else                % Otherwise, use the values to iterate upon
                    rhoCV = pdata(indx,13)/(RCV(indx)*TCV(indx));
                    YfaCV = YCV(indx,NM+1);
                    YflCV = YCV(indx,NM+2);
                    YO2CV = YCV(indx,1);
                    VolCV = pdata(indx,5);
                    TCVi = TCV(indx);
                    itit = itit + 1;
                end
                %--- Check for the ignition of the fuels
                if (indx < IGNfa) && (indx < IGNfl)     % neither fuel has started burning
                    iofa = 0;
                    iofl = 0;
                elseif (indx >= IGNfa) && (indx < IGNfl) % added gas fuel has started burning, but not the direct injected fuel
                    iofa = 1;
                    iofl = 0;
                elseif (indx < IGNfa) && (indx >= IGNfl) % added gas fuel is not burning, but the direct injected fuel is
                    iofa = 0;
                    iofl = 1;
                else  % both fuels are burning
                    iofa = 1;
                    iofl = 1;
                end
                %--- Run the ode solver
                % Set up the ODE solver options
                options = odeset('RelTol', 1e-5, ...    % A relative error tolerance that applies to all components
                                'AbsTol', 1e-10, ...    % The absolute error tolerance that applies to all components
                                'MaxStep', dt, ...      % An upper bound on the magnitude of the step size
                                'Refine', 1);           % Produces smoother output
                % Let us solve for all masses (unburned, burned, added
                % gas and liquid fuel) in the combustion routine over
                % the time-step
                % index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
                for kk=1:NM
                    Zo(kk) = mu(indx-1,kk);
                    Zo(kk+NM) = mb(indx-1,kk);
                end
                Zo(2*NM+1) = mfa(indx-1);
                Zo(2*NM+2) = mfl(indx-1);
                % Variables to call into ode
                varode = [iS(8), iofa, iofl, NM, fuel(20), fuel(21), fuel(22), fuel(23), ...
                    fuel(1), fuel(2), fuel(3), fuel(4), ...
                    Kfa, Efa, Mfa, Kfl, Efl, Mfl, rhoCV, VolCV, Runiv, TCVi, MM ...
                    delta, epsilon, g2, g3, g4, g5, g6, g7, g8, g9, fuelinj(indx), inp(88), inp(89)];
                %[t, Z] = ode15s(@combustfafl, [0 dt], Zo, options, varode);
                [t, Z] = ode15s(@combustfafl, [0 dt], Zo, options);
                % Calculate the size of the arrays of the output
                ln = size(Z);
                % Fill the output masses
                for kk=1:NM
                    mu(indx,kk) = Z(ln(1),kk);
                    mb(indx,kk) = Z(ln(1),kk+NM);
                end
                mfa(indx) = Z(ln(1),2*NM+1);
                mofa(indx) = mfa(indx)/Mfa;
                mfl(indx) = Z(ln(1),2*NM+2);
                mofl(indx) = mfl(indx)/Mfl;
                % Fill the output molar values
                % index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
                % Unburned side
                mou_O2(indx) = mu(indx,1)/MO2;      % [kg] / [kg/mol] = [mol]
                mou_N2(indx) = mu(indx,2)/MN2;
                mou_Ar(indx) = mu(indx,3)/MAr;
                mou_CO2(indx) = mu(indx,4)/MCO2;
                mou_H2O(indx) = mu(indx,5)/MH2O;
                mou_CO(indx) = mu(indx,6)/MCO;
                mou_H2(indx) = mu(indx,7)/MH2;
                mou_OH(indx) = mu(indx,8)/MOH;
                mou_H(indx) = mu(indx,9)/MH;
                mou_O(indx) = mu(indx,10)/MO;
                mu(indx,NM+1) = mu(indx,1)+mu(indx,2)+mu(indx,3)+mu(indx,4)+mu(indx,5)+...
                    mu(indx,6)+mu(indx,7)+mu(indx,8)+mu(indx,9)+mu(indx,10);
                % Burned side
                mob_O2(indx) = mb(indx,1)/MO2;      % [kg] / [kg/mol] = [mol]
                mob_N2(indx) = mb(indx,2)/MN2;
                mob_Ar(indx) = mb(indx,3)/MAr;
                mob_CO2(indx) = mb(indx,4)/MCO2;
                mob_H2O(indx) = mb(indx,5)/MH2O;
                mob_CO(indx) = mb(indx,6)/MCO;
                mob_H2(indx) = mb(indx,7)/MH2;
                mob_OH(indx) = mb(indx,8)/MOH;
                mob_H(indx) = mb(indx,9)/MH;
                mob_O(indx) = mb(indx,10)/MO;
                mb(indx,NM+1) = mb(indx,1)+mb(indx,2)+mb(indx,3)+mb(indx,4)+mb(indx,5)+...
                    mb(indx,6)+mb(indx,7)+mb(indx,8)+mb(indx,9)+mb(indx,10);
                % Total control volume
                [mCV(indx), YCV, XCV, RCV(indx), TCV(indx), cpCV(indx), cvCV(indx), UCV(indx), HCV(indx), SCV(indx), ECV(indx)] = ...
                    cvprop(indx, inp, Runiv, pref, pdata(indx,13), pdata(indx,5), YCV, XCV, MM, Mfa, Mfl, ...
                    mu, mfa(indx), mfl(indx), mb, O2c, N2c, Arc, CO2c, H2Oc, COc, H2c, OHc, Hc, Oc, fuel, fuelchem);
                % To balance the ideal gas law and the internal energy,
                % the temperatures all need to be the same
                Tu(indx) = TCV(indx);
                Tfa(indx) = TCV(indx);
                Tfl(indx) = TCV(indx);
                Tb(indx) = TCV(indx);
                % Check internal energy and ideal gas law
                % Unburned mixture
                [mu, Xu, Yu, Mu(indx), cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou, Ru(indx), Vu(indx)] = ...
                    unburnprop(indx, Runiv, pref, pdata(indx,13), Tu(indx), MM, ...
                    mou_O2(indx), mou_N2(indx), mou_Ar(indx), mou_CO2(indx), mou_H2O(indx), mou_CO(indx), mou_H2(indx), mou_OH(indx), mou_H(indx), mou_O(indx), ...
                    O2c, N2c, Arc, CO2c, H2Oc, COc, H2c, OHc, Hc, Oc, ...
                    mu, Xu, Yu, cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou);
                Uu = mu(indx,NM+1)*uu(indx,NM+1);
                % Added gas
                [mfa(indx), Rfa, Vfa(indx), cpmfa(indx), cvmfa(indx), hmfa(indx), umfa(indx), somfa(indx), ...
                    cpfa(indx), cvfa(indx), hfa(indx), ufa(indx), sofa(indx)] = faprop(Runiv, pref, pdata(indx,13), Tfa(indx), Mfa, mofa(indx), fuelchem);
                Ufa = mfa(indx)*ufa(indx);
                % Liquid fuel
                [mfl(indx), Rfl, Vfl(indx), cpmfl(indx), cvmfl(indx), hmfl(indx), umfl(indx), somfl(indx), ...
                    cpfl(indx), cvfl(indx), hfl(indx), ufl(indx), sofl(indx)] = flprop(Runiv, pref, pdata(indx,13), Tfl(indx), Mfl, mofl(indx), fuelchem);
                Ufl = mfl(indx)*ufl(indx);
                % Find the internal energy of the burned zone
                % We can use the same routine as the unburned properties
                % Only need to calculate if the fuel(s) have started
                % burning
                if (mb(indx,NM+1) > 0)
                    [mb, Xb, Yb, Mb(indx), cpmb, cvmb, hmb, umb, somb, cpb, cvb, hb, ub, sob, Rb(indx), Vb(indx)] = ...
                        unburnprop(indx, Runiv, pref, pdata(indx,13), Tb(indx), MM, ...
                        mob_O2(indx), mob_N2(indx), mob_Ar(indx), mob_CO2(indx), mob_H2O(indx), mob_CO(indx), mob_H2(indx), mob_OH(indx), mob_H(indx), mob_O(indx), ...
                        O2c, N2c, Arc, CO2c, H2Oc, COc, H2c, OHc, Hc, Oc, ...
                        mb, Xb, Yb, cpmb, cvmb, hmb, umb, somb, cpb, cvb, hb, ub, sob);
                    Ub = mb(indx,NM+1)*ub(indx,NM+1);
                else
                    Ub = 0;
                    Rb(indx) = 0;
                end
                % Checks
                g1 = UCV(indx) - Uu - Ufa - Ufl - Ub;
                g2 = mCV(indx)*RCV(indx)*TCV(indx) - mu(indx,NM+1)*Ru(indx)*Tu(indx) - mfa(indx)*Rfa*Tfa(indx) - mfl(indx)*Rfl*Tfl(indx) - mb(indx,NM+1)*Rb(indx)*Tb(indx);
                % Heat transfer of liquid fuel
                dQHTfldt(indx) = fuelinj(indx)*(fuel(10)*(inp(25)-fuel(17)) - fuel(6) + cpfl(indx)*(fuel(17)-Tfl(indx)));
                % Entropy due to this heat transfer
                Tflavg = (0.5*(inp(25)+fuel(17))*fuel(10)*(inp(25)-fuel(17)) - fuel(17)*fuel(6) + ...
                    0.5*(fuel(17)+Tfl(indx))*(cpfl(indx)*(fuel(17)-Tfl(indx))))/ ...
                    (fuel(10)*(inp(25)-fuel(17)) - fuel(6) + cpfl(indx)*(fuel(17)-Tfl(indx)));
                dSHTfldt(indx) = dQHTfldt(indx)/Tflavg;     % [W/K]
                % Enthalpy addition into the control volume
                [~, ~, ~, ~, ~, ~, ~, somflin, ~, ~, hflin, ~, ~] = flprop(Runiv, pref, pdata(indx,13), fuel(17), Mfl, mofl(indx), fuelchem);
                dHindt(indx) = fuelinj(indx)*hflin;
                smflin = somflin - Runiv*log(fuel(18)/pref);      % [J/molK]
                sflin = smflin/Mfl;         % [J/kgK] = [J/molK]*[mol/kg]
                dSindt(indx) = fuelinj(indx)*sflin;     % [W] = [J/kgK]*[kg/s]
                % Exergy calculations
                % First is the exergy due to heat transfer
                dEHTfldt(indx) = dQHTfldt(indx)*(1 - To/Tflavg);        % [W]
                % Second is the specific exergy flow in at the pressure and
                % temperature of the standard state
                [~, ~, ~, ~, ~, ~, ~, somflref, ~, ~, hflref, ~, ~] = flprop(Runiv, pref, po, To, Mfl, mofl(indx), fuelchem);
                smflref = somflref - Runiv*log(po/pref);
                sflref = smflref/Mfl;
                % Need to add in the chemical exergy component
                % fuel(11) = fuel chemical exergy [J/mol]
                % fuel(5) = molecular mass [kg/mol]
                echfl = fuel(11)/fuel(5);       % [J/kg] = [J/mol] * [mol/kg]
                dEindt(indx) = fuelinj(indx)*((hflin-hflref) - To*(sflin-sflref) + echfl);  % [W]
                %--- Check for convergence
                epsit = abs(YO2CV - YCV(indx,1));
            end
            % Update indx
            indx = indx + 1;
        end
    else   % We do NOT have an added gas
        indx = IVC+1;
        % Run until the first injection process occurs
        while (indx < INJ1)     % No liquid fuel, no burned gases, no added gas fuel
            % Since there is no mass entering or exiting - unburned gases first
            mou_O2(indx) = mou_O2(indx-1);
            mou_N2(indx) = mou_N2(indx-1);
            mou_Ar(indx) = mou_Ar(indx-1);
            mou_CO2(indx) = mou_CO2(indx-1);
            mou_H2O(indx) = mou_H2O(indx-1);
            mu(indx,:) = mu(indx-1,:);
            % Added gas does not exist
            mofa(indx) = 0;
            mfa(indx) = 0;
            % Total control volume
            [mCV(indx), YCV, XCV, RCV(indx), TCV(indx), cpCV(indx), cvCV(indx), UCV(indx), HCV(indx), SCV(indx), ECV(indx)] = ...
                cvprop(indx, inp, Runiv, pref, pdata(indx,13), pdata(indx,5), YCV, XCV, MM, Mfa, Mfl, ...
                mu, mfa(indx), 0, mb, O2c, N2c, Arc, CO2c, H2Oc, fuel, fuelchem);
            % We only have unubrned mixture
            Tu(indx) = TCV(indx);
            % Check internal energy and ideal gas law
            % Unburned mixture
            [mu, Xu, Yu, Mu(indx), cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou, Ru(indx), Vu(indx)] = ...
                unburnprop(indx, Runiv, pref, pdata(indx,13), Tu(indx), MM, ...
                mou_O2(indx), mou_N2(indx), mou_Ar(indx), mou_CO2(indx), mou_H2O(indx), O2c, N2c, Arc, CO2c, H2Oc, ...
                mu, Xu, Yu, cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou);
            Uu = mu(indx,6)*uu(indx,6);
            % Checks - yeah, there's no real need to do this, but just
            % being consistent with the other routine
            g1 = UCV(indx) - Uu;
            g2 = mCV(indx)*RCV(indx)*TCV(indx) - mu(indx,6)*Ru(indx)*Tu(indx);
            indx = indx + 1;    % Move to the next step
        end
        % Now, we add the fuel until ignition happens
        while (indx < IGNfl)
            % Unburned gas will not change in mass
            mou_O2(indx) = mou_O2(indx-1);
            mou_N2(indx) = mou_N2(indx-1);
            mou_Ar(indx) = mou_Ar(indx-1);
            mou_CO2(indx) = mou_CO2(indx-1);
            mou_H2O(indx) = mou_H2O(indx-1);
            mu(indx,:) = mu(indx-1,:);
            % Added gas does not exist
            mofa(indx) = 0;
            mfa(indx) = 0;
            % Injection of the liquid fuel
            mfl(indx) = mfl(indx-1) + fuelinj(indx)*dt;     % Total mass of liquid fuel [kg]
            mofl(indx) = mfl(indx)/Mfl;                 % Total moles of liquid fuel [mol]
            % Total control volume
            [mCV(indx), YCV, XCV, RCV(indx), TCV(indx), cpCV(indx), cvCV(indx), UCV(indx), HCV(indx), SCV(indx), ECV(indx)] = ...
                cvprop(indx, inp, Runiv, pref, pdata(indx,13), pdata(indx,5), YCV, XCV, MM, Mfa, Mfl, ...
                mu, mfa(indx), mfl(indx), mb, O2c, N2c, Arc, CO2c, H2Oc, fuel, fuelchem);
            % To balance the ideal gas law and the internal energy,
            % the temperatures all need to be the same
            Tu(indx) = TCV(indx);
            Tfl(indx) = TCV(indx);
            % Check internal energy and ideal gas law
            % Unburned mixture
            [mu, Xu, Yu, Mu(indx), cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou, Ru(indx), Vu(indx)] = ...
                unburnprop(indx, Runiv, pref, pdata(indx,13), Tu(indx), MM, ...
                mou_O2(indx), mou_N2(indx), mou_Ar(indx), mou_CO2(indx), mou_H2O(indx), O2c, N2c, Arc, CO2c, H2Oc, ...
                mu, Xu, Yu, cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou);
            Uu = mu(indx,6)*uu(indx,6);
            % Liquid fuel
            [mfl(indx), Rfl, Vfl(indx), cpmfl(indx), cvmfl(indx), hmfl(indx), umfl(indx), somfl(indx), ...
                cpfl(indx), cvfl(indx), hfl(indx), ufl(indx), sofl(indx)] = flprop(Runiv, pref, pdata(indx,13), Tfl(indx), Mfl, mofl(indx), fuelchem);
            Ufl = mfl(indx)*ufl(indx);
            % Checks
            g1 = UCV(indx) - Uu - Ufl;
            g2 = mCV(indx)*RCV(indx)*TCV(indx) - mu(indx,6)*Ru(indx)*Tu(indx) - mfl(indx)*Rfl*Tfl(indx);
            % Heat transfer of liquid fuel
            % The liquid fuel heats up until it vaporizes
            dQHTfldt(indx) = fuelinj(indx)*(fuel(10)*(inp(25)-fuel(17)) - fuel(6) + cpfl(indx)*(fuel(17)-Tfl(indx)));
            % Entropy transfer due to heat transfer
            % Consider the system and surroundings. The system is our
            % gas that is losing heat to the liquid fuel subsequently
            % heating it up. Thus, the liquid fuel is our surroundings.
            % Now, what temperature do we use for the denominator. As
            % an example, the fuel is initially at 310 K, it heats up
            % to 341 K to vaporize, vaporizes at that level, and then
            % heats up to the temperature of the control volume at
            % around 900 K. So, it is not a constant temperature to
            % which this heat transfer occurs. I would suggest then we
            % use a weighted value based on the energy it takes for
            % each facet.
            % Qheat1 = cp*(310 K - 341 K)
            % Qheat2 = hfg @ 341 K
            % Qheat3 = cp,g*(341 K- ~900K)
            Tflavg = (0.5*(inp(25)+fuel(17))*fuel(10)*(inp(25)-fuel(17)) - fuel(17)*fuel(6) + ...
                0.5*(fuel(17)+Tfl(indx))*(cpfl(indx)*(fuel(17)-Tfl(indx))))/ ...
                (fuel(10)*(inp(25)-fuel(17)) - fuel(6) + cpfl(indx)*(fuel(17)-Tfl(indx)));
            dSHTfldt(indx) = dQHTfldt(indx)/Tflavg;     % [W/K]
            % Enthalpy addition into the control volume
            % The control volume contains the gas and not the liquid.
            % The liquid fuel mass technically is only going to be
            % added to the control volume once it vaporizes. Let's find
            % the enthalpy of the fuel that is being added at the fuel
            % volume temperature
            [~, ~, ~, ~, ~, ~, ~, somflin, ~, ~, hflin, ~, ~] = flprop(Runiv, pref, pdata(indx,13), fuel(17), Mfl, mofl(indx), fuelchem);
            dHindt(indx) = fuelinj(indx)*hflin;
            % OK, we have entropy coming in. It is coming in at the
            % fuel vaporization temperature (fuel(17)). In other words,
            % it only adds to our control volume (i.e., gas) when it is
            % a gas. As for its pressure? When it is injected, it is
            % injected at the nozzle pressure (inp(24)). However, it
            % expands into the control volume that will change its
            % pressure. Here, we assume it is at the vapor pressure
            smflin = somflin - Runiv*log(fuel(18)/pref);      % [J/molK]
            sflin = smflin/Mfl;         % [J/kgK] = [J/molK]*[mol/kg]
            dSindt(indx) = fuelinj(indx)*sflin;     % [W] = [J/kgK]*[kg/s]
            % Exergy calculations
            % First is the exergy due to heat transfer
            dEHTfldt(indx) = dQHTfldt(indx)*(1 - To/Tflavg);        % [W]
            % Second is the specific exergy flow in at the pressure and
            % temperature of the standard state
            [~, ~, ~, ~, ~, ~, ~, somflref, ~, ~, hflref, ~, ~] = flprop(Runiv, pref, po, To, Mfl, mofl(indx), fuelchem);
            smflref = somflref - Runiv*log(po/pref);
            sflref = smflref/Mfl;
            % Need to add in the chemical exergy component
            % fuel(11) = fuel chemical exergy [J/mol]
            % fuel(5) = molecular mass [kg/mol]
            echfl = fuel(11)/fuel(5);       % [J/kg] = [J/mol] * [mol/kg]
            dEindt(indx) = fuelinj(indx)*((hflin-hflref) - To*(sflin-sflref) + echfl);  % [W]
            indx = indx + 1;
        end
        % Ignition of the liquid fuel and run until EVO
        while (indx <= EVO)
            % added gas fuel
            mofa(indx) = 0;
            mfa(indx) = 0;
            % --- Liquid fuel --- %
            % Fuel could be added through fuel injection while burning at
            % the same time.
            dmfli = fuelinj(indx)*dt;       % Mass of liquid fuel being injected
            YflCV = YCV(indx-1,7);          % Mass fraction of liquid fuel in the cylinder
            dmflbdt = -Kfl*(rhoCV^2)*YflCV*(YO2CV^5)*pdata(indx-1,5)*exp(-Efl/(Runiv*TCV(indx-1))); % [kg/s]
            dmflb = dmflbdt*dt;             % Mass [kg] change due to combustion
            dmfl = dmfli + dmflb;           % Total mass of liquid fuel change [kg]
            % Now, we convert this to a molar change so that we can use the
            % local burn reaction
            dmofl = dmfl/Mfl;           % [mol]
            mofl(indx) = mofl(indx-1) + dmofl;
            % Make sure it does not go negative
            if (mofl(indx) < 0)
                mofl(indx) = 0;
                dmofl = mofl(indx) - mofl(indx-1);
            end
            % New mass of liquid fuel
            mfl(indx) = mofl(indx)*Mfl;
            % --- Total change of it all
            dmou_O2 = l2*dmofl;      % [mol] change in O2, unburned
            dmob_CO2 = -l3*dmofl;    % [mol] change in CO2, burned
            dmob_H2O = -l4*dmofl;    % [mol] change in H2O, burned
            dmob_N2 = -l5*dmofl;     % [mol] change in N2, burned
            % Note, we should never run out of oxygen since we are always
            % running lean for a CI engine. If we assume that the nitrogen
            % and argon are with oxygen in the air, as oxygen converts,
            % shouldn't we also move that nitrogen and argon with the
            % oxygen to the other side. Nitrogen and argon do not convert
            % to things like CO2 and water, but they should show up in
            % the products side. Same thing with the residual and EGR, if
            % that is embedded with the air, as oxygen burns, the residual
            % and EGR should find their way to the burned side. So, let's
            % find the change in moles of the other unburned species
            dmou_N2 = dmou_O2*(delta*g3 + epsilon*g8)/(delta*g2 + epsilon*g5);
            dmou_CO2 = dmou_O2*epsilon*g6/(delta*g2 + epsilon*g5);
            dmou_H2O = dmou_O2*epsilon*g7/(delta*g2 + epsilon*g5);
            dmou_Ar = dmou_O2*(delta*g4 + epsilon*g9)/(delta*g2 + epsilon*g5);
            % Update the arrays of the species
            % Unburned side - mol
            mou_O2(indx) = mou_O2(indx-1) + dmou_O2;
            mou_N2(indx) = mou_N2(indx-1) + dmou_N2;
            mou_Ar(indx) = mou_Ar(indx-1) + dmou_Ar;
            mou_CO2(indx) = mou_CO2(indx-1) + dmou_CO2;
            mou_H2O(indx) = mou_H2O(indx-1) + dmou_H2O;
            % Unburned side - kg
            mu(indx,1) = mou_O2(indx)*MO2;          % [mol] * [kg/mol] = [kg]
            mu(indx,2) = mou_N2(indx)*MN2;
            mu(indx,3) = mou_Ar(indx)*MAr;
            mu(indx,4) = mou_CO2(indx)*MCO2;
            mu(indx,5) = mou_H2O(indx)*MH2O;
            mu(indx,6) = mu(indx,1)+mu(indx,2)+mu(indx,3)+mu(indx,4)+mu(indx,5);
            % Burned side
            % Through my logic, the moles of nitrogen we lose from the
            % unburned side should also show up on the burned side
            mob_O2(indx) = 0;
            % We are gaining nitrogen on the burned side from that coming
            % from the fuel nitrogen and the moving of unburned N2 (air,
            % residual, EGR)dm
            mob_N2(indx) = mob_N2(indx-1) + dmob_N2 - dmou_N2;
            mob_Ar(indx) = mob_Ar(indx-1) - dmou_Ar;
            mob_CO2(indx) = mob_CO2(indx-1) + dmob_CO2 - dmou_CO2;
            mob_H2O(indx) = mob_H2O(indx-1) + dmob_H2O - dmou_H2O;
            % Burned side - kg
            mb(indx,1) = mob_O2(indx)*MO2;          % [mol] * [kg/mol] = [kg]
            mb(indx,2) = mob_N2(indx)*MN2;
            mb(indx,3) = mob_Ar(indx)*MAr;
            mb(indx,4) = mob_CO2(indx)*MCO2;
            mb(indx,5) = mob_H2O(indx)*MH2O;
            mb(indx,6) = mb(indx,1)+mb(indx,2)+mb(indx,3)+mb(indx,4)+mb(indx,5);
            % Total control volume
            [mCV(indx), YCV, XCV, RCV(indx), TCV(indx), cpCV(indx), cvCV(indx), UCV(indx), HCV(indx), SCV(indx), ECV(indx)] = ...
                cvprop(indx, inp, Runiv, pref, pdata(indx,13), pdata(indx,5), YCV, XCV, MM, Mfa, Mfl, ...
                mu, mfa(indx), mfl(indx), mb, O2c, N2c, Arc, CO2c, H2Oc, fuel, fuelchem);
            % To balance the ideal gas law and the internal energy,
            % the temperatures all need to be the same
            Tu(indx) = TCV(indx);
            Tfl(indx) = TCV(indx);
            Tb(indx) = TCV(indx);
            % Check internal energy and ideal gas law
            % Unburned mixture
            [mu, Xu, Yu, Mu(indx), cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou, Ru(indx), Vu(indx)] = ...
                unburnprop(indx, Runiv, pref, pdata(indx,13), Tu(indx), MM, ...
                mou_O2(indx), mou_N2(indx), mou_Ar(indx), mou_CO2(indx), mou_H2O(indx), O2c, N2c, Arc, CO2c, H2Oc, ...
                mu, Xu, Yu, cpmu, cvmu, hmu, umu, somu, cpu, cvu, hu, uu, sou);
            Uu = mu(indx,6)*uu(indx,6);
            % Liquid fuel
            [mfl(indx), Rfl, Vfl(indx), cpmfl(indx), cvmfl(indx), hmfl(indx), umfl(indx), somfl(indx), ...
                cpfl(indx), cvfl(indx), hfl(indx), ufl(indx), sofl(indx)] = flprop(Runiv, pref, pdata(indx,13), Tfl(indx), Mfl, mofl(indx), fuelchem);
            Ufl = mfl(indx)*ufl(indx);
            % Find the internal energy of the burned zone
            % We can use the same routine as the unburned properties
            [mb, Xb, Yb, Mb(indx), cpmb, cvmb, hmb, umb, somb, cpb, cvb, hb, ub, sob, Rb(indx), Vb(indx)] = ...
                unburnprop(indx, Runiv, pref, pdata(indx,13), Tb(indx), MM, ...
                mob_O2(indx), mob_N2(indx), mob_Ar(indx), mob_CO2(indx), mob_H2O(indx), O2c, N2c, Arc, CO2c, H2Oc, ...
                mb, Xb, Yb, cpmb, cvmb, hmb, umb, somb, cpb, cvb, hb, ub, sob);
            Ub = mb(indx,6)*ub(indx,6);
            % Checks
            g1 = UCV(indx) - Uu - Ufl - Ub;
            g2 = mCV(indx)*RCV(indx)*TCV(indx) - mu(indx,6)*Ru(indx)*Tu(indx) - mfl(indx)*Rfl*Tfl(indx) - mb(indx,6)*Rb(indx)*Tb(indx);
            % Heat transfer of liquid fuel
            dQHTfldt(indx) = fuelinj(indx)*(fuel(10)*(inp(25)-fuel(17)) - fuel(6) + cpfl(indx)*(fuel(17)-Tfl(indx)));
            % Entropy due to this heat transfer
            Tflavg = (0.5*(inp(25)+fuel(17))*fuel(10)*(inp(25)-fuel(17)) - fuel(17)*fuel(6) + ...
                0.5*(fuel(17)+Tfl(indx))*(cpfl(indx)*(fuel(17)-Tfl(indx))))/ ...
                (fuel(10)*(inp(25)-fuel(17)) - fuel(6) + cpfl(indx)*(fuel(17)-Tfl(indx)));
            dSHTfldt(indx) = dQHTfldt(indx)/Tflavg;     % [W/K]
            % Enthalpy addition into the control volume
            [~, ~, ~, ~, ~, ~, ~, somflin, ~, ~, hflin, ~, ~] = flprop(Runiv, pref, pdata(indx,13), fuel(17), Mfl, mofl(indx), fuelchem);
            dHindt(indx) = fuelinj(indx)*hflin;
            smflin = somflin - Runiv*log(fuel(18)/pref);      % [J/molK]
            sflin = smflin/Mfl;         % [J/kgK] = [J/molK]*[mol/kg]
            dSindt(indx) = fuelinj(indx)*sflin;     % [W] = [J/kgK]*[kg/s]
            % Exergy calculations
            % First is the exergy due to heat transfer
            dEHTfldt(indx) = dQHTfldt(indx)*(1 - To/Tflavg);        % [W]
            % Second is the specific exergy flow in at the pressure and
            % temperature of the standard state
            [~, ~, ~, ~, ~, ~, ~, somflref, ~, ~, hflref, ~, ~] = flprop(Runiv, pref, po, To, Mfl, mofl(indx), fuelchem);
            smflref = somflref - Runiv*log(po/pref);
            sflref = smflref/Mfl;
            % Need to add in the chemical exergy component
            % fuel(11) = fuel chemical exergy [J/mol]
            % fuel(5) = molecular mass [kg/mol]
            echfl = fuel(11)/fuel(5);       % [J/kg] = [J/mol] * [mol/kg]
            dEindt(indx) = fuelinj(indx)*((hflin-hflref) - To*(sflin-sflref) + echfl);  % [W]
            indx = indx + 1;
        end
    end

    % --- Calculate combustion efficiency ---%
    Ymb_fl = (1 - mfl(indx-1)/max(mfl));
    epsetac1 = abs(Ymb_fl - etac_flt);
    if (iS(2) == 1)
        Ymb_fa = (1 - mfa(indx-1)/mfa(IVC));
        epsetac2 = abs(Ymb_fa - etac_fat);
        if (itetac == 1)
            Kfaold = Kfa;
            Ymb_faold = Ymb_fa;
            Kfa = Kfa*1.0001;
            Kflold = Kfl;
            Ymb_flold = Ymb_fl;
            Kfl = Kfl*1.0001;
            itetac = 2;
        else
            % Newton Raphson
            Kfanew = Kfa - (Ymb_fa-etac_fat)/((Ymb_fa-Ymb_faold)/(Kfa - Kfaold))/2;
            % Check to make sure > 0
            if (Kfanew < 0)
                Kfanew = Kfa/(50*rand());
            end
            Kfaold = Kfa;
            Kfa = Kfanew;
            Ymb_faold = Ymb_fa;
            Kflnew = Kfl - (Ymb_fl-etac_flt)/((Ymb_fl-Ymb_flold)/(Kfl - Kflold))/2;
            if (Kflnew < 0)
                Kflnew = Kfl/(50*rand());
            end
            Kflold = Kfl;
            Kfl = Kflnew;
            Ymb_flold = Ymb_fl;
        end
        epsetac = max(epsetac1,epsetac2)
    else
        epsetac = 0;
    end
    Kfa
    Kfl
end
% Update the K values for the next time through to go faster
inp(73) = Kfa;
inp(75) = Kfl;
