%% Heat Release 2022
% Authors: Dr. Christopher Depcik, Dr. Jonathan Mattson, Dr. Shah Saud Alam
% Version: 2022a
% Date: 1/4/2023
clearvars;
% --- Data files --- %
pdatafile = '18.0NM_DME_ULSD_18ESRRaw.txt';                % Raw pressure data [bar]
mdatafile = 'YanmarL100vMotoringv2014.txt';                 % Raw motoring pressure data [bar]
bdatafile = 'SampleBiodiesel.txt';                          % Biodiesel GC data

% --- Switches --- %
iS(1) = 1;              % Number of injection events [-]
iS(2) = 1;              % Added fuel to intake: 0-No, 1-Yes (e.g., DME assisted biodiesel)
iS(3) = 1;              % First-order finite difference: 1-central difference O(2); 2-forward diff O(1); 3-backward diff O(1)
iS(4) = 1;              % Second-order finite difference: 1-central difference O(2); 2-forward diff O(1); 3-backward diff O(1)
iS(5) = 2;              % Liquid fuel: 1-biodiesel; 2-diesel
iS(6) = 1;              % Specify ignition timing: 0-No, 1-Yes
% When the ignition timing is not specified, it will be computed using the
% second derivative of the pressure trace. However, there are situations
% where it might be best to investigate this second derivative pressure trace
% and specify the ignition timing instead of computing it. Right now, there
% is the potential for two ignition points. One for the gaseous fuel and
% one for the liquid fuel.
iS(7) = 0;              % Remove injection artifact (0-No, 1-Yes)
% The model currently assumes that the fuel injected immediately vaporizes and does not go through
% an evaporation process. Running the model without showing this artifact
% can provide better insight into the heat release as the artifact can get
% in the way.
iS(8) = 0;              % Combustion model: 0-single reaction, 1-expanded

% --- Input --- %
% Why put everthing in an array? It makes it easier to add additional input
% and send it into subroutines. It also makes it easier to print the output
% Geometry
inp(1) = 0.086;                     % Engine bore [m]
inp(2) = 0.075;                     % Engine stroke [m]
inp(3) = 0.118;                     % Connecting rod length [m]
inp(4) = 0.000435;                  % Engine displacement [m3] - this could be computed from the bore and stroke
inp(5) = 21.2;                      % Compression ratio [-]
inp(6) = 0;                         % Piston bowl volume [m3]
inp(7) = 4;                         % Two(2) or four(4)-stroke cycle [strokes/cycle]
inp(8) = inp(2)/2;                  % Crank length [m]
inp(9) = inp(4)/(inp(5)-1);         % Clearance volume [m3]
inp(10) = pi()*inp(1)*inp(1)/4;     % Piston surface area [m2]
inp(11) = -122;                     % Intake valve closing BTDC [deg]
inp(12) = 144;                      % Exhaust valve opening ATDC [deg]
% Chemical species
% 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
inp(13) = 10;                       % Number of chemical species

% Operating conditions
% The program is written to be able to handle gas added in the intake and
% up to 5 injection events. It can be expanded to handle more injection
% events if need be
inp(20) = 1799.90256;               % Engine speed [rev/min]
inp(21) = 0.005814897;              % Mass flow rate of air [kg/s]
inp(22) = -9.6;                     % Injection timing (BTDC) of first injection event [deg]
inp(23) = 0.000193256;              % Mass flow rate of fuel added during first injection event [kg/s]
inp(24) = 47000000;                 % Injection pressure for first injection event [Pa]
inp(25) = 310;                      % Temperature of injected fuel [K]
% If one is to use multiple injection events, then one needs to consider
% what data has been taken. It seems logical that the fuel flow rate
inp(26) = 0;                        % Injection timing (BTDC) of second injection event [deg]
inp(27) = 0;                        % Mass flow rate of fuel added during second injection event [kg/s]
inp(28) = 0;                        % Injection pressure for second injection event [Pa]
inp(29) = 0;                        % Temperature of second injected fuel [K]
inp(30) = 0;                        % Injection timing (BTDC) of third injection event [deg]
inp(31) = 0;                        % Mass flow rate of fuel added during third injection event [kg/s]
inp(32) = 0;                        % Injection pressure for third injection event [Pa]
inp(33) = 0;                        % Temperature of third injected fuel [K]
inp(34) = 0;                        % Injection timing (BTDC) of fourth injection event [deg]
inp(35) = 0;                        % Mass flow rate of fuel added during fourth injection event [kg/s]
inp(36) = 0;                        % Injection pressure for fourth injection event [Pa]
inp(37) = 0;                        % Temperature of fourth injected fuel [K]
inp(38) = 0;                        % Injection timing (BTDC) of fifth injection event [deg]
inp(39) = 0;                        % Mass flow rate of fuel added during fifth injection event [kg/s]
inp(40) = 0;                        % Injection pressure for fifth injection event [Pa]
inp(41) = 0;                        % Temperature of fifth injected fuel [K]
inp(42) = 1;                        % Number of injectors [-]
inp(43) = 6;                        % Number of injector holes [-]
inp(44) = 0.00017;                  % Injector hole diameter [m]
inp(45) = pi()*inp(44)^2/4;         % Injector hole area [m2]
inp(46) = 0.000065004;              % Added gas mass flow rate [kg/s]
% In reality, the next two items simply indicate when to start computing
% the combustion reactions
inp(47) = -40;                      % Ignition point of added gas fuel (BTDC) [deg]
inp(48) = -40;                      % Ignition point of liquid fuel (BTDC) [deg]
% Ran without specifying the ignition timing; the second derivative tells
% you that ignition happens at -0.8 deg BTDC. Then, I reset to specify the
% ignition timing due to the added gas fuel
% Other input
inp(60) = 307.503;                  % Estimate of IVC temperature [K]
inp(61) = 60;                       % Number of pressure layers in raw file [-]
inp(62) = 0.2;                      % Degree resolution of pressure datafile [deg]
inp(63) = 0.0;                      % Exhaust Gas Recirculation volume percentage [-]
inp(64) = 0.979919;                 % Combustion efficiency of added gas [-]
inp(65) = 0.979919;                 % Combustion efficiency of liquid fuel [-]
inp(66) = 0.4;                      % Thermodynamic Loss Angle shift [deg]
inp(67) = 700;                      % Estimate of exhaust temperature [K]
inp(68) = 101325;                   % Intake pressure [Pa]
inp(69) = 101325;                   % Exhaust pressure [Pa]
inp(70) = 1/inp(5)*inp(69)/inp(68)*inp(60)/inp(67); % Residual fraction estimate [-]
inp(71) = 0.37;                     % Wall absorptivity [-]
inp(72) = 0.39;                     % Constant in fuel injection model [-]

if (iS(8) == 0) % Single reaction
    inp(73) = 2.076e3;                    % Kfa
    inp(74) = 22210;                    % Efa for fuel air burn rate [J/mol]
    inp(75) = 5.511e7;                    % Kfl
    inp(76) = 124473;                   % Efl for liquid fuel burn rate [J/mol]
    inp(87) = 1.64257E-01;
    inp(89) = 9.30162E-02;
    inp(88) = 1.00753E+00;
else
    % Expanded reaction set
    inp(73) = 1.4e8;                    % Kfa
    inp(74) = 66120;                    % Efa for fuel air burn rate [J/mol]
    inp(75) = 7.6e18;                    % Kfl
    inp(76) = 279291;                   % Efl for liquid fuel burn rate [J/mol]
    inp(87) = 0.16425;                     % Evaporation term for fuel injection array [-]
    inp(89) = 0.092000;                     % Power factor on oxygen for liquid injected fuel [-]
    inp(88) = 1.0000;                        % Power factor on oxygen for added gas fuel [-]
end
inp(77) = 1e-3;                     % Conversion criteria on combustion efficiency [-]
inp(78) = 1e-8;                     % Conversion criteria on energies [J]
inp(79) = 2;                        % Coefficient in front of heat transfer expression [-]
inp(80) = 3/4;                      % Coefficient on Reynolds number for heat transfer [-]
inp(81) = 1/3;                      % Coefficient on Prandtl number of heat transfer [-]
inp(82) = 400;                      % Wall temperature [K]
inp(83) = 100000;                   % Chemkin reference pressure [Pa]
inp(84) = 25+273.15;                % Exergy reference temperature To [K]
inp(85) = 101325;                   % Exergy reference pressure po [Pa]
inp(86) = 300;                      % Ambient temperature [K] - for Carnot efficiency calculation
% Constants, air, etc.
% Note: program was written without water and carbon dioxide in air - this
% can be added later
inp(100) = 8.31446261815324;    % Universal gas constant [J/(mol*K)]
inp(101) = 12.011/1000;         % Molecular mass of carbon [kg/mol]
inp(102) = 15.999/1000;         % Molecular mass of oxygen [kg/mol]
inp(103) = 14.0067/1000;        % Molecular mass of nitrogen [kg/mol]
inp(104) = 1.00784/1000;        % Molecular mass of hydrogen [kg/mol]
inp(105) = 39.948/1000;         % Molecular mass of argon [kg/mol]
inp(106) = 4.002602/1000;       % Molecular mass of helium [kg/mol]
inp(110) = 0.2095;              % Mole fraction of oxygen in air [-]
inp(111) = 0.7808;              % Mole fraction of nitrogen in air [-]
inp(112) = 0.0097;              % Mole fraction of argon in air [-]
inp(113) = inp(110)*2*inp(102) + inp(111)*2*inp(103) + inp(112)*inp(105);   % Molecular mass of air [kg/mol]
inp(114) = 5.670374419*10^-8;   % Stefan-Boltzman constant [W/(m2*K^4)]
% Added gas specification - mole fraction percentage
% More gases can be added onto the end of this array if desired
% Do not overwrite any gases
agas(1) = 1;                    % Dimethyl ether
agas(2) = 0;                    % Carbon dioxide
agas(3) = 0;                    % Carbon monoxide
agas(4) = 0;                    % Ethane
agas(5) = 0;                    % Ethylene
agas(6) = 0;                    % Helium
agas(7) = 0;                    % Hydrogen
agas(8) = 0;                    % Isobutane
agas(9) = 0;                    % Methane
agas(10) = 0;                   % Nitric oxide
agas(11) = 0;                   % Nitrogen
agas(12) = 0;                   % Nitrogen dioxide
agas(13) = 0;                   % Nitrous oxide
agas(14) = 0;                   % Oxygen
agas(15) = 0;                   % Ozone
agas(16) = 0;                   % Propane
agas(17) = 0;                   % Water
% Langness - 2016 CNG example
% agas(1) = 0; agas(9) = 0.87; agas(4) = 0.0510; agas(16) = 0.0150; agas(8) = 0.0029;
% agas(11) = 0.0560; agas(2) = 0.0051;

%% Some Calculations
if (inp(7) == 4)                % Four stroke
    calc(1) = 720;              % Total number of degrees
else
    calc(1) = 360;
end
calc(2) = calc(1)/inp(62);                      % Total number of datapoints in pressure data
calc(3) = round((inp(11)+180)/inp(62) + 1);     % Index of IVC in datafile
calc(4) = round((inp(12)+180)/inp(62) + 1);     % Index of EVO in datafile
calc(5) = round((inp(22)+180)/inp(62) + 1);     % Index of INJ1
calc(6) = round((inp(26)+180)/inp(62) + 1);     % Index of INJ2
calc(7) = round((inp(30)+180)/inp(62) + 1);     % Index of INJ3
calc(8) = round((inp(34)+180)/inp(62) + 1);     % Index of INJ4
calc(9) = round((inp(38)+180)/inp(62) + 1);     % Index of INJ5
calc(10) = round((0 + 180)/inp(62) + 1);        % Index of TDC
calc(11) = inp(20)*360/60;                      % Degrees per second [deg/s] = [rev/min]*[360 deg/rev]*[min/60 s]
calc(12) = inp(62)/calc(11);                    % Time-step of data [s] = [deg] / [deg/s]
calc(13) = inp(20)/inp(7)/60*2;                 % Cycles per second = [rev/min]*[cycle/stroke]*[min/60 s]*[2 stroke/rev]
calc(14) = inp(66)/inp(62);                     % How many index points to shift TDC due to thermodynamic loss angle [-]
% calc(15) => Pressure at IVC [Pa]
% calc(16) => Volume at IVC [m3]
% calc(17) => Index of maximum 2nd derivative of pressure [-]
% calc(18) => Degree of maximum 2nd derivative of pressure [deg]
% calc(19) => Ignition point index [-] gaseous fuel
% calc(20) => Ignition point [deg] gaseous fuel
% calc(21) => Ignition point index [-] liquid fuel
% calc(22) => Ignition point [deg] liquid fuel
calc(23) = 2*inp(2)*inp(20)/60;                 % Mean piston speed [m/s] <= 2 strokes per revolution
calc(24) = calc(11)*pi()/180;                   % Radians per second [rad/s] = [deg/s] * 3.14156 rad/180 deg
% calc(30) => Index of end of first injection process (in fuelinject.m)
% calc(31) => Degree at end of first injection process (in fuelinject.m)
% calc(32) => Index of end of second injection process (in fuelinject.m)
% calc(33) => Degree at end of second injection process (in fuelinject.m)
% calc(34) => Index of end of third injection process (in fuelinject.m)
% calc(35) => Degree at end of third injection process (in fuelinject.m)
% calc(36) => Index of end of fourth injection process (in fuelinject.m)
% calc(37) => Degree at end of fourth injection process (in fuelinject.m)
% calc(38) => Index of end of fifth injection process (in fuelinject.m)
% calc(39) => Degree at end of fith injection process (in fuelinject.m)
% calc(40) => Maximum pressure rise rate [Pa/s]

%% Load the pressure data and perform the filtering
[pdata, imep] = ptraceeval(iS, pdatafile, inp, calc);
% imep(1) => imep of smoothed data [Pa]
% m = inp(61) - Number of pressure layers
% imep(2 to m+1) => individual imep values of each pressure layer [Pa]
% imep(m+2) => mean imep,net of individual pressure layers [Pa]
% imep(m+3) => standard deviation of imep,net of individual pressure layers [Pa]
% imep(m+4) => Coefficient of Variance of imep,net - individual pressures [-]
% Find the pressure and volume at IVC
calc(15) = pdata(calc(3),13);       % Pressure at IVC [Pa]
calc(16) = pdata(calc(3),5);        % Volume at IVC [m3]
% Motoring trace
[mdata, mmep] = ptraceeval(iS, mdatafile, inp, calc);
% Shift the motoring trace to fit the pressure data
% This is done at the time of injection as the compression process should
% be the same between the two datasets providing there is no combustion
% happening. This might not be the case when there is a gas like DME
% combusting earlier. So, this might need to be revised.
prat = pdata(calc(5),13)/mdata(calc(5),13);
% Scale the motoring curves
mdata(:,10) = prat*mdata(:,10);
mdata(:,11) = prat*mdata(:,11);
mdata(:,12) = prat*mdata(:,12);
mdata(:,13) = prat*mdata(:,13);
mdata(:,14) = prat*mdata(:,14);

%% Find the maximum pressure rise rate [Pa/s]
% Over the number of data points in the pressure data
calc(40) = 0;  % Maximum pressure rise rate
for i=2:calc(2)
    prise = (pdata(i,13)-pdata(i-1,13))/calc(12);   % [Pa/s]
    if (prise > calc(40))
        calc(40) = prise;
    end
end

%% Compute or specify the ignition timing
if (iS(6) == 0) % Compute using the second derivative of pressure
    % Previous experience has stated it is more stable to find the maximum
    % of the second derivative of pressure and then work backwards to find
    % the inflection point. It will be assumed that if there is both added
    % gas and liquid fuel, they will ignite at the same time according to
    % the second derivative of pressure
    p2max = 0;
    for i=1:calc(2)
        if (pdata(i,34) > p2max)
            calc(17) = i;           % Index of maximum 2nd derivative of pressure
            p2max = pdata(i,34);
        end
    end
    calc(18) = (calc(17)-1)*inp(62) - 180;      % Degree of maximum 2nd derivative
    % Step backward to find the inflection point
    ichk = 0;
    j = calc(17);
    while (ichk == 0)
        j = j-1;
        if (pdata(j,34)) < 0
            ichk = 1;
            calc(19) = j;                           % Ignition point index [-] - gaseous fuel
            calc(20) = (calc(19)-1)*inp(62) - 180;  % Ignition point [deg] - gaseous fuel
            calc(21) = calc(19);                    % Ignition point index [-] - liquid fuel
            calc(22) = calc(20);                    % Ignition point [deg] - liquid fuel
        end
    end
else     % Specify ignition timing
    calc(17) = 0;
    calc(18) = 0;
    calc(20) = inp(47);
    calc(19) = (calc(20) + 180)/inp(62) + 1;
    calc(22) = inp(48);
    calc(21) = (calc(22) + 180)/inp(62) + 1;
end
%% Fuel Specification - returns both the liquid and assisted gas values
[fuel, fuelchem] = fuelspecs(bdatafile,iS,inp,agas);

% Iterate over residual fraction
epsres = 100;
while (epsres > 1e-7)
    %% Combustion reactions - does global and both local reactions when there is an added gas
    cr = combreact(iS, inp, calc, fuel);

    %% Fuel Injection Array
    [fuelinj, calc] = fuelinject(iS, inp, calc, fuel, cr, pdata);
    % calc(30) = index of end of first injection process
    % calc(31) = degree at end of first injection process

    %% Mass and Energy Calculations
    [inp, mu, mb, mfa, mfl, mCV, YCV, XCV, RCV, TCV, cpCV, cvCV, UCV, HCV, SCV, ECV, ...
        Tfl, cpfl, dQHTfldt, dHindt, dSHTfldt, dSindt, dEHTfldt, dEindt] = massenergy(iS, inp, calc, fuel, fuelchem, fuelinj, cr, pdata);

    % End residual fraction iteration
    TEVO = TCV(calc(4));    % Temperature of gas at EVO
    TIVC = TCV(calc(3));    % Temperature of gas at IVC
    % The residual fraction is computed using the temperature of the intake and
    % of the exhaust. These are not the temperatures at IVC and EVO as the
    % temperature of the intake is likely less than IVC and the temperature of
    % the exhaust is less than EVO. If we assume they both are off by the same
    % amount, we'll get something comparable. So, let us just use:
    resold = inp(70);
    res = (1/inp(5))*(inp(69)/inp(68))*TIVC/TEVO;
    epsres = abs(resold-res);
    inp(70) = res;
end

%% 1st Law Heat Release
HR1out = firstlawHR(iS, pdata, TCV, cpCV, RCV, XCV, YCV, UCV, dQHTfldt, dHindt, ...
    mu, mfa, mfl, mCV, inp, calc, cr, fuel, fuelchem);

%% 2nd Law Heat Release
HR2out = secondlawHR(pdata, TCV, SCV, XCV, YCV, cpCV, RCV, dSHTfldt, dSindt, HR1out, mu, mfa, mfl, mCV, inp, calc, fuel, fuelchem);

%% Exergy Heat Release
EXRGYout = exergyHR(inp, calc, pdata, XCV, YCV, TCV, ECV, HR1out, HR2out, mu, mfa, mfl, mCV, dEHTfldt, dEindt, fuel, fuelchem);

% Flip the output from row to column format for plotting
TCV=TCV';

%% Output parameters
EVO = calc(4);
eta1 = HR1out(EVO,20)/HR1out(EVO,17);       % First law fuel conversion efficiency
etac = 1 - inp(86)/max(TCV);                % Carnot efficiency
eta2 = eta1/etac;                           % Second law efficiency
etae = EXRGYout(EVO,28)/EXRGYout(EVO,34);   % Exergetic efficiency
eta1t = eta1/inp(65);                       % First law thermal efficiency - based on liquid fuel combustion efficiency
% Energy substitution rate
if (iS(2) == 1) % We have an added gas
    ESR = inp(46)*fuel(26)/(inp(46)*fuel(26) + (inp(23)+inp(27)+inp(31)+inp(35)+inp(39))*fuel(7));
    eta1t = eta1/(0.5*(inp(65)+inp(64)));     % First law thermal efficiency - based on average of liquid and gaseous fuel combustion efficiencies
end
% Gross values
Wig = HR1out(EVO,20);                   % gross indicated work [J]
imepg = Wig/inp(4);                     % gross imep = Wig/Vd [Pa]
tig = Wig/(2*pi()*(inp(7)/2));          % gross torque = Wig/(2*pi()*nR) [N-m]
Pig = tig*inp(20)*2*pi()/60;            % [J/cycle]*[rev/min]*[2*pi rad/rev]*[min/3600 s] = [W]
bsfcg = ((inp(23)+inp(27)+inp(31)+inp(35)+inp(39))+inp(46))/Pig*1000*3600*1000;        % bsfc = mdotf/Pig = [kg/(W*s)]*[1000 W/kW]*[3600 s/hr]*[1000 g/kg]
bsec = ((inp(23)+inp(27)+inp(31)+inp(35)+inp(39))*fuel(7) + inp(46)*fuel(26))/Pig*3.6;               
% [kg/s]*[J/kg] = [W]/[W] = [-]
% Units are [MJ/kW*hr] = 1000000/3600000, so we multiply by 3.6if (iS(2) == 1)
% Brake specific energy consumption

%% Plots for review
subplot(2,1,1);
plot(pdata(:,3),HR1out(:,16),pdata(:,3),HR1out(:,26));
xlim([-80 100]);
ylim([-20 50]);
xlabel('Crank Angle [deg]');
ylabel('Rate of Heat Release [J/deg]');
subplot(2,1,2);
plot(pdata(:,3),HR1out(:,17),pdata(:,3),HR1out(:,24));
xlim([-80 100]);
xlabel('Crank Angle [deg]');
ylabel('Cumulative Rate of Heat Release [J]');
%% Compute myfun
myfun = 0;
% -40 to 80 deg ATDC
for i=701:1301
    myfun = myfun + (HR1out(i,16)-HR1out(i,26))^2;
    %myfun = myfun + ((HR1out(i,17)-HR1out(i,24))^2)/100;
end
myfun = sqrt(myfun)
% plot(pdata(:,3),fuelinj);
% xlim([-80 50]);
% plot(pdata(:,3),HR1out(:,17),pdata(:,3),HR1out(:,24));
% xlim([inp(11) inp(12)]);
% xlabel('Crank Angle [deg]');
% ylabel('Cumulative Heat Release [J]');