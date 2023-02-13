function [fuel, fuelchem] = fuelspecs(bdata,iS,inp,agas)

% --- Routine specific variables
MC = inp(101);      % Molecular mass of carbon [kg/mol]
MO = inp(102);      % Molecular mass of hydrogen [kg/mol]
MN = inp(103);      % Molecular mass of nitrogen [kg/mol]
MH = inp(104);      % Molecular mass of hydrogen [kg/mol]
MAr = inp(105);     % Molecular mass of argon [kg/mol]
MHe = inp(106);     % Molecular mass of helium [kg/mol]

%% Direct injected fuel
% Choose the fuel and create the properties
if (iS(5) == 1)     % Biodiesel - uses input GC file
    % Get the CHEMKIN data curve fits for the different constituents
    cd 'CHEMKIN III Data';
    fidcpo=fopen('mh.txt','r');
    fidcpy=fopen('mo.txt','r');
    fidcpa=fopen('md.txt','r');
    fidcplm = fopen('mlaurm.txt','r');
    fidcpm2 = fopen('myrist.txt','r');
    fidcpm3 = fopen('mpalmit.txt','r');
    fidcpm4 = fopen('mpalmito.txt','r');
    fidcpm5 = fopen('mstearate.txt','r');
    fidcpm6 = fopen('moleate.txt','r');
    fidcpm7 = fopen('melaidate.txt','r');
    fidcpm8 = fopen('mlinole.txt','r');
    fidcpm9 = fopen('mlinolo.txt','r');
    fidcpm10 = fopen('malplin.txt','r');
    fidcpm11 = fopen('mgamlin.txt','r');
    fidcpm12 = fopen('marachid.txt','r');
    fidcpm13 = fopen('meico.txt','r');
    fidcpm14 = fopen('mbehen.txt','r');
    fdcpo=textscan(fidcpo,'%f');
    fdcpy=textscan(fidcpy,'%f');
    fdcpa=textscan(fidcpa,'%f');
    fdcplm=textscan(fidcplm,'%f');
    fdcpm2=textscan(fidcpm2,'%f');
    fdcpm3=textscan(fidcpm3,'%f');
    fdcpm4=textscan(fidcpm4,'%f');
    fdcpm5=textscan(fidcpm5,'%f');
    fdcpm6=textscan(fidcpm6,'%f');
    fdcpm7=textscan(fidcpm7,'%f');
    fdcpm8=textscan(fidcpm8,'%f');
    fdcpm9=textscan(fidcpm9,'%f');
    fdcpm10=textscan(fidcpm10,'%f');
    fdcpm11=textscan(fidcpm11,'%f');
    fdcpm12=textscan(fidcpm12,'%f');
    fdcpm13=textscan(fidcpm13,'%f');
    fdcpm14=textscan(fidcpm14,'%f');
    caproate=fdcpo{1};
    caprylate=fdcpy{1};
    caprate=fdcpa{1};
    mlaurm=fdcplm{1};
    myrist=fdcpm2{1};
    mpalmit=fdcpm3{1};
    mpalmito=fdcpm4{1};
    mstearate=fdcpm5{1};
    moleate=fdcpm6{1};
    melaidate=fdcpm7{1};
    mlinole=fdcpm8{1};
    mlinolo=fdcpm9{1};
    malplin=fdcpm10{1};
    mgamlin=fdcpm11{1};
    marachid=fdcpm12{1};
    meico=fdcpm13{1};
    mbehen=fdcpm14{1};
    fclose(fidcpo);
    fclose(fidcpy);
    fclose(fidcpa);
    fclose(fidcplm);
    fclose(fidcpm2);
    fclose(fidcpm3);
    fclose(fidcpm4);
    fclose(fidcpm5);
    fclose(fidcpm6);
    fclose(fidcpm7);
    fclose(fidcpm8);
    fclose(fidcpm9);
    fclose(fidcpm10);
    fclose(fidcpm11);
    fclose(fidcpm12);
    fclose(fidcpm13);
    fclose(fidcpm14);
    cd ..
    % --- Set methyl ester species molar masses [kg/mol] and estimated heat
    % of formation
    M(1)=(MC*7+14*MH+2*MO);            %methyl caproate C6:0
    M(2)=(MC*9+18*MH+2*MO);            %methyl caprylate C8:0
    M(3)=(MC*11+22*MH+2*MO);           %methyl caprate C10:0
    M(4)=(MC*13+26*MH+2*MO);           %methyl laurate C12:0
    %hmof(4) = 0.5*(-714.6-693.6)*1000;  % https://webbook.nist.gov/cgi/cbook.cgi?ID=C111820&Mask=2#Thermo-Condensed
    M(5)=(MC*15+30*MH+2*MO);           %methyl myristate C14:0
    M(6)=(MC*17+34*MH+2*MO);           %methyl palmitate C16:0
    M(7)=(MC*17+32*MH+2*MO);           %methyl palmitoleate C16:1
    M(8)=(MC*19+38*MH+2*MO);           %methyl stearate C18:0
    M(9)=(MC*19+36*MH+2*MO);           %methyl oleate C18:1
    M(10)=(MC*19+36*MH+2*MO);          %methyl elaidate C18:1
    M(11)=(MC*19+34*MH+2*MO);          %methyl linoleate C18:2
    M(12)=(MC*19+34*MH+2*MO);          %methyl linloelaidate C18:2
    M(13)=(MC*19+32*MH+2*MO);          %methyl alpha-linolenate C18:3
    M(14)=(MC*19+32*MH+2*MO);          %methyl gamma-linolenate C18:3
    M(15)=(MC*21+42*MH+2*MO);          %methyl arachidate C20:0
    M(16)=(MC*21+40*MH+2*MO);          %methyl eicosenoate C20:1
    M(17)=(MC*23+46*MH+2*MO);          %methyl behenate C22:0
    % Read in gas-chromatography file
    gcfile=sprintf(bdata);
    fidgc=fopen(gcfile,'r');
    fgc=textscan(fidgc,'%s%f');
    Cnum=fgc{2};
    fclose(fidgc);
    % Calculate molecular mass of biodiesel mixture
    for i=1:17
        Cweight(i)=Cnum(i)/M(i);
    end
    Mbio=1/sum(Cweight);
    % Calculate the mole fraction of the species
    for i=1:17
        Cmol(i)=Cnum(i)*Mbio/M(i);
    end
    % Calculate the equivalent formula
    fuel(1) = Cmol(1)*7+Cmol(2)*9+Cmol(3)*11+Cmol(4)*13+Cmol(5)*15+Cmol(6)*17+Cmol(7)*17+...
        Cmol(8)*19+Cmol(9)*19+Cmol(10)*19+Cmol(11)*19+Cmol(12)*19+Cmol(13)*19+Cmol(14)*19 + ...
        Cmol(15)*21+Cmol(16)*21+Cmol(17)*23;
    fuel(2) = Cmol(1)*14+Cmol(2)*18+Cmol(3)*22+Cmol(4)*26+Cmol(5)*30+Cmol(6)*34+Cmol(7)*32 + ...
        Cmol(8)*38+Cmol(9)*36+Cmol(10)*36+Cmol(11)*34+Cmol(12)*34+Cmol(13)*32+Cmol(14)*32 + ...
        Cmol(15)*42+Cmol(16)*40+Cmol(17)*46;
    fuel(3) = 2;
    fuel(4) = 0;
    % Rest of the fuel properties
    % We are just approximating the overall biodiesel here. We could get
    % fancier and figure out mixture values for them
    fuel(5) = Mbio;         % Molecular mass [kg/mol]
    fuel(6) = 356.25*1000;  % Heat of vaporization [J/kg] - https://www.researchgate.net/publication/267492494_Combustion_Emissions_Modeling_and_Testing_of_Neat_Biodiesel_Fuels/figures?lo=1
    fuel(7) = 39798*1000;   % Lower Heating Value [J/kg]
    fuel(8) = 884;          % Density of biodiesel at ambient conditions [kg/m3]
    fuel(9) = 422.15;       % Not used: Flash point of biodiesel [K]
    fuel(10) = 1200;        % Liquid specific heat capacity [J/kgK]
    % https://www.sciencedirect.com/science/article/pii/S0961953416303609
    % Mass fractional components of biodiesel
    yc = fuel(1)*MC/(fuel(1)*MC + fuel(2)*MH + fuel(3)*MO + fuel(4)*MN);
    yh = fuel(2)*MH/(fuel(1)*MC + fuel(2)*MH + fuel(3)*MO + fuel(4)*MN);
    yo = fuel(3)*MO/(fuel(1)*MC + fuel(2)*MH + fuel(3)*MO + fuel(4)*MN);
    yn = fuel(4)*MN/(fuel(1)*MC + fuel(2)*MH + fuel(3)*MO + fuel(4)*MN);
    ys = 0;
    r = 1.0401 + 0.1728*yh/yc + 0.0432*yo/yc + 0.2169*ys/yc*(1-2.0628*yh/yc);
    fuel(11) = r*fuel(7)*Mbio;       % Chemical exergy [J/mol] = [J/kg]*[kg/mol]
    % https://www.researchgate.net/publication/267492494_Combustion_Emissions_Modeling_and_Testing_of_Neat_Biodiesel_Fuels/figures?lo=1
    % From the lower heating value, molecular mass, and the combustion
    % reaction, I calculated the heat of formation of biodiesel to be
    fuel(12) = (fuel(7)*fuel(5) + fuel(1)*-393522 + fuel(2)/2*(-285830))/fuel(5); % [J/kg]
    fuel(13) = 40*1e-6;     % Not used: Diameter of fuel droplets [m]
    fuel(14) = 4*1e-6;      % Not used: Kinematic viscosity of fuel [m2/s]
    fuel(15) = 18;          % Not used: Velocity of fuel droplets [m/s]
    fuel(16) = 0.14;        % Not used: Thermal conductivity of fuel [W/(mK)]
    fuel(17) = 341;         % Vaporization temperature [K] - https://www.mdpi.com/1996-1073/13/14/3637
    fuel(18) = 5.91*1000;   % Vapor pressure [Pa] - https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=910403
    % Calculate lumped CHEMKIN III coefficients for biodiesel fuel
    % There are fourteen total coefficients. Each coefficient becomes a
    % mole fraction average
    for i=1:14
        fuelchem(i)=Cmol(1)*caproate(i)+Cmol(2)*caprylate(i)+Cmol(3)*caprate(i)+Cmol(4)*mlaurm(i) + ...
            Cmol(5)*myrist(i)+Cmol(6)*mpalmit(i)+Cmol(7)*mpalmito(i)+Cmol(8)*mstearate(i) + ...
            Cmol(9)*moleate(i)+Cmol(10)*melaidate(i)+Cmol(11)*mlinole(i)+Cmol(12)*mlinolo(i) + ...
            Cmol(13)*malplin(i)+Cmol(14)*mgamlin(i)+Cmol(15)*marachid(i)+Cmol(16)*meico(i)+Cmol(17)*mbehen(i);
    end
elseif (iS(5) == 2) % Diesel fuel #1
    cd 'CHEMKIN III Data';
    fidd1=fopen('C14.4H24.9.txt','r');
    fd1=textscan(fidd1,'%f');
    diesel1=fd1{1};
    fclose(fidd1);
    cd ..
    % Equivalent formula
    fuel(1) = 14.4;     % C
    fuel(2) = 24.9;     % H
    fuel(3) = 0;        % O
    fuel(4) = 0;        % N
    fuel(5) = fuel(1)*MC + fuel(2)*MH + fuel(3)*MO + fuel(4)*MN; % Molecular mass [kg/mol]
    % https://www.sciencedirect.com/science/article/pii/S0360544207000527
    fuel(6) = 250*1000;         %heat of vaporization (J/kg)
    fuel(7) = 43000*1000;                   % Lower Heating Value [J/kg]
    fuel(8) = 837;              % Density of diesel at ambient conditions [kg/m3]
    fuel(9) = 335;              % Flash point of diesel [K]
    fuel(10) = 1850;            % Liquid specific heat capacity [J/kgK]
    fuel(11) = 1.065*fuel(7)*fuel(5);  % Chemical exergy [J/mol]
    fuel(12) = (fuel(7)*fuel(5) + fuel(1)*-393522 + fuel(2)/2*(-285830))/fuel(5); % [J/kg] - heat of formation
    fuel(13) = 0;                   % Not used: Diameter of fuel droplets [m]
    fuel(14) = 2.6/1000/1000;       % Kinematic viscosity of fuel [m2/s]
    fuel(15) = 0;                   % Not used: Velocity of fuel droplets [m/s]
    fuel(16) = 0;               % Not used: Thermal conductivity of fuel [W/(mK)]
    fuel(17) = 341;             % Vaporization temperature [K] - https://www.mdpi.com/1996-1073/13/14/3637
    fuel(18) = 0.40*133.322;    % Vapor pressure [Pa] - https://www.atlasoil.com/media/documents/safety-data-sheets/Conoco-Phillips/No-2_Diesel-Fuel.pdf
    for i=1:14
        fuelchem(i) = diesel1(i);
    end
end

%% Gaseous Added Fuel
if (iS(2) == 1)
    %Get the CHEMKIN data curve fits for the different constituents
    cd 'CHEMKIN III Data';
    fidDME = fopen('CH3OCH3.txt','r');
    fidCO2=fopen('CO2.txt','r');
    fidCO=fopen('CO.txt','r');
    fidC2H6=fopen('C2H6.txt','r');
    fidC2H4=fopen('C2H4.txt','r');
    fidHe=fopen('HE.txt','r');
    fidH2=fopen('H2.txt','r');
    fidC4H10=fopen('C4H10.txt','r');
    fidCH4=fopen('CH4.txt','r');
    fidNO=fopen('NO.txt','r');
    fidN2=fopen('N2.txt','r');
    fidNO2=fopen('NO2.txt','r');
    fidN2O=fopen('N2O.txt','r');
    fidO2=fopen('O2.txt','r');
    fidO3=fopen('O3.txt','r');
    fidC3H8=fopen('C3H8.txt','r');
    fidH2O=fopen('H2O.txt','r');
    fDME=textscan(fidDME,'%f');
    fO2=textscan(fidO2,'%f');
    fO3=textscan(fidO3,'%f');
    fCO2=textscan(fidCO2,'%f');
    fCO=textscan(fidCO,'%f');
    fHe=textscan(fidHe,'%f');
    fH2=textscan(fidH2,'%f');
    fH2O=textscan(fidH2O,'%f');
    fN2=textscan(fidN2,'%f');
    fNO=textscan(fidNO,'%f');
    fNO2=textscan(fidNO2,'%f');
    fN2O=textscan(fidN2O,'%f');
    fCH4=textscan(fidCH4,'%f');
    fC2H4=textscan(fidC2H4,'%f');
    fC2H6=textscan(fidC2H6,'%f');
    fC3H8=textscan(fidC3H8,'%f');
    fC4H10=textscan(fidC4H10,'%f');
    DMEc=fDME{1};
    O2c=fO2{1};
    O3c=fO3{1};
    CO2c=fCO2{1};
    COc=fCO{1};
    Hec=fHe{1};
    H2c=fH2{1};
    H2Oc=fH2O{1};
    N2c=fN2{1};
    NOc=fNO{1};
    NO2c=fNO2{1};
    N2Oc=fN2O{1};
    CH4c=fCH4{1};
    C2H4c=fC2H4{1};
    C2H6c=fC2H6{1};
    C3H8c=fC3H8{1};
    C4H10c=fC4H10{1};
    fclose(fidDME);
    fclose(fidO2);
    fclose(fidO3);
    fclose(fidCO2);
    fclose(fidCO);
    fclose(fidHe);
    fclose(fidH2);
    fclose(fidH2O);
    fclose(fidN2);
    fclose(fidNO);
    fclose(fidNO2);
    fclose(fidN2O);
    fclose(fidCH4);
    fclose(fidC2H4);
    fclose(fidC2H6);
    fclose(fidC3H8);
    fclose(fidC4H10);
    cd ..
    % --- Set the molar masses [kg/mol]
    Ma(1) = 2*MC + 6*MH + MO;
    Ma(2) = MC + 2*MO;
    Ma(3) = MC + MO;
    Ma(4) = 2*MC + 6*MH;
    Ma(5) = 2*MC + 4*MH;
    Ma(6) = MHe;
    Ma(7) = 2*MH;
    Ma(8) = 4*MC + 10*MH;
    Ma(9) = MC + 4*MH;
    Ma(10) = MN + MO;
    Ma(11) = 2*MN;
    Ma(12) = MN + 2*MO;
    Ma(13) = 2*MN + MO;
    Ma(14) = 2*MO;
    Ma(15) = 3*MO;
    Ma(16) = 3*MC + 8*MH;
    Ma(17) = 2*MH + MO;
    % --- Set the lower heating values [J/kg]
    % https://link.springer.com/content/pdf/bbm%3A978-981-287-212-8%2F1.pdf
    % https://link.springer.com/content/pdf/bbm%3A978-1-4419-7943-8%2F1.pdf
    Qlhva(1) = 28882*1000;      % Dimethyl ether  % C2H6O - https://doi.org/10.1016/j.jpowsour.2005.05.082
    emch(1) = 30.75*10^6*46.07/1000;  % Molar chemical exergy [J/mol] = [MJ/kg] * [10^6 J/MJ] * [kg/kmol] * [kmol/1000 mol]
    Qlhva(2)=0;                 % Carbon dioxide: CO2
    emch(2) = 19.48*1000;
    Qlhva(3)=10112*1000;        % Carbon monoxide: CO
    emch(3) = 274.71*1000;
    Qlhva(4)=47611*1000;        % Ethane: C2H6
    emch(4) = 1495.0*1000;
    Qlhva(5)=47132*1000;        % Ethylene: C2H4
    emch(5) = 1360.3*1000;
    Qlhva(6)=0;                 % Helium: He
    emch(6) = 30.37*1000;     % He - https://www.exergoecology.com/exergoecology/szargut2005
    Qlhva(7)=119960*1000;       % Hydrogen: H2
    emch(7) = 236.09*1000;
    Qlhva(8)=45577*1000;        % Isobutane: C4H1O
    emch(8) = 2804.2*1000;
    Qlhva(9)=50048*1000;        % Methane: CH4
    emch(9) = 831.2*1000;
    Qlhva(10)=0;                % Nitric Oxide: NO
    emch(10) = 88.9*1000;
    Qlhva(11)=0;                % Nitrogen: N2
    emch(11) = 0.72*1000;
    Qlhva(12)=0;                % Nitrogen Dioxide: NO2
    emch(12) = 55.6*1000;
    Qlhva(13)=0;                % Nitrous Oxide: N2O
    emch(13) = 106.9*1000;
    Qlhva(14)=0;                % Oxygen: O2
    emch(14) = 3.97*1000;
    Qlhva(15)=0;                % Ozone: O3
    emch(15) = 168.1*1000;
    Qlhva(16)=46330*1000;       % Propane: C3H8
    emch(16) = 2152.8*1000;
    Qlhva(17)=0;                % Water: H2O
    emch(17) = 9.5*1000; % H2O (gas) value is 0.9 when it is a liquid
    % Calculate molecular mass of mixture
    % Calculate chemical exergy of mixture
    Mmixa = 0;
    emchmixa = 0;
    for i=1:17
        Mmixa = Mmixa + agas(i)*Ma(i);
        emchmixa = emchmixa + agas(i)*emch(i);
    end
    % Calculate the equivalent formula
    fuel(20) = agas(1)*2 + agas(2)*1 + agas(3)*1 + agas(4)*2 + agas(5)*2 + ...
        agas(6)*0 + agas(7)*0 + agas(8)*4 + agas(9)*1 + agas(10)*0 + agas(11)*0 + ...
        agas(12)*0 + agas(13)*0 + agas(14)*0 + agas(15)*0 + agas(16)*3 + agas(17)*0;
    fuel(21) = agas(1)*6 + agas(2)*0 + agas(3)*0 + agas(4)*6 + agas(5)*4 + ...
        agas(6)*0 + agas(7)*2 + agas(8)*10 + agas(9)*4 + agas(10)*0 + agas(11)*0 + ...
        agas(12)*0 + agas(13)*0 + agas(14)*0 + agas(15)*0 + agas(16)*8 + agas(17)*2;
    fuel(22) = agas(1)*1 + agas(2)*2 + agas(3)*1 + agas(4)*0 + agas(5)*0 + ...
        agas(6)*0 + agas(7)*0 + agas(8)*0 + agas(9)*0 + agas(10)*1 + agas(11)*0 + ...
        agas(12)*2 + agas(13)*1 + agas(14)*2 + agas(15)*3 + agas(16)*0 + agas(17)*1;
    fuel(23) = agas(1)*0 + agas(2)*0 + agas(3)*0 + agas(4)*0 + agas(5)*0 + ...
        agas(6)*0 + agas(7)*0 + agas(8)*0 + agas(9)*0 + agas(10)*1 + agas(11)*2 + ...
        agas(12)*1 + agas(13)*3 + agas(14)*0 + agas(15)*0 + agas(16)*0 + agas(17)*0;
    % Rest of the added gas fuel properties
    fuel(24) = Mmixa;       % Molecular mass [kg/mol]
    fuel(25) = 0;           % heat of vaporization [J/mol] - not used
    Qlhvmmix = 0;
    for i=1:17
        Qlhvmmix = Qlhvmmix + agas(i)*Qlhva(i)*Ma(i);
    end
    fuel(26) = Qlhvmmix/Mmixa;  % Lower heating value [J/kg]
    fuel(27) = 0;           % Density at ambient conditions [kg/m3] - not used
    fuel(28) = 0;           % Flash point of added gas [K] - not used
    fuel(29) = 0;           % Liquid specific heat capacity [J/kgK] - not used
    fuel(30) = emchmixa;  % Molar chemical exergy [J/mol]
    % C2H6O - https://doi.org/10.1016/j.jpowsour.2005.05.082;

    % Calculate lumped CHEMKIN III coefficients for biodiesel fuel
    % There are fourteen total coefficients. Each coefficient becomes a
    % mole fraction average
    for i=1:14
        fuelchem(19+i) = agas(1)*DMEc(i) + agas(2)*CO2c(i) + agas(3)*COc(i) + ...
            agas(4)*C2H6c(i) + agas(5)*C2H4c(i) + agas(6)*Hec(i) + agas(7)*H2c(i) + ...
            agas(8)*C4H10c(i) + agas(9)*CH4c(i) + agas(10)*NOc(i) + agas(11)*N2c(i) + ...
            agas(12)*NO2c(i) + agas(13)*N2Oc(i) + agas(14)*O2c(i) + agas(15)*O3c(i) + ...
            agas(16)*C3H8c(i) + agas(17)*H2Oc(i);
    end
end

%% Molar chemical exergy of the species [kJ/mol]
%http://web.mit.edu/2.813/www/readings/APPENDIX.pdf
% exm(1) = 1416.6525;    % C2H6O - https://doi.org/10.1016/j.jpowsour.2005.05.082
% exm(2) = 19.48;         % CO2
% exm(3) = 274.71;        % CO
% exm(4) = 1495.0;        % C2H6
% exm(5) = 1360.3;        % C2H4
% exm(6) = 30.37;         % He - https://www.exergoecology.com/exergoecology/szargut2005
% exm(7) = 236.09;        % H2
% exm(8) = 2804.2;        % C4H10
% exm(9) = 831.2;         % CH4
% exm(10) = 88.9;          % NO
% exm(11) = 0.72;         % N2
% exm(12) = 55.6;         % NO2
% exm(13) = 106.9;        % N2O
% exm(14) = 3.97;         % O2
% exm(15) = 168.1;       % O3
% exm(16) = 2152.8;       % C3H8
% exm(17) = 9.5;          % H2O (gas) value is 0.9 when it is a liquid
% exm(18) = 12;         % Argon

%% CHEMKIN Biodiesel fits
% mh         1/24/ 7 thermc   7h  14o   2    0g   300.000  5000.000 1383.000    71
%  2.53098086E+01 3.24987915E-02-1.10223928E-05 1.69923778E-09-9.81167630E-14    2
% -7.15533311E+04-1.03851509E+02 2.54011733E+00 8.11387825E-02-4.96793859E-05    3
%  1.53233872E-08-1.90205405E-12-6.31288049E+04 2.01562116E+01                   4
% mo         1/24/ 7 thermc   9h  18o   2    0g   300.000  5000.000 1382.000    91
%  3.23164209E+01 4.06680764E-02-1.38575381E-05 2.14498069E-09-1.24215444E-13    2
% -8.00487771E+04-1.38892250E+02 2.15464837E+00 1.05249188E-01-6.52504339E-05    3
%  2.02445809E-08-2.51113060E-12-6.89196322E+04 2.52953802E+01                   4
% md         1/24/ 7 thermc  11h  22o   2    0g   300.000  5000.000 1382.000    11
%  3.93230373E+01 4.88368389E-02-1.66923510E-05 2.59065840E-09-1.50309877E-13    2
% -8.85441006E+04-1.73932688E+02 1.76901386E+00 1.29360919E-01-8.08243357E-05    3
%  2.51676921E-08-3.12062272E-12-7.47104475E+04 3.04352079E+01                   4
% 	! FAMEs (Fatty Acid Methyl Esters) (see sections 1, 2, 3, 3.1.1., 3.1.2., 3.2.1., 3.2.2.)
% https://pubs.acs.org/doi/abs/10.1021/jp904896r
% Species c13h26o2_laur_m_est Fit:
% c13h26o2_laur_m_        H  26C  13O   2     G   300.000  5000.000 1000.00      1
%  2.72501751E+01 8.47694402E-02-3.40659424E-05 6.26244016E-09-4.32293588E-13    2
% -8.91162999E+04-1.04824931E+02 9.57485028E+00 8.29372061E-02 7.74827088E-05    3
% -1.40636861E-07 5.44259149E-11-8.19545582E+04-1.41832275E+00                   4
% Species c15h30o2_myrist_m_est Fit:
% c15h30o2_myrist_        H  30C  15O   2     G   300.000  5000.000 1000.00      1
%  3.16691246E+01 9.70485101E-02-3.90176912E-05 7.17503271E-09-4.95405495E-13    2
% -9.61391160E+04-1.26743816E+02 1.06291320E+01 9.52713616E-02 9.25537095E-05    3
% -1.66476353E-07 6.44017206E-11-8.76342616E+04-3.75373409E+00                   4
% Species c17h32o2_palmitole_m_est Fit:
% c17h32o2_palmito        H  32C  17O   2     G   300.000  5000.000 1000.00      1
%  3.33606142E+01 1.06719124E-01-4.29500189E-05 7.90457654E-09-5.46125152E-13    2
% -8.68236054E+04-1.26816968E+02 1.14310960E+01 1.02101194E-01 1.02480878E-04    3
% -1.81385357E-07 6.98603588E-11-7.78209017E+04 2.06428096E+00                   4
% Species c17h34o2_palmit_m_est Fit:
% c17h34o2_palmit_        H  34C  17O   2     G   300.000  5000.000 1000.00      1
%  3.55305296E+01 1.09746228E-01-4.40981772E-05 8.10571172E-09-5.59469785E-13    2
% -1.02923616E+05-1.45292419E+02 1.17629996E+01 1.08027697E-01 1.03662595E-04    3
% -1.87190120E-07 7.24616510E-11-9.33306777E+04-6.43066315E+00                   4
% Species c18h34o2_margarole_m_est Fit:
% c18h34o2_margaro        H  34C  18O   2     G   300.000  5000.000 1000.00      1
%  3.50485489E+01 1.13308344E-01-4.55923811E-05 8.38989274E-09-5.79617589E-13    2
% -9.01480161E+04-1.32752863E+02 1.00539604E+01 1.19793418E-01 8.49199298E-05    3
% -1.72111595E-07 6.79190745E-11-8.04744343E+04 1.12048989E+01                   4
% Species c18h36o2_margar_m_est Fit:
% c18h36o2_margar_        H  36C  18O   2     G   300.000  5000.000 1000.00      1
%  3.72861457E+01 1.16329013E-01-4.67604775E-05 8.59835934E-09-5.93685219E-13    2
% -1.06231428E+05-1.52474293E+02 8.19916026E+00 1.38983265E-01 5.97986793E-05    3
% -1.56134768E-07 6.40130192E-11-9.57293460E+04 1.12770202E+01                   4
% Species c19h32o2_alp_linolen_m_est Fit:
% c19h32o2_alp_lin        H  32C  19O   2     G   300.000  5000.000 1000.00      1
%  3.58895203E+01 1.10175658E-01-4.43745442E-05 8.17216614E-09-5.64924900E-13    2
% -6.16590179E+04-1.38785637E+02 1.03444238E+01 1.18705837E-01 8.33054991E-05    3
% -1.70598070E-07 6.75401860E-11-5.18674882E+04 7.86723855E+00                   4
% Species c19h34o2_linole_m_est Fit:
% c19h34o2_linole_        H  34C  19O   2     G   300.000  5000.000 1000.00      1
%  3.65096642E+01 1.14873436E-01-4.62435721E-05 8.51265063E-09-5.88247571E-13    2
% -7.75019016E+04-1.41570033E+02 1.13703260E+01 1.16880160E-01 9.85722841E-05    3
% -1.86581882E-07 7.28230427E-11-6.75465025E+04 4.35039881E+00                   4
% Species c19h36o2_ole_m_est Fit:
% c19h36o2_ole_m_e        H  36C  19O   2     G   300.000  5000.000 1000.00      1
%  3.65188489E+01 1.20153901E-01-4.83472752E-05 8.89684509E-09-6.14639562E-13    2
% -9.33691225E+04-1.37669517E+02 1.14337130E+01 1.21980228E-01 9.66845577E-05    3
% -1.86305260E-07 7.28144407E-11-8.34263845E+04 7.98031762E+00                   4
% Species c19h38o2_stear_m_est Fit:
% c19h38o2_stear_m        H  38C  19O   2     G   300.000  5000.000 1000.00      1
%  3.96632507E+01 1.22304300E-01-4.91657272E-05 9.04051604E-09-6.24181172E-13    2
% -1.09835419E+05-1.65539291E+02 1.16440324E+01 1.28610900E-01 1.00029782E-04    3
% -1.96193430E-07 7.71268737E-11-9.89430613E+04-3.92019089E+00                   4
% Species c21h42o2_arachid_m_est Fit:
% c21h42o2_arachid        H  42C  21O   2     G   300.000  5000.000 1000.00      1
%  4.38052601E+01 1.34942383E-01-5.42885910E-05 9.98822800E-09-6.89907657E-13    2
% -1.16780396E+05-1.85960433E+02 1.43150898E+01 1.31672243E-01 1.32462853E-04    3
% -2.35743556E-07 9.10507438E-11-1.04820821E+05-1.33797027E+01                   4
%
% 	! H-abstraction on carbon 2 from saturated FAMEs (CH3-O-(C=O)-CH.-CH2-R species)(see 3.1.1., 3.1.1.1. sections)
% Species c13h25o2_laur_m_est_2j Fit:
% c13h25o2_laur_m_        H  25C  13O   2     G   300.000  5000.000 1000.00      1
%  2.75449350E+01 8.21814020E-02-3.32104533E-05 6.13045790E-09-4.24502731E-13    2
% -6.97530550E+04-1.03555711E+02 7.74925646E+00 9.79494516E-02 3.82594691E-05    3
% -1.04930822E-07 4.31944833E-11-6.26231861E+04 7.80066121E+00                   4
% Species c15h29o2_myrist_m_est_2j Fit:
% c15h29o2_myrist_        H  29C  15O   2     G   300.000  5000.000 1000.00      1
%  3.26253256E+01 9.30510763E-02-3.73695431E-05 6.86667022E-09-4.73851506E-13    2
% -7.69811947E+04-1.30204783E+02 1.09319777E+01 9.87183525E-02 7.57887158E-05    3
% -1.49678284E-07 5.89389160E-11-6.85872194E+04-5.27039047E+00                   4
% Species c17h33o2_palmit_m_est_2j Fit:
% c17h33o2_palmit_        H  33C  17O   2     G   300.000  5000.000 1000.00      1
%  3.47969412E+01 1.08174887E-01-4.36560662E-05 8.05036686E-09-5.57001848E-13    2
% -8.33532968E+04-1.36711346E+02 1.18429649E+01 1.08281798E-01 9.37470576E-05    3
% -1.75260710E-07 6.81980157E-11-7.41770519E+04-3.04443036E+00                   4
% Species c18h35o2_margar_m_est_2j Fit:
% c18h35o2_margar_        H  35C  18O   2     G   300.000  5000.000 1000.00      1
%  3.73611467E+01 1.13903983E-01-4.59591485E-05 8.47503001E-09-5.86429211E-13    2
% -8.69280667E+04-1.51416645E+02 1.28009370E+01 1.14599835E-01 9.93145536E-05    3
% -1.85919092E-07 7.23983479E-11-7.71387754E+04-8.54158312E+00                   4
% Species c19h37o2_stear_m_est_2j Fit:
% c19h37o2_stear_m        H  37C  19O   2     G   300.000  5000.000 1000.00      1
%  3.88001836E+01 1.20970409E-01-4.88425039E-05 9.00921828E-09-6.23447490E-13    2
% -9.01343860E+04-1.56022686E+02 1.25388297E+01 1.26045945E-01 9.34990125E-05    3
% -1.85855006E-07 7.30850786E-11-7.98836212E+04-4.33436425E+00                   4
% Species c21h41o2_arachid_m_est_2j Fit:
% c21h41o2_arachid        H  41C  21O   2     G   300.000  5000.000 1000.00      1
%  4.44177597E+01 1.31220602E-01-5.27329358E-05 9.69428587E-09-6.69211431E-13    2
% -9.75250445E+04-1.85644170E+02 1.28302810E+01 1.43260808E-01 1.00671318E-04    3
% -2.06884926E-07 8.20530187E-11-8.54920633E+04-4.67541735E+00                   4
%
% 	! H-abstraction from monounsaturated FAMEs (11j:on carbon 11;2j:on carbon 2;8j:on carbon 8)(see 3.1.2., 3.1.2.1., 3.1.2.2. sections)
% Species c17h31o2_palmitole_m_est_11j Fit:
% c17h31o2_palmito        H  31C  17O   2     G   300.000  5000.000 1000.00      1
%  3.42860856E+01 1.03351981E-01-4.17149510E-05 7.69246466E-09-5.32219654E-13    2
% -7.27264352E+04-1.32679107E+02 1.11206200E+01 1.04980126E-01 9.23934076E-05    3
% -1.72746825E-07 6.73360322E-11-6.35416562E+04 1.83930341E+00                   4
% Species c17h31o2_palmitole_m_est_2j Fit:
% c17h31o2_palmito        H  31C  17O   2     G   300.000  5000.000 1000.00      1
%  3.32908324E+01 1.12584694E-01-4.55839826E-05 8.42611891E-09-5.84037805E-13    2
% -6.83963296E+04-1.36807242E+02 8.18188826E+00 1.03010496E-01 1.33792275E-04    3
% -2.21168027E-07 8.43169922E-11-5.78740421E+04 1.18313933E+01                   4
% Species c17h31o2_palmitole_m_est_8j Fit:
% c17h31o2_palmito        H  31C  17O   2     G   300.000  5000.000 1000.00      1
%  3.38189966E+01 1.03803861E-01-4.17698167E-05 7.68682122E-09-5.31067680E-13    2
% -7.25156376E+04-1.30777258E+02 1.05799819E+01 1.09245899E-01 8.13381567E-05    3
% -1.61899182E-07 6.37439382E-11-6.34921337E+04 3.21606013E+00                   4
% Species c18h33o2_margarole_m_est_11j Fit:
% c18h33o2_margaro        H  33C  18O   2     G   300.000  5000.000 1000.00      1
%  3.60767842E+01 1.09867595E-01-4.43359807E-05 8.17496750E-09-5.65573398E-13    2
% -7.61029792E+04-1.39550378E+02 1.06830434E+01 1.15364283E-01 9.15364018E-05    3
% -1.78484897E-07 7.01189615E-11-6.62203173E+04 6.97931134E+00                   4
% Species c18h33o2_margarole_m_est_2j Fit:
% c18h33o2_margaro        H  33C  18O   2     G   300.000  5000.000 1000.00      1
%  3.58708970E+01 1.09624794E-01-4.41248211E-05 8.12211837E-09-5.61247536E-13    2
% -7.08891179E+04-1.35023651E+02 1.28405420E+01 1.08712030E-01 9.67955999E-05    3
% -1.78859013E-07 6.94425811E-11-6.16313377E+04-6.56954973E-01                   4
% Species c18h33o2_margarole_m_est_8j Fit:
% c18h33o2_margaro        H  33C  18O   2     G   300.000  5000.000 1000.00      1
%  3.55702553E+01 1.10375329E-01-4.44156545E-05 8.17390073E-09-5.64725703E-13    2
% -7.58704865E+04-1.35781913E+02 1.06266307E+01 1.17379688E-01 8.42330169E-05    3
% -1.70362020E-07 6.72617893E-11-6.62432546E+04 7.74919187E+00                   4
% Species c19h35o2_ole_m_est_11j Fit:
% c19h35o2_ole_m_e        H  35C  19O   2     G   300.000  5000.000 1000.00      1
%  3.80837111E+01 1.16163722E-01-4.67371030E-05 8.60036801E-09-5.94163979E-13    2
% -7.94133887E+04-1.49455861E+02 9.55411685E+00 1.34933086E-01 6.81323702E-05    3
% -1.63328294E-07 6.62252547E-11-6.89400192E+04 1.20201930E+01                   4
% Species c19h35o2_ole_m_est_2j Fit:
% c19h35o2_ole_m_e        H  35C  19O   2     G   300.000  5000.000 1000.00      1
%  3.85345397E+01 1.15581406E-01-4.66096518E-05 8.59185996E-09-5.94367839E-13    2
% -7.44747428E+04-1.51538218E+02 1.14785115E+01 1.29062681E-01 7.52826938E-05    3
% -1.67412542E-07 6.70924423E-11-6.43263952E+04 2.67718831E+00                   4
% Species c19h35o2_ole_m_est_8j Fit:
% c19h35o2_ole_m_e        H  35C  19O   2     G   300.000  5000.000 1000.00      1
%  3.83418393E+01 1.15872570E-01-4.66047845E-05 8.57355733E-09-5.92166553E-13    2
% -7.95171154E+04-1.51905055E+02 1.10715515E+01 1.26541641E-01 8.50097296E-05    3
% -1.77581532E-07 7.05496259E-11-6.91424538E+04 4.26134009E+00                   4
