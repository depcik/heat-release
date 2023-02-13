function RHS = combustfafl(t, Z)
global varode
NM = varode(4);
RMLT = [1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1];
%% Masses being solved
% index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
muO2 = Z(1);
muN2 = Z(2);
muAr = Z(3);
muCO2 = Z(4);
muH2O = Z(5);
muCO = Z(6);
muH2 = Z(7);
muOH = Z(8);
muH = Z(9);
muO = Z(10);
mbO2 = Z(NM+1);
mbN2 = Z(NM+2);
mbAr = Z(NM+3);
mbCO2 = Z(NM+4);
mbH2O = Z(NM+5);
mbCO = Z(NM+6);
mbH2 = Z(NM+7);
mbOH = Z(NM+8);
mbH = Z(NM+9);
mbO = Z(NM+10);
mfa = Z(2*NM+1);
mfl = Z(2*NM+2);
%% Variables
iS = varode(1);    % Combustion model: 0-single reaction, 1-expanded
iofa = varode(2);  % 0-added gas is not burning, 1-added gas is burning
iofl = varode(3);  % 0-direct injected fuel is not burning, 1-direct injected fuel is burning
% Added gas fuel properties
a = varode(5);
b = varode(6);
c = varode(7);
d = varode(8);
w = varode(9);
x = varode(10);
y = varode(11);
z = varode(12);
Kfa = varode(13);
Efa = varode(14);
Mfa = varode(15);
Kfl = varode(16);
Efl = varode(17);
Mfl = varode(18);
rhoCV = varode(19);
VCV = varode(20);           % [m3]
Runiv = varode(21);
TCV = varode(22);
for kk=1:NM
    MM(kk) = varode(22+kk);     % Molecular masses of the species
end
% Global combustion reaction parameters
delta = varode(23+NM);
epsilon = varode(24+NM);
g2 = varode(25+NM);
g3 = varode(26+NM);
g4 = varode(27+NM);
g5 = varode(28+NM);
g6 = varode(29+NM);
g7 = varode(30+NM);
g8 = varode(31+NM);
g9 = varode(32+NM);
fuelinj = varode(33+NM); % Liquid fuel being injected [kg/s]
nYO2fa = varode(34+NM); % Power term on oxygen for added gas fuel [-]
nYO2fl = varode(35+NM); % Power term on oxygen for liquid fuel [-]
% Set the derivatives all to zero first
for kk=1:(2*NM+2)
    RHS(kk) = 0;
end
RHS = RHS'; % To return in the proper Matlab column/row format
% Fill the molecular masses
% index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
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
% Compute the derivatives
% Determine the mass fractions of added gas and oxygen in the CV
mtot = 0;
for j=1:NM
    mtot = mtot + Z(j);
    mtot = mtot + Z(j+NM);
end
mtot = mtot + Z(2*NM+1);
mtot = mtot + Z(2*NM+2);
YfaCV = mfa/mtot;
YflCV = mfl/mtot;
YO2CV = (mbO2+muO2)/mtot;
%% Check which combustion reactions to solve
if (iS == 0) % Single reaction
    %% Single combustion reactions to completion
    % CaHbOcNd + n2*O2 => n3*CO2 + n4*H2O + n5*N2 (added gas)
    % CwHxOyNz + l2*O2 => l3*CO2 + l4*H2O + l5*N2 (direct injected fuel)
    % unburned => burned
    % CO2, H2O, and N2 are our burned products (completed combustion)
    % Added gas fuel
    % Compute n2, n3, n4, and n5 from the fuel formula
    n3 = a;
    n4 = b/2;
    n2 = (2*n3 + n4 - c)/2;
    n5 = d/2;
    % Liquid fuel
    l3 = w;
    l4 = x/2;
    l2 = (2*l3 + l4 - y)/2;
    l5 = z/2;
    % Added gas burning
    % We will assume that the density and temperature of the control volume
    % does not change significantly over the time-step. In reality, it is set
    % up in an implicit manner so it will recalculate these properties at the
    % end and then go back through using them until convergence
    if (iofa == 1)
        dmfabdt = -Kfa*(rhoCV^2)*YfaCV*(YO2CV^nYO2fa)*VCV*exp(-Efa/(Runiv*TCV)); % [kg/s]
    else
        dmfabdt = 0;
    end
    % Liquid fuel burning & being injected
    if (iofl == 1)
        %dmflbdt = -Kfl*(rhoCV^2)*YflCV*(YO2CV^5)*VCV*exp(-Efl/(Runiv*TCV)); % [kg/s]
        dmflbdt = -Kfl*(rhoCV^2)*YflCV*(YO2CV^nYO2fl)*VCV*exp(-Efl/(Runiv*TCV)); % [kg/s]
    else
        dmflbdt = 0;  % [kg/s]
    end
    % Add any fuel that has been injected
    dmfldt = dmflbdt + fuelinj;
    % Calculate all derivatives
    % We can now update our derivatives
    % Added gas fuel
    RHS(2*NM+1) = dmfabdt;  % [kg/s]
    % Liquid fuel - includes what has burned and what has been injected
    RHS(2*NM+2) = dmfldt;  % [kg/s]
    % To do so requires putting everything in molar format
    dnfabdt = dmfabdt/Mfa;  % [kg/s] * [mol/kg] = [mol,fa/s]
    dnflbdt = dmflbdt/Mfl;   % [mol,fl/s]
    % The localized combustion reaction for the added gas fuel
    % Unburned oxygen is being converted
    dnuO2 = n2*dnfabdt + l2*dnflbdt;
    RHS(1) = dnuO2*MO2;     % [mol/s]*[kg/mol] = [kg unburned O2/s]
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
    % index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
    % As the oxygen converts to CO2 & H2O, it transforms nitrogen from the
    % unburned side to the burned side. Since the global combustion
    % reaction is different than the localized reaction, converting
    % nitrogen from unburned to burned via the global does not hold.
    % Instead, for every 1 gram of oxygen converted, there is muN2/muO2
    % grams of nitrogen converted. The nitrogen converts in the same mass
    % ratio to start.
    % These are only the case before combustion happens
    % mO2chk = (delta*g2 + epsilon*g5)*MO2;
    % mN2chk = (delta*g3 + epsilon*g8)*MN2;
    % nO2chk = (delta*g2 + epsilon*g5);
    % nN2chk = (delta*g3 + epslion*g8);
    % This is what we had originally
    % dnuN2 = (delta*g3 + epsilon*g8)/(delta*g2 + epsilon*g5)*dnuO2;
    % Convert to mass format
    % dmuN2 = (delta*g3 + epsilon*g8)/(delta*g2 + epsilon*g5)*dnuO2*MN2;
    % dmuN2 = (delta*g3 + epsilon*g8)/(delta*g2 + epsilon*g5)*dmuO2*MN2/MO2;
    % If we just write the mass ratios, we get the same result
    % mN2chk/mO2chk = [(delta*g3 + epsilon*g8)*MN2]/[(delta*g2 + epsilon*g5)*MO2];
    % So, my original calculation is correct right at the beginning.
    % However, as the oxygen converts, delta*g2 and eps*g5 decrease since
    % there are less moles of unburned oxygen
    % Holding to this original ratio will not work as it will overpredict
    % the conversion of N2, Ar, CO2, and H2O. Instead, the updated mass ratios need to be
    % calculated to ensure that they convert at the same rate
    RHS(2) = muN2/muO2*RHS(1);
    dnuN2 = RHS(2)/MN2;
    RHS(3) = muAr/muO2*RHS(1);
    dnuAr = RHS(3)/MAr;
    RHS(4) = muCO2/muO2*RHS(1);
    dnuCO2 = RHS(4)/MCO2;
    RHS(5) = muH2O/muO2*RHS(1);
    dnuH2O = RHS(5)/MH2O;
    %RHS(2) = dnuN2*MN2;     % [kg unburned N2/s]
    %dnuAr = (delta*g4 + epsilon*g9)/(delta*g2 + epsilon*g5)*dnuO2;
    %RHS(3) = dnuAr*MAr;     % [kg unburned Ar/s]
    %dnuCO2 = epsilon*g6/(delta*g2 + epsilon*g5)*dnuO2;
    %RHS(4) = dnuCO2*MCO2;   % [kg unburned CO2/s]
    %dnuH2O = epsilon*g7/(delta*g2 + epsilon*g5)*dnuO2;
    %RHS(5) = dnuH2O*MH2O;   % [kg unburned H2O/s]
    % Now, we are forming burned products of CO2, H2O, and N2 according to the
    % localized combustion reaction + some of the unburned side is moving onto
    % the burned side
    dnbCO2 = -n3*dnfabdt - l3*dnflbdt - dnuCO2;
    RHS(NM+4) = dnbCO2*MCO2;    % [kg burned CO2/s]
    dnbH2O = -n4*dnfabdt - l4*dnflbdt - dnuH2O;
    RHS(NM+5) = dnbH2O*MH2O;    % [kg burned H2O/s]
    dnbN2 = -n5*dnfabdt - l5*dnflbdt - dnuN2;
    RHS(NM+2) = dnbN2*MN2;      % [kg burned N2/s]
    dnbAr = -dnuAr;
    RHS(NM+3) = dnbAr*MAr;      % [kg burned Ar/s]
else % expanded combustion through CO and H2
    % CaHbOcNd + n2*O2 => n3*CO + n4*H2 + n5*N2 (added gas)
    % CwHxOyNz + l2*O2 => l3*CO + l4*H2 + l5*N2 (direct injected fuel)
    % Added gas
    n3 = a;
    n4 = b/2;
    n2 = 0.5*(n3-c);
    n5 = d/2;
    % Liquid fuel
    l3 = w;
    l4 = x/2;
    l2 = 0.5*(l3-y);
    l5 = z/2;
    % --- Concentrations are on a [mol/cm3] basis
    % rhoCV => [kg/m3] - the density of the entire control volume. We need
    % to convert this to [mol/m3] using the mixture molecular mass. So,
    % let's compute the overall control volume mass fractions that we will
    % convert into mole fractions
    for j=1:NM
        YCV(j) = Z(j)/mtot;             % unburned mass fractions
        YCV(j+NM) = Z(j+NM)/mtot;       % burned mass fractions
    end
    YCV(2*NM+1) = Z(2*NM+1)/mtot;       % added gas fuel
    YCV(2*NM+2) = Z(2*NM+2)/mtot;       % liquid fuel
    % mole fractions now
    btm = 0;
    for j=1:NM
        btm = btm + YCV(j)/MM(j) + YCV(j+NM)/MM(j);
    end
    btm = btm + YCV(2*NM+1)/Mfa + YCV(2*NM+2)/Mfl;
    for j=1:NM
        XCV(j) = YCV(j)/MM(j)/btm;
        XCV(j+NM) = YCV(j+NM)/MM(j)/btm;
    end
    XCV(2*NM+1) = YCV(2*NM+1)/Mfa/btm;
    XCV(2*NM+2) = YCV(2*NM+2)/Mfl/btm;
    % Mixture molecular mass [kg/mol]
    Mmix = 0;
    for j=1:NM
        Mmix = Mmix + XCV(j)*MM(j) + XCV(j+NM)*MM(j);
    end
    Mmix = Mmix + XCV(2*NM+1)*Mfa + XCV(2*NM+2)*Mfl;
    % Thus the molar density is 
    % [mol/cm3] = [kg/m3]*[mol/kg]*[m/100 cm]^3
    rhomCVcm3 = rhoCV/Mmix/(100^3);
    % The cylinder volume on a [cm3] basis
    % [cm3] = [m3]*[100 cm/m]^3
    VCVcm3 = VCV*(100^3);
    % Added gas burning
    if (iofa == 1)
        dmfabdt = -Kfa*(rhoCV^2)*YfaCV*(YO2CV^5)*VCV*exp(-Efa/(Runiv*TCV)); % [kg/s]
    else
        dmfabdt = 0;
    end
    % Liquid fuel burning & being injected
    if (iofl == 1)
        dmflbdt = -Kfl*(rhoCV^2)*YflCV*(YO2CV^5)*VCV*exp(-Efl/(Runiv*TCV)); % [kg/s]
    else
        dmflbdt = 0;  % [kg/s]
    end
    % Add any fuel that has been injected
    dmfldt = dmflbdt + fuelinj;
    % Calculate all derivatives
    % We can now update our derivatives
    % Added gas fuel
    RHS(2*NM+1) = dmfabdt;          % [kg/s]
    % Liquid fuel - includes what has burned and what has been injected
    RHS(2*NM+2) = dmfldt;           % [kg/s]
    % To do so requires putting everything in molar format
    dnfabdt = dmfabdt/Mfa;          % [kg/s] * [mol/kg] = [mol,fa/s]
    dnflbdt = dmflbdt/Mfl;          % [mol,fl/s]
    % Unburned oxygen is being converted through the fuels
    dnuO2_fuel = n2*dnfabdt + l2*dnflbdt;        % We are converting "unburned" O2
    dnuCO_fuel = -n3*dnfabdt - l3*dnflbdt;       % We are producing "unburned" CO
    dnuH2_fuel = -n4*dnfabdt - l4*dnflbdt;       % We are producing "unburned" H2
    dnbN2_fuel = -n5*dnfabdt - l5*dnflbdt;       % We are producing "burned" N2
    % The rest of the reactions follow more conventional homogeneous
    % reaction kinetics. A good example as to how to write these reactions
    % can be found here: https://doi.org/10.1115/IMECE2019-10028
    % Example with methane
    % Overall localized combustion reaction:
    % CH4 + 2O2 => CO2 + 2H2O 
    % Now, we break it up:
    % CH4 + 0.5O2 => CO + 2H2
    % Thus, to complete the reaction we need
    % CO + 0.5O2 => CO2
    % 2*(H2 + 0.5O2 => H2O)
    % We then write the homogeneous set of reactions to achieve these
    % combustion reactions
    % --- H2 + 0.5*O2 => H2O
    % R1: 1.5 *     (H + O2 => O + OH)
    % R2: 0.25 *    (O + H2 => H + OH)
    % R3: 0.25 *    (OH + H2 => H + H2O)
    % R4: 1 *       (O + H2O => OH + OH)
    % R5: 0.5 *     (H2 => 2H)
    % R6: 1 *       (2O => O2)
    % R7: -1.75 *   (O + H => OH)
    % R8: 1.75 *    (H + OH => H2O)
    % --- CO + 0.5*O2 => CO2
    % R9: 0.5 *     (CO + OH => CO2 + H)
    % R10: 0.5 *    (CO + O2 => O + CO2)
    % R11: 0.5 *    (O + H => OH)
    % The only species that show up in the burned products are CO2 and H2O
    % The others remain "unburned". In addition, for this first example, we
    % are only concerning ourselves with the forward reactions indicated.
    % Now, we have: (n3+l3)*CO + (n4+l4)*H2
    f = n3+l3;
    g = n4+l4;
    % Let's call this: fCO + gH2 and these need to go to CO2 and H2O. As a
    % result, we multiply R1-R8 by g and R9-R11 by f
    % Originally, we have fuel, O2, N2, Ar, CO2, and H2O from air, residual, and EGR
    % fuel and O2 => CO and H2 : we put CO and H2 in the unburned mixture as it is not "burnt" yet
    % CO & H2 turn into H, O, OH, CO2, and H2O : we put H, O, and OH into
    % the unburned pile. However, CO2 and H2O we put in the burned pile.
    % But, H2O does contribute to the chain reactions. So, to compensate
    % for our concentration of CO2 and H2O that impacts the combustion
    % reaction, we will include both the unburned and burned piles
    % Now, the concentrations of the unburned species on a [mol/cm3] basis
    % Only the unburned species contribute to the combustion reactions
    CuO2 = rhomCVcm3*XCV(1);
    CuN2 = rhomCVcm3*XCV(2);
    CuAr = rhomCVcm3*XCV(3);
    CubCO2 = rhomCVcm3*(XCV(4)+XCV(4+NM));
    CubH2O = rhomCVcm3*(XCV(5)+XCV(5+NM));
    CuCO = rhomCVcm3*XCV(6);
    CuH2 = rhomCVcm3*XCV(7);
    CuOH = rhomCVcm3*XCV(8);
    CuH = rhomCVcm3*XCV(9);
    CuO = rhomCVcm3*XCV(10);
    % Universal gas constant units [J/(mol*K)]
    % Reaction 1
    % R1: 1.5 *     (H + O2 => O + OH)
    dg1 = [-8.594962E-18	-2.480531E-13	2.966302E-10	8.659237E-07	1.808994E-04	-2.593436E+01	7.049391E+04];
    DG1 = dg1(1)*TCV^6  + dg1(2)*TCV^5 + dg1(3)*TCV^4 + dg1(4)*TCV^3 + dg1(5)*TCV^2 + dg1(6)*TCV + dg1(7); % [J/mol]
    KpK1 = exp(-DG1/(Runiv*TCV));
    KcK1 = KpK1*(101325/(Runiv*TCV))^(1+1-1-1);       % Concentration units
    Ea1 = 16.44*4186.8;     % [J/mol]
    A1 = 1.91*10^14;        % [cm3/mol/s]
    n1 = 0;
    k1f = A1*(TCV^n1)*exp(-Ea1/(Runiv*TCV));
    k1r = k1f/KcK1;
    R1f = RMLT(1)*k1f*CuH*CuO2;
    R1r = RMLT(2)*k1r*CuO*CuOH;
    R1 = g*1.5*(R1f-R1r);
    dnuH_R1 = -R1*VCVcm3;       % [mol/s]
    dnuO2_R1 = -R1*VCVcm3;      % [mol/s]
    dnuO_R1 = R1*VCVcm3;        % [mol/s]
    dnuOH_R1 = R1*VCVcm3;       % [mol/s]
    % Reaction 2
    % R2: 0.25 *    (O + H2 => H + OH)
    dg2 = [-2.107744E-18	1.923065E-13	-9.929900E-11	-1.303400E-06	2.356975E-03	-7.982448E+00	7.991874E+03];
    DG2 = dg2(1)*TCV^6  + dg2(2)*TCV^5 + dg2(3)*TCV^4 + dg2(4)*TCV^3 + dg2(5)*TCV^2 + dg2(6)*TCV + dg2(7); % [J/mol]
    KpK2 = exp(-DG2/(Runiv*TCV));
    KcK2 = KpK2*(101325/(Runiv*TCV))^(1+1-1-1);       % Concentration units
    Ea2 = 6.29*4186.8;      % [J/mol]
    A2 = 5.08*10^4;         % [cm3/mol/s]
    n2 = 2.67;
    k2f = A2*(TCV^n2)*exp(-Ea2/(Runiv*TCV));
    k2r = k2f/KcK2;
    R2f = RMLT(3)*k2f*CuO*CuH2;   % [mol/s]
    R2r = RMLT(4)*k2r*CuH*CuOH;
    R2 = g*0.25*(R2f-R2r);
    dnuO_R2 = -R2*VCVcm3;
    dnuH2_R2 = -R2*VCVcm3;
    dnuH_R2 = R2*VCVcm3;
    dnuOH_R2 = R2*VCVcm3;
    % Reaction 3
    % R3: 0.25 *    (OH + H2 => H + H2O)
    dg3 = [1.680276E-17	-8.648821E-14	1.548565E-09	-6.239563E-06	9.519679E-03	7.173586E+00	-6.245794E+04];
    DG3 = dg3(1)*TCV^6  + dg3(2)*TCV^5 + dg3(3)*TCV^4 + dg3(4)*TCV^3 + dg3(5)*TCV^2 + dg3(6)*TCV + dg3(7); % [J/mol]
    KpK3 = exp(-DG3/(Runiv*TCV));
    KcK3 = KpK3*(101325/(Runiv*TCV))^(1+1-1-1);       % Concentration units
    Ea3 = 3.43*4186.8;
    A3 = 2.16*10^8;
    n3 = 1.51;
    k3f = A3*(TCV^n3)*exp(-Ea3/(Runiv*TCV));
    k3r = k3f/KcK3;
    R3f = RMLT(5)*k3f*CuOH*CuH2;
    R3r = RMLT(6)*k3r*CuH*CubH2O;
    R3 = g*0.25*(R3f-R3r);
    dnuOH_R3 = -R3*VCVcm3;
    dnuH2_R3 = -R3*VCVcm3;
    dnuH_R3 = R3*VCVcm3;
    dnbH2O_R3 = R3*VCVcm3;          % H2O shows up in the burned pile
    % Reaction 4
    % R4: 1 *       (O + H2O => OH + OH)
    dg4 = [-1.891050E-17	2.787948E-13	-1.647864E-09	4.936163E-06	-7.162704E-03	-1.515603E+01	7.044981E+04];
    DG4 = dg4(1)*TCV^6  + dg4(2)*TCV^5 + dg4(3)*TCV^4 + dg4(4)*TCV^3 + dg4(5)*TCV^2 + dg4(6)*TCV + dg4(7); % [J/mol]
    KpK4 = exp(-DG4/(Runiv*TCV));
    KcK4 = KpK4*(101325/(Runiv*TCV))^(1+1-1-1);       % Concentration units
    Ea4 = 13.4*4186.8;
    A4 = 2.97*10^6;
    n4 = 2.02;
    k4f = A4*(TCV^n4)*exp(-Ea4/(Runiv*TCV));
    k4r = k4f/KcK4;
    R4f = RMLT(7)*k4f*CuO*CubH2O; 
    R4r = RMLT(8)*k4r*CuOH^2;
    R4 = g*1*(R4f-R4r);
    dnuO_R4 = -R4*VCVcm3;
    dnbH2O_R4 = -R4*VCVcm3;
    dnuOH_R4 = 2*R4*VCVcm3;
    % Reaction 5
    % R5: 0.5 *     (H2 + M => 2H + M)
    dg5 = [-2.590861E-17	4.106902E-13	-2.725284E-09	1.003505E-05	-2.326487E-02	-8.800201E+01	4.346409E+05];
    DG5 = dg5(1)*TCV^6  + dg5(2)*TCV^5 + dg5(3)*TCV^4 + dg5(4)*TCV^3 + dg5(5)*TCV^2 + dg5(6)*TCV + dg5(7); % [J/mol]
    KpK5 = exp(-DG5/(Runiv*TCV));
    KcK5 = KpK5*(101325/(Runiv*TCV))^(2-1);       % Concentration units
    Ea5 = 104.38*4186.8;
    A5 = 4.58*10^19;
    n5 = -1.4;
    k5f = A5*(TCV^n5)*exp(-Ea5/(Runiv*TCV));
    k5r = k5f/KcK5;
    % Efficiency factors need to be included
    CM = CuO2 + CuN2 + 0.75*CuAr + CubCO2 + 12*CubH2O + CuCO + 2.5*CuH2 + CuOH + CuH;
    R5f = RMLT(9)*k5f*CuH2*CM;
    R5r = RMLT(10)*k5r*CuH^2*CM;
    R5 = g*0.5*(R5f-R5r);
    dnuH2_R5 = -R5*VCVcm3;
    dnuH_R5 = 2*R5*VCVcm3;
    % Reaction 6
    % R6: 1 *       (2O + M => O2 + M)
    dg6 = [4.162680E-17	-6.223858E-13	3.799118E-09	-1.234049E-05	2.363014E-02	1.073359E+02	-4.974186E+05];
    DG6 = dg6(1)*TCV^6  + dg6(2)*TCV^5 + dg6(3)*TCV^4 + dg6(4)*TCV^3 + dg6(5)*TCV^2 + dg6(6)*TCV + dg6(7); % [J/mol]
    KpK6 = exp(-DG6/(Runiv*TCV));
    KcK6 = KpK6*(101325/(Runiv*TCV))^(1-2);       % Concentration units
    Ea6 = 0;
    A6 = 6.16*10^15;
    n6 = -0.5;
    k6f = A6*(TCV^n6)*exp(-Ea6/(Runiv*TCV));
    k6r = k6f/KcK6;
    R6f = RMLT(11)*k6f*CuO^2*CM;
    R6r = RMLT(12)*k6r*CuO2*CM;
    R6 = g*1*(R6f-R6r);
    dnuO_R6 = -2*R6*VCVcm3;
    dnuO2_R6 = R6*VCVcm3;
    % Reaction 7
    % The reaction rate is written as: O + H + M => OH + M
    % But, we need 1.75*(OH + M => O + H + M)
    % So, we compute the forward rate constant
    % Then, we use the equilibrium constant in concentration units to
    % determine the reverse rate constant
    % The curve-fit developed is for DG [J/mol] - Gibbs free energy of the
    % expression as written
    dg7 = [2.895933E-17	-4.575052E-13	3.012610E-09	-1.087605E-05	2.433332E-02	8.084047E+01	-4.268032E+05];
    DG7 = dg7(1)*TCV^6  + dg7(2)*TCV^5 + dg7(3)*TCV^4 + dg7(4)*TCV^3 + dg7(5)*TCV^2 + dg7(6)*TCV + dg7(7); % [J/mol]
    KpK7 = exp(-DG7/(Runiv*TCV));
    KcK7 = KpK7*(101325/(Runiv*TCV))^(1-1-1);       % Concentration units
    Ea7 = 0;
    A7 = 4.71*10^18;
    n7 = -1;
    k7f = A7*(TCV^n7)*exp(-Ea7/(Runiv*TCV));
    k7r = k7f/KcK7;             % Reverse rate constant
    % OK, but we want 1.75*(OH + M => O + H + M). So, we will write the
    % reactions in this manner
    R7f = RMLT(14)*k7r*CuOH*CM;
    R7r = RMLT(13)*k7f*CuO*CuH*CM;
    R7 = g*1.75*(R7f - R7r);
    dnuOH_R7 = -R7*VCVcm3;
    dnuO_R7 = R7*VCVcm3;
    dnuH_R7 = R7*VCVcm3;
    % Reaction 8
    % R8: 1.75 *    (H + OH + M => H2O + M)
    dg8 = [4.786983E-17	-7.363000E-13	4.660474E-09	-1.581221E-05	3.149603E-02	9.599650E+01	-4.972530E+05];
    DG8 = dg8(1)*TCV^6  + dg8(2)*TCV^5 + dg8(3)*TCV^4 + dg8(4)*TCV^3 + dg8(5)*TCV^2 + dg8(6)*TCV + dg8(7); % [J/mol]
    KpK8 = exp(-DG8/(Runiv*TCV));
    KcK8 = KpK8*(101325/(Runiv*TCV))^(1-1-1);       % Concentration units
    Ea8 = 0;
    A8 = 2.21*10^22;
    n8 = -2;
    k8f = A8*(TCV^n8)*exp(-Ea8/(Runiv*TCV));
    k8r = k8f/KcK8;
    R8f = RMLT(15)*k8f*CuH*CuOH*CM;
    R8r = RMLT(16)*k8r*CubH2O*CM;
    R8 = g*1.75*(R8f-R8r);
    dnuH_R8 = -R8*VCVcm3;
    dnuOH_R8 = -R8*VCVcm3;
    dnbH2O_R8 = R8*VCVcm3;
    % Reaction 9
    % GRI-mech
    % with concentration units mol/cm3. The units of A are cm3/mol/s, T is in K, and E is in cal/mol. For termolecular recombination reactions, the units of A are cm6/mol2/s.
    % http://combustion.berkeley.edu/gri-mech/data/k_form.html
    % http://combustion.berkeley.edu/gri-mech/version30/files30/grimech30.dat
    % OH+CO<=>H+CO2                            4.760E+07    1.228      70.00
    % R9: 0.5 *     (CO + OH => CO2 + H)
    dg9 = [3.115726E-17	4.115625E-13	-1.823856E-09	3.917183E-06	-1.079668E-02	6.001400E+01	-1.053410E+05];
    DG9 = dg9(1)*TCV^6  + dg9(2)*TCV^5 + dg9(3)*TCV^4 + dg9(4)*TCV^3 + dg9(5)*TCV^2 + dg9(6)*TCV + dg9(7); % [J/mol]
    KpK9 = exp(-DG9/(Runiv*TCV));
    KcK9 = KpK9*(101325/(Runiv*TCV))^(1+1-1-1);       % Concentration units
    Ea9 = 70*4186.8/1000;
    A9 = 4.760*10^7;
    n9 = 1.228;
    k9f = A9*(TCV^n9)*exp(-Ea9/(Runiv*TCV));
    k9r = k9f/KcK9;
    R9f = RMLT(17)*k9f*CuCO*CuOH;
    R9r = RMLT(18)*k9r*CubCO2*CuH;
    R9 = f*0.5*(R9f-R9r);
    dnuCO_R9 = -R9*VCVcm3;
    dnuOH_R9 = -R9*VCVcm3;
    dnbCO2_R9 = R9*VCVcm3;
    dnuH_R9 = R9*VCVcm3;
    % Reaction 10
    % R10: 0.5 *    (CO + O2 => O + CO2)
    % O2+CO<=>O+CO2                            2.500E+12     .000   47800.00
    dg10 = [1.848979E-17	5.764431E-13	-2.610364E-09	5.381624E-06	-1.009349E-02	3.351859E+01	-3.472562E+04];
    DG10 = dg10(1)*TCV^6  + dg10(2)*TCV^5 + dg10(3)*TCV^4 + dg10(4)*TCV^3 + dg10(5)*TCV^2 + dg10(6)*TCV + dg10(7); % [J/mol]
    KpK10 = exp(-DG10/(Runiv*TCV));
    KcK10 = KpK10*(101325/(Runiv*TCV))^(1+1-1-1);       % Concentration units
    Ea10 = 47800*4186.8/1000;
    A10 = 2.5*10^12;
    n10 = 0;
    k10f = A10*(TCV^n10)*exp(-Ea10/(Runiv*TCV));
    k10r = k10f/KcK10;
    R10f = RMLT(19)*k10f*CuO2*CuCO;
    R10r = RMLT(20)*k10r*CuO*CubCO2;
    R10 = f*0.5*(R10f-R10r);
    dnuO2_R10 = -R10*VCVcm3;
    dnuCO_R10 = -R10*VCVcm3;
    dnuO_R10 = R10*VCVcm3;
    dnbCO2_R10 = R10*VCVcm3;
    % Reaction 11 - this is a repeat of Reaction 7
    % R11: 0.5 *    (O + H + M => OH + M)
    % Since there is already the forward and reverse reaction rates
    % computed
    R11f = RMLT(13)*k7f*CuO*CuH*CM;
    R11r = RMLT(14)*k7r*CuOH*CM;
    R11 = f*0.5*(R11f-R11r);
    dnuO_R11 = -R11*VCVcm3;
    dnuH_R11 = -R11*VCVcm3;
    dnuOH_R11 = R11*VCVcm3;
    % --- Sum up all the unburned species conversions
    % index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
    dnuO2 = dnuO2_fuel + dnuO2_R1 + dnuO2_R6 + dnuO2_R10;      % [mol/s]
    dnuCO = dnuCO_fuel + dnuCO_R9 + dnuCO_R10;
    dnuH2 = dnuH2_fuel + dnuH2_R2 + dnuH2_R3 + dnuH2_R5;
    dnuOH = dnuOH_R1 + dnuOH_R2 + dnuOH_R3 + dnuOH_R4 + dnuOH_R7 + dnuOH_R8 + dnuOH_R9 + dnuOH_R11;
    dnuH = dnuH_R1 + dnuH_R2 + dnuH_R3 + dnuH_R5 + dnuH_R7 + dnuH_R8 + dnuH_R9 + dnuH_R11;
    dnuO = dnuO_R1 + dnuO_R2 + dnuO_R4 + dnuO_R6 + dnuO_R7 + dnuO_R10 + dnuO_R11;
    % --- As before, nitrogen, argon, carbon dioxide, and water that is
    % originally mixed with the oxygen moves to the burned side
    dnuN2 = muN2/muO2*MO2/MN2*dnuO2;
    dnuAr = muAr/muO2*MO2/MAr*dnuO2;
    dnuCO2 = muCO2/muO2*MO2/MCO2*dnuO2;
    dnuH2O = muH2O/muO2*MO2/MH2O*dnuO2;
    %dnuN2 = (delta*g3 + epsilon*g8)/(delta*g2 + epsilon*g5)*dnuO2;
    %dnuAr = (delta*g4 + epsilon*g9)/(delta*g2 + epsilon*g5)*dnuO2;
    %dnuCO2 = epsilon*g6/(delta*g2 + epsilon*g5)*dnuO2;
    %dnuH2O = epsilon*g7/(delta*g2 + epsilon*g5)*dnuO2;
    % --- Sum up all the burned species conversions
    dnbCO2 = -dnuCO2 + dnbCO2_R10 + dnbCO2_R9;
    dnbH2O = -dnuH2O + dnbH2O_R8 + dnbH2O_R4 + dnbH2O_R3;
    dnbN2 = dnbN2_fuel - dnuN2;
    dnbAr = -dnuAr;
    % Fill the right hand side terms
    % index: 1-O2, 2-N2, 3-Ar, 4-CO2, 5-H2O, 6-CO, 7-H2, 8-OH, 9-H, 10-O
    % unburned
    RHS(1) = dnuO2*MO2;
    RHS(2) = dnuN2*MN2;
    RHS(3) = dnuAr*MAr;
    RHS(4) = dnuCO2*MCO2;
    RHS(5) = dnuH2O*MH2O;
    RHS(6) = dnuCO*MCO;
    RHS(7) = dnuH2*MH2;
    RHS(8) = dnuOH*MOH;
    RHS(9) = dnuH*MH;
    RHS(10) = dnuO*MO;
    % burned
    RHS(NM+2) = dnbN2*MN2;
    RHS(NM+3) = dnbAr*MAr;
    RHS(NM+4) = dnbCO2*MCO2;
    RHS(NM+5) = dnbH2O*MH2O;
end