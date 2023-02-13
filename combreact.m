function [cr] = combreact(iS, inp, calc, fuel)
% inp(21) = Mass flow rate of air [kg/s]
% inp(113) = Molecular mass of air [kg/mol]
% inp(23) = Mass flow rate of fuel added during first injection event [kg/s]
% inp(27), inp(31), inp(35), inp(39) = Mass flow rates of different fuel
% injection events [kg/s]
% fuel(5) = Molecular mass of fuel [kg/mol]
% inp(46) = Added gas mass flow rate [kg/s]
% fuel(24) = Molecular mass of added gas [kg/mol]
% inp(64) = Combustion efficiency of added gas [%]
% inp(65) = Combustion efficiency of liquid fuel [%]
% calc(13) = Cycles per second
%--- Mass and moles of air, liquid fuel, gaseous fuel added per cycle
cr(1) = inp(21)/calc(13);       % Air mass per cycle [kg/cycle]
cr(2) = cr(1)/inp(113);         % Air moles per cycle [mol/cycle]
cr(3) = (inp(23)+inp(27)+inp(31)+inp(35)+inp(39))/calc(13);  % Total liquid fuel mass per cycle [kg/cycle]
cr(4) = cr(3)/fuel(5);          % Liquid fuel moles per cycle [mol/cycle]
if (iS(2) == 1) % Added gas
    cr(5) = inp(46)/calc(13);   % Added gas fuel mass per cycle [kg/cycle]
    cr(6) = cr(5)/fuel(24);     % Added gas molar fuel per cycle [kg/cycle]
    mmgas = fuel(24);           % Molecular mass of fuel
else
    cr(5) = 0;
    cr(6) = 0;
    mmgas = 0;
end
%--- Let's make our lives easier by using variables that conform to our global model reactions
PIVC = calc(15);            % Pressure at IVC [Pa]
VIVC = calc(16);            % Volume at IVC [m3]
Xip = cr(4);                % [mol/cycle]
Xig = cr(6);                % [mol/cycle]
delta = cr(2);              % [mol/cycle]
Zetap = (1-inp(65))*Xip;    % [mol/cycle]
Zetag = (1-inp(64))*Xig;    % [mol/cycle]
w = fuel(1);
x = fuel(2);
y = fuel(3);
z = fuel(4);
if (iS(2) == 1)
    a = fuel(20);
    b = fuel(21);
    c = fuel(22);
    d = fuel(23);
else
    a = 0;
    b = 0;
    c = 0;
    d = 0;
end
g2 = inp(110);
g3 = inp(111);
g4 = inp(112);
mmfuel = fuel(5);
Runiv = inp(100);
MAr = inp(105);
MCO2 = inp(101) + 2*inp(102);
MN2 = 2*inp(103);
MH2O = 2*inp(104) + inp(102);
MO2 = 2*inp(102);
egr = inp(63);
res = inp(70);
% Note: EGR is on a volume basis and residual fraction is on a mass basis
g6 = (Xig-Zetag)*a + (Xip-Zetap)*w;
g7 = 0.5*((Xig-Zetag)*b + (Xip-Zetap)*x);
g5 = 0.5*(2*g2*delta + (Xig-Zetag)*c +((Xip-Zetap)*y - (g7 + 2*g6)));     % Slightly off from hand calcs
g8 = 0.5*(2*g3*delta + (Xig-Zetag)*d + (Xip-Zetap)*z);
g9 = delta*g4;
% Estimate our mass of residual
mres = res*(cr(1)+cr(5));   
% Get our first estimate of alpha
alp = mres/(Zetag*mmgas + MAr*g9 + MN2*g8 + MH2O*g7 + MCO2*g6 + MO2*g5);
% Then from EGR, we can solve for beta
beta = ((-g4-g3-g2)*egr*delta - (Xig + alp*Zetag + (g9+g8+g7+g6+g5)*alp)*egr)/...
    ((Zetag+g9+g8+g7+g6+g5)*egr - (Zetag+g9+g8+g7+g6+g5));
% We will converge on the combined values
eps = alp + beta;
% --- Now, we iterate until convergence
epsc = 1000;
while (epsc > 1e-10)
    epsold = eps;
    g6n = (-Zetag*a*eps + (Zetag-Xig)*a + (Zetap-Xip)*w)/(eps-1);
    g7n = (-Zetag*b*eps + (Zetag-Xig)*b + (Zetap-Xip)*x)/(2*eps - 2);
    g5n = ((-Zetag*c - (g7n + 2*g6n))*eps + ((Zetag-Xig)*c + (Zetap-Xip)*y + g7n + 2*g6n - 2*g2*delta))/(2*eps - 2);
    g8n = (-(Zetag*d*eps) + ((Zetag-Xig)*d + (Zetap-Xig)*z - 2*g3*delta))/(2*eps - 2);
    g9n = -g4*delta/(eps-1);
    % Compute the mass of EGR at IVC
    megr = Zetag*beta*mmgas + beta*(g5n*MO2 + g6n*MCO2 + g7n*MH2O + g8n*MN2 + g9n*MAr);
    % We have enough information to solve for the mass of residual
    mres = res*(-megr - cr(5) - cr(1))/(res-1);
    % The total intake mass is
    mint = cr(1) + cr(5) + mres + megr;
    % We can now find alpha and beta
    alp = mres/(Zetag*mmgas + MAr*g9n + MN2*g8n + MH2O*g7n + MCO2*g6n + MO2*g5n);
    % Then from EGR, we can solve for beta
    beta = ((-g4-g3-g2)*egr*delta - (Xig + alp*Zetag + (g9n+g8n+g7n+g6n+g5n)*alp)*egr)/...
        ((Zetag+g9n+g8n+g7n+g6n+g5n)*egr - (Zetag+g9n+g8n+g7n+g6n+g5n));
    eps = alp + beta;
    epsc = abs(epsold-eps);
end
% Our global combustion reaction has now been defined. Enough information
% has been provided to determine the temperature at IVC. 
Xfueln = 0;
Xgasn = (Xig + Zetag*eps)/((Xig+Zetag*eps) + delta*(g2+g3+g4) + eps*(g5n+g6n+g7n+g8n+g9n));
XO2 = (delta*g2 + eps*g5n)/((Xig+Zetag*eps) + delta*(g2+g3+g4) + eps*(g5n+g6n+g7n+g8n+g9n));
XN2 = (delta*g3 + eps*g8n)/((Xig+Zetag*eps) + delta*(g2+g3+g4) + eps*(g5n+g6n+g7n+g8n+g9n));
XAr = (delta*g4 + eps*g9n)/((Xig+Zetag*eps) + delta*(g2+g3+g4) + eps*(g5n+g6n+g7n+g8n+g9n));
XCO2 = eps*g6n/((Xig+Zetag*eps) + delta*(g2+g3+g4) + eps*(g5n+g6n+g7n+g8n+g9n));
XH2O = eps*g7n/((Xig+Zetag*eps) + delta*(g2+g3+g4) + eps*(g5n+g6n+g7n+g8n+g9n));
mmint = Xfueln*mmfuel + Xgasn*mmgas + XO2*MO2 + XN2*MN2 + XAr*MAr + XCO2*MCO2 + XH2O*MH2O;
Rint = Runiv/mmint;
TIVC = PIVC*VIVC/(mint*Rint);
% --- Fill output --- %
cr(7) = mres;               % Mass of residual per cycle [kg/cycle]
cr(8) = megr;               % Mass of egr per cycle [kg/cycle]
cr(9) = mint;               % Total mass at IVC [kg/cycle]
cr(10) = TIVC;              % Intake temperature at IVC [K]
cr(11) = Rint;              % Intake gas constant at IVC [K]
% Overall global reaction parameters
% A few of these are repeats, but I thought it would be wise to keep
% them all together in order
cr(20) = Xip;               % Moles of liquid fuel [mol/cycle]
cr(21) = Xig;               % Moles of gaseous fuel [mol/cycle]
cr(22) = Zetap;             % Moles of liquid fuel left after combustion [mol/cycle]
cr(23) = Zetag;             % Moles of gaseous fuel left after combustion [mol/cycle]
cr(24) = alp;               % Residual component [mol]
cr(25) = beta;              % EGR component [mol]
cr(26) = alp+beta;          % Epsilon component [mol]
cr(27) = delta;             % Moles of air [mol/cycle]
cr(28) = g2;                % Oxygen component in air [-]
cr(29) = g3;                % Nitrogen component in air [-]
cr(30) = g4;                % Argon component in air [-]
cr(31) = g5n;               % Oxygen component in products
cr(32) = g6n;               % Carbon dioxide component in products
cr(33) = g7n;               % Water component in products
cr(34) = g8n;               % Nitrogen component in products
cr(35) = g9n;               % Argon component in products

%% Local reaction for liquid fuel
% l3 = w;
% l4 = x/2;
% l2 = (2*l3 + l4 - y)/2;
% l5 = z/2;
% cr(40) = l2;
% cr(41) = l3;
% cr(42) = l4;
% cr(43) = l5;
% 
% % Local reaction for gaseous fuel
% if (iS(2) == 1)
%     n3 = a;
%     n4 = b/2;
%     n2 = (2*n3 + n4 - c)/2;
%     n5 = d/2;
%     cr(50) = n2;
%     cr(51) = n3;
%     cr(52) = n4;
%     cr(53) = n5;
% else
%     cr(50) = 0;
%     cr(51) = 0;
%     cr(52) = 0;
%     cr(53) = 0;
% end