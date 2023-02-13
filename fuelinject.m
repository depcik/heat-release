function [fuelinj, calc] = fuelinject(iS, inp, calc, fuel, cr, pdata)

% iS(1) - number of fuel injection events
% calc(2) - total number of data points
% calc(12) - time-step of data [sec]
% inp(23) - first injection mass flow rate [kg/s]
% cr(3) - total liquid fuel mass per cycle [kg/cycle]
% Our model will inject fuel until we reach the total liquid fuel mass. If
% there are multiple injections, each injection will end when the mass of
% fuel injected during that event has been reached.
% calc(5) - index of INJ1
% calc(6) - index of INJ2
% calc(7) - index of INJ3
% calc(8) - index of INJ4
% calc(9) - index of INJ5
% inp(72) - constant in fuel injection model
% fuel(8) - density of liquid fuel [kg/m3]
% inp(24) - pressure of fuel injector [Pa]
% inp(42) - number of fuel injectors [-]
% inp(43) - number of fuel injector holes
% inp(45) - individual injector nozzle area [m2]
% pdata(:,13) - Average smoothed pressure [Pa]

% We will create a fuel injection profile that fills the entire number of
% data points
fuelinj=zeros(calc(2),1);     % Mass flow rate of fuel at each point [kg/s]

% --- Local variables --- %
Cd = inp(72);
rhof = fuel(8);
ninj = inp(42);
nh = inp(43);
An = inp(45);
dt = calc(12);  
evap = inp(87);
% Find the total mass injected in each injection event
mtot1 = inp(23)/calc(13);   % First injection event [kg]
mtot2 = inp(27)/calc(13);   % Second injection event [kg]
mtot3 = inp(31)/calc(13);   % Third injection event [kg]
mtot4 = inp(35)/calc(13);   % Fourth injection event [kg]
mtot5 = inp(39)/calc(13);   % Fifth injection event [kg]
% Find the injection pressures in each injection event
pnoz1 = inp(24);    % [Pa]
pnoz2 = inp(28);    % [Pa]
pnoz3 = inp(32);    % [Pa]
pnoz4 = inp(36);    % [Pa]
pnoz5 = inp(40);    % [Pa]
% Indexes of injection events
INJ1 = calc(5);
INJ2 = calc(6);
INJ3 = calc(7);
INJ4 = calc(8);
INJ5 = calc(9);

% --- Create injection profile for first injection
mtot1c = 0;                                         % Check of total mass of fuel injected in first injection [kg]
i = INJ1;                                           % When the first injection starts
iEND1 = 0;
while ((mtot1c < mtot1) && (iEND1 == 0))
    pcyl = pdata(i,13);                             % Pressure in the cylinder at the time injection starts [Pa]
    dp = (pnoz1-pcyl);                               % Pressure difference seen by injector and cylinder [Pa]
    fuelinj(i) = fuelinj(i) + evap*Cd*An*nh*ninj*sqrt(2*rhof*dp);     % Mass flow rate of fuel in [kg/s]
    dm = fuelinj(i)*dt;                                        % Mass added during time-step [kg]
    % Check to make sure does not exceed total amount
    if ((mtot1c+dm) > mtot1)
        dm = mtot1-mtot1c;                              % Mass added during last time-step [kg]
        iEND1 = 1;
        fuelinj(i) = dm/dt;                         % Mass flow rate during last time-step [kg/s]
    end
    mtot1c = mtot1c + dm;
    i = i+1;
end
% Find the end of the injection process
INJ1END = i-1; 
INJ1deg = (INJ1END-1)*inp(62) - 180;            % Crank angle degree
% Fill output
calc(30) = INJ1END;
calc(31) = INJ1deg;

% --- Create injection profile for second injection
if (iS(1) > 1)
    mtot2c = 0;                                         % Check of total mass of fuel injected in first injection [kg]
    i = INJ2;                                           % When the first injection starts
    iEND2 = 0;
    while ((mtot2c < mtot2) && (iEND2 == 0))
        pcyl = pdata(i,13);                             % Pressure in the cylinder at the time injection starts [Pa]
        dp = (pnoz2-pcyl);                               % Pressure difference seen by injector and cylinder [Pa]
        fuelinj(i) = fuelinj(i) + evap*Cd*An*nh*ninj*sqrt(2*rhof*dp);     % Mass flow rate of fuel in [kg/s]
        dm = fuelinj(i)*dt;                                        % Mass added during time-step [kg]
        % Check to make sure does not exceed total amount
        if ((mtot2c+dm) > mtot2)
            dm = mtot2-mtot2c;                              % Mass added during last time-step [kg]
            iEND2 = 1;
            fuelinj(i) = dm/dt;                         % Mass flow rate during last time-step [kg/s]
        end
        mtot2c = mtot2c + dm;
        i = i+1;
    end
    % Find the end of the injection process
    INJ2END = i-1;
    INJ2deg = (INJ2END-1)*inp(62) - 180;            % Crank angle degree
    % Fill output
    calc(32) = INJ2END;
    calc(33) = INJ2deg;
end

% --- Create injection profile for third injection
if (iS(1) > 2)
    mtot3c = 0;                                         % Check of total mass of fuel injected in first injection [kg]
    i = INJ3;                                           % When the first injection starts
    iEND3 = 0;
    while ((mtot3c < mtot3) && (iEND3 == 0))
        pcyl = pdata(i,13);                             % Pressure in the cylinder at the time injection starts [Pa]
        dp = (pnoz3-pcyl);                               % Pressure difference seen by injector and cylinder [Pa]
        fuelinj(i) = fuelinj(i) + evap*Cd*An*nh*ninj*sqrt(2*rhof*dp);     % Mass flow rate of fuel in [kg/s]
        dm = fuelinj(i)*dt;                                        % Mass added during time-step [kg]
        % Check to make sure does not exceed total amount
        if ((mtot3c+dm) > mtot3)
            dm = mtot3-mtot3c;                              % Mass added during last time-step [kg]
            iEND3 = 1;
            fuelinj(i) = dm/dt;                         % Mass flow rate during last time-step [kg/s]
        end
        mtot3c = mtot3c + dm;
        i = i+1;
    end
    % Find the end of the injection process
    INJ3END = i-1;
    INJ3deg = (INJ3END-1)*inp(62) - 180;            % Crank angle degree
    % Fill output
    calc(34) = INJ3END;
    calc(35) = INJ3deg;
end

% --- Create injection profile for fourth injection
if (iS(1) > 3)
    mtot4c = 0;                                         % Check of total mass of fuel injected in first injection [kg]
    i = INJ4;                                           % When the first injection starts
    iEND4 = 0;
    while ((mtot4c < mtot4) && (iEND4 == 0))
        pcyl = pdata(i,13);                             % Pressure in the cylinder at the time injection starts [Pa]
        dp = (pnoz4-pcyl);                               % Pressure difference seen by injector and cylinder [Pa]
        fuelinj(i) = fuelinj(i) + evap*Cd*An*nh*ninj*sqrt(2*rhof*dp);     % Mass flow rate of fuel in [kg/s]
        dm = fuelinj(i)*dt;                                        % Mass added during time-step [kg]
        % Check to make sure does not exceed total amount
        if ((mtot4c+dm) > mtot4)
            dm = mtot4-mtot4c;                              % Mass added during last time-step [kg]
            iEND4 = 1;
            fuelinj(i) = dm/dt;                         % Mass flow rate during last time-step [kg/s]
        end
        mtot4c = mtot4c + dm;
        i = i+1;
    end
    % Find the end of the injection process
    INJ4END = i-1;
    INJ4deg = (INJ4END-1)*inp(62) - 180;            % Crank angle degree
    % Fill output
    calc(36) = INJ4END;
    calc(37) = INJ4deg;
end

% --- Create injection profile for fifth injection
if (iS(1) > 4)
    mtot5c = 0;                                         % Check of total mass of fuel injected in first injection [kg]
    i = INJ5;                                           % When the first injection starts
    iEND5 = 0;
    while ((mtot5c < mtot5) && (iEND5 == 0))
        pcyl = pdata(i,13);                             % Pressure in the cylinder at the time injection starts [Pa]
        dp = (pnoz5-pcyl);                               % Pressure difference seen by injector and cylinder [Pa]
        fuelinj(i) = fuelinj(i) + evap*Cd*An*nh*ninj*sqrt(2*rhof*dp);     % Mass flow rate of fuel in [kg/s]
        dm = fuelinj(i)*dt;                                        % Mass added during time-step [kg]
        % Check to make sure does not exceed total amount
        if ((mtot5c+dm) > mtot5)
            dm = mtot5-mtot5c;                              % Mass added during last time-step [kg]
            iEND5 = 1;
            fuelinj(i) = dm/dt;                         % Mass flow rate during last time-step [kg/s]
        end
        mtot5c = mtot5c + dm;
        i = i+1;
    end
    % Find the end of the injection process
    INJ5END = i-1;
    INJ5deg = (INJ5END-1)*inp(62) - 180;            % Crank angle degree
    % Fill output
    calc(38) = INJ5END;
    calc(39) = INJ5deg;
end