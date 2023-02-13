function [data, imep] = ptraceeval(iS, datafile, inp, calc)
% Output
% data(:,1) = data points
% data(:,2) = time [s]
% data(:,3) = crank angle degrees [deg]
% data(:,4) = crank angle radians [rad]
% data(:,5) = volume [m3]
% data(:,6) = change in volume [m3]
% data(:,7) = dV/d(theta) [m3/rad]
% data(:,8) = surface area [m2]
% 10+ Pressure
% data(:,10)= TLA adjusted average pressure [Pa]
% data(:,11)= Spline smoothed pressure [Pa]
% data(:,12)= Filtered smoothed pressure [Pa]
% data(:,13)= Average smoothed pressure [Pa]
% data(:,14)= Matlab filtered TLA adjusted average pressure [Pa]
% 20+ First Derivative
% data(:,20)= Raw first derivative of pressure [Pa/deg] - uses data(:,13)
% data(:,21)= Spline fit [Pa/deg]
% data(:,22)= Filtered spline fit [Pa/deg]
% data(:,23)= Filtered/rolled fit [Pa/deg]
% data(:,24)= Matlab filtered data(:,20)
% 30+ Second Derivative
% data(:,30)= Raw second derivative of pressure [Pa/deg^2] - uses data(:,13)
% data(:,31)= Finite difference of first derivative results [Pa/deg^2] - uses data(:,23)
% Note - the use of data(:,31) for the second derivative provides better looking results
% data(:,32)= Spline fit [Pa/deg^2]
% data(:,33)= Filtered spline fit [Pa/deg^2]
% data(:,34)= Filtered/rolled fit [Pa/deg^2]
% data(:,35)= Matlab filtered data(:,30)

% --- Switches
% The user can add more higher order finite difference options if they wish
% https://ocw.snu.ac.kr/sites/default/files/NOTE/Lecture%2008_0.pdf
ifd(1) = iS(3);         % First-order finite difference: 1-central difference O(2); 2-forward diff O(1); 3-backward diff O(1)
ifd(2) = iS(4);         % Second-order finite difference: 1-central difference O(2); 2-forward diff O(1); 3-backward diff O(1)

% --- Input
nR = inp(7);            % four-stroke cycle [-]
dres = inp(62);         % degree resolution of data [deg]
Nrpm = inp(20);         % Engine speed [rev/min]
br = inp(1);            % Engine bore [m]
st = inp(2);            % Engine stroke [m]
rod = inp(3);           % Connecting rod length [m]
Vd = inp(4);            % Engine displacement [m3]
rc = inp(5);            % Compression ratio [-]
volbwl = inp(6);        % Piston bowl volume [m3]
TLAshift = inp(66);     % Thermodynamic Loss Angle [deg]
m = inp(61);            % Layers of pressure data [-]
cycle = calc(1);        % Total number of degrees 
n = calc(2);            % Total number of data points [-]
dps = calc(11);         % Degrees per second: [rev/min]*[360 deg/rev]*[1 min/60 sec]
rps = calc(24);         % Radians per second [rad/s] = [deg/s] * 3.14156 rad/180 deg
dt = calc(12);          % Time-step between data points [s]
tdcfit = calc(14);      % How many data points to shift TDC [-]
cl = inp(8);            % Crank length [m]
Vc = inp(9);            % Clearance volume [m3]
SApiston = inp(10);     % Piston surface area [m2]

% --- Matlab filtering parameters
% https://www.mathworks.com/help/matlab/data_analysis/filtering-data.html
% Create the filter coefficient vectors
af = 1;
bf = [1/4 1/4 1/4 1/4];

% --- Output matrix
s=0;                                    %creates sum variable, sets to 0
data=zeros(n,40);                       %creates output matrix

% --- Read the raw data
fid=fopen(datafile,'r');
a=textscan(fid,'%f');
dat=a{1};               % array of raw pressure data: 1 to n (crank angle)
fclose(fid);

% --- Preliminary data and filtering
for i=1:n                               %runs for all datapoints
    data(i,1)=i;                        %indexes run datapoint number
    data(i,2)=dt*(i-1);                 %indexes time
    data(i,3)=dat((m+1)*(i-1)+1);       %indexes current angle, in degrees
    data(i,4)=data(i,3)*pi()/180;       %indexes current angle, in radians
    pos=cl*cos(data(i,4))+(rod^2-(cl*sin(data(i,4)))^2)^(1/2);
    %calculates current volume within cylinder
    data(i,5)=Vc+volbwl+(pi()*br^2)/4*(rod+cl-pos);
    if i==1
        data(i,6)=0;                    %first change in volume is 0
     else
        data(i,6)=data(i,5)-data(i-1,5);%calculates change in volume [m3]
    end
    % calculate dV/d(theta) - m3/rad
    data(i,7)=(pi()*br*br*cl*cl*sin(2*data(i,4)))/(8*sqrt(rod*rod - cl*cl*sin(data(i,4))*sin(data(i,4)))) + pi()*br*br*cl*sin(data(i,4))/4;
    %calculate clearance volume height
    hclear=Vc/(pi()*br^2/4);
    %calculate clearance surface area (head of cylinder included)
    SAclear=pi()*br^2/4+pi()*br*hclear;
    %calculate piston sweep surface area
    SAmoving=pi()*br*(rod+cl-pos);
    %calculate total surface area
    data(i,8)=SApiston+SAclear+SAmoving;
    %sum pressures for a point
    for j=1:m
        s=s+dat((m+1)*(i-1)+1+j);       %sums pressures for a data point
    end
    %calculates average of m pressures at a data point, accounts for TLA
    %pressures are converted to Pa
    if i+tdcfit<1                       %rolls first pressures to end
        data(n-(i+tdcfit-1),10)=s/m*100000;
    elseif i+tdcfit>n                   %rolls last pressures to beginning
        data(i+tdcfit-n,10)=s/m*100000;
    else                                %shifts pressures by TLA
        data(i+tdcfit,10)=s/m*100000;
    end
    s=0;                                %resets sum to 0 for next step
    % Note: there was an option here to check for too high of pressures
    for j=1:m
        if dat((m+1)*(i-1)+1+j)<=0
            %removes negative pressures in case there are any
            dat((m+1)*(i-1)+1+j)=dat((m+1)*(i-1)+j);
        else
            s=s+dat((m+1)*(i-1)+1+j);   %sums pressure readings at a point
        end
    end
    if i+tdcfit<1                       %rolls first pressures to end
        data(n-(i+tdcfit-1),11)=s/m*100000;
    elseif i+tdcfit>n                   %rolls last pressures to beginning
        data(i+tdcfit-n,11)=s/m*100000;
    else                                %shifts pressures by TLA
        data(i+tdcfit,11)=s/m*100000;
    end
    s=0;                                %resets sum to 0 for next step
end

%NOTE: At this point, we only have the average crank angle pressure data
%that has been adjusted according to TLA. No filtering or smoothing has
%been accomplished yet

%% --- Mattson Filtered pressure calculation
for i=1:n
    if i<3
        data(i,12)=data(i,11);            %first pressures are unfiltered
    elseif i>(n-2)
        data(i,12)=data(i,11);            %last pressures are unfiltered
    else
        %calculates local average value
        lav=(data(i-2,11)+data(i-1,11)+data(i+1,11)+data(i+2,11))/4;
        %if pressure is too large or too small, reset to local average
        if data(i,11)>lav*1.001
            data(i,12)=lav;
        elseif data(i,11)<lav*0.999
            data(i,12)=lav;
        else
            data(i,12)=data(i,11);        %if in range, set to previous value
        end            
    end
end

% --- Rolling average pressure calculation
for i=1:n
    if i<3                              %first pressures are unfiltered
        data(i,13)=data(i,12);
    elseif i>(n-2)                      %last pressures are unfiltered
        data(i,13)=data(i,12);
    else
        %averages filtered values for a local value
        data(i,13)=(data(i-2,12)+data(i-1,12)+data(i,12)+data(i+1,12)+data(i+2,12))/5;
    end
end

% Filter the average pressure data that was shifted by TLA
data(:,14) = filter(bf,af,data(:,10));

%% --- Change in pressure with respect to crank angle
for i=1:n
    if i==1
        %first change in pressure is set to 0
        data(i,20)=0;
        data(i,30)=0;
    else
        %calculates 1st derivative of pressure (Pa/deg)
        % Euler explicit
        if (ifd(1) == 1) % central difference
            if (i ~= n)
                data(i,20) = (data(i+1,13) - data(i-1,13))/(2*dres);
            else
                data(i,20) = data(i-1,20);
            end
        elseif (ifd(1) == 2) % forward difference
            data(i,20)=(data(i,13)-data(i-1,13))/dres;
        else                % backward difference
            if (i ~= n)
                data(i,20)=(data(i+1,13)-data(i,13))/dres;
            else 
                data(i,20)=data(i-1,20);
            end
        end
        % Calculates 2nd derivative of pressure (Pa/deg^2)
        if (ifd(2) == 1)        % Central difference
            if (i ~= n)
                data(i,30)=(data(i+1,13) - 2*data(i,13) + data(i-1,13))/(dres*dres);
            else
                data(i,30) = data(i-1,30);
            end
        elseif (ifd(2) == 2)    % forward difference
            if (i < n-2)
                data(i,30) = (data(i,13) - 2*data(i+1,13) + data(i+2,13))/(dres*dres);
            else
                data(i,30) = data(i-1,30);
            end
        else            % backward difference
            if (i > 2)
                data(i,30) = (data(i-2,13) - 2*data(i-1,13) + data(i,13))/(dres*dres);
            else
                data(i,30) = data(i-1,30);
            end
        end
    end
end

% --- Pressure filtering operation, 1st derivative of pressure
cspline = spline(data(:,3),data(:,20));
for i=1:n
    data(i,21)=ppval(cspline,data(i,3));       %solves for spline pressure
    if i<3
        data(i,22)=data(i,21);                %first pressures are unfiltered
    elseif i>(n-2)
        data(i,22)=data(i,21);                %last pressures are unfiltered
    else
        %calculates local average value
        lav=(data(i-2,21)+data(i-1,21)+data(i+1,21)+data(i+2,21))/4;
        %if pressure is too large or too small, reset to local average
        if data(i,3)<=-60                   %values before combustion
            if data(i,21)>lav*1.01
                data(i,22)=lav;
            elseif data(i,21)<lav*0.99
                data(i,22)=lav;
            else
                data(i,22)=data(i,1);        %if in range, set to previous value
            end
        elseif data(i,3)>=150               %values after combustion
            if data(i,21)>lav*1.01
                data(i,22)=lav;
            elseif data(i,21)<lav*0.99
                data(i,22)=lav;
            else
                data(i,22)=data(i,21);        %if in range, set to previous value
            end
        else                            %combustion values are unfiltered
            data(i,22)=data(i,21);
        end
    end
end
for i=1:n
    if i<3                              %first pressures are unfiltered
        data(i,23)=data(i,22);
    elseif i>(n-2)                      %last pressures are unfiltered
        data(i,23)=data(i,22);
    else
        %averages filtered values for a local value
        data(i,23)=(data(i-2,22)+data(i-1,22)+data(i,22)+data(i+1,22)+data(i+2,22))/5;
    end
end
% Matlab filtered
data(:,24) = filter(bf,af,data(:,20));

% --- Second derivatives
for i=3:n-1
    if (ifd(1) == 1) % central difference
        if (i ~= n) 
            data(i,31) = (data(i+1,23) - data(i-1,23))/(2*dres);
        else
            data(i,31) = data(i-1,31);
        end
    elseif (ifd(1) == 2) % forward difference
        data(i,31) = (data(i+1,23)-data(i,23))/dres;
    else
        if (i ~= n)
            data(i,31)=(data(i+1,23)-data(i,23))/dres;
        else
            data(i,31)=data(i-1,31);
        end
    end
end

cspline2=spline(data(:,3),data(:,31));            %calculates spline equation

for i=1:n
    data(i,32)=ppval(cspline2,data(i,3));
    if i<3
        data(i,33)=data(i,32);                %first pressures are unfiltered
    elseif i>(n-2)
        data(i,33)=data(i,32);                %last pressures are unfiltered
    else
        %calculates local average value
        lav=(data(i-2,32)+data(i-1,32)+data(i+1,32)+data(i+2,32))/4;
        %if pressure is too large or too small, reset to local average
        if data(i,3)<=-60                   %values before combustion
            if data(i,32)>lav*1.2
                data(i,33)=lav;
            elseif data(i,32)<lav*0.8
                data(i,33)=lav;
            else
                data(i,33)=data(i,32);        %if in range, set to previous value
            end   
        elseif data(i,3)>=150               %values after combustion
            if data(i,32)>lav*1.2
                data(i,33)=lav;
            elseif data(i,32)<lav*0.8
                data(i,33)=lav;
            else
                data(i,33)=data(i,32);        %if in range, set to previous value
            end     
        else                            %combustion values are unfiltered
            data(i,33)=data(i,32);
        end
    end
end
for i=1:n
    if i<3                              %first pressures are unfiltered
        data(i,34)=data(i,33);
    elseif i>(n-2)                      %last pressures are unfiltered
        data(i,34)=data(i,33);
    else
        %averages filtered values for a local value
        data(i,34)=(data(i-2,33)+data(i-1,33)+data(i,33)+data(i+1,33)+data(i+2,33))/5;
    end
end
% Matlab filtered
data(:,35) = filter(bf,af,data(:,30));

%% --- Third derivative of pressure
% There are a few papers out there that indicate the third deriviative of
% pressure might be helpful. That has not been my experience; however, feel
% free to add it to the code and start indexing it at 40+
% https://www.sae.org/publications/technical-papers/content/861216/
% https://www.sae.org/publications/technical-papers/content/920808/

%% --- IMEP,net calculations
dataD = load(datafile);     % Will load it as a data array
% imepn = Win/Vd
% Using the smoothed pressure
Win = 0;
for i=2:n
    Win = Win + (data(i,13)+data(i-1,13))*data(i,7)*rps*dt;     % [Pa]*[m3/rad]*[rad/s]*[s] = [J]
end
imep(1) = Win/Vd; % [J]/[m3] = [Pa]
% Find the imepn of each pressure trace
for j=1:m       % over the layers of pressure data
    Winm(j) = 0;
    for i=2:n
        Winm(j) = Winm(j) + (dataD(i,j+1)+dataD(i-1,j+1))*100000*data(i,7)*rps*dt;
    end
    imep(1+j) = Winm(j)/Vd; % imep,net [Pa] of the raw data
end
imepnm = Winm/Vd;                       % imep,net [Pa] of the raw data
%imepnmean = mean(imepnm);              % Mean imep,net [Pa]
imep(m+2) = mean(imepnm);               % Mean imep,net [Pa]
%imepnS = std(imepnm);                  % Standard deviation of imep [Pa]
imep(m+3) = std(imepnm);                % Standard deviation of imep [Pa]
%imepCOV = imepnS/imepnmean;            % Coefficient of variation [-]
imep(m+4) = imep(m+3)/imep(m+2);        % Coefficient of variation [-]