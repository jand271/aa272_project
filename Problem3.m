clc
clear
close all

%% Problem 1
xdiff = 1:20;
ydiff = 1:20;
zdiff = 1:20;
bdiff = 1:20;
count = 0;
for rownum = 1:20
    count = count + 1;
%rownum = 1;
gnss_t = readtable('gnss_log');
gnss = gnss_t{rownum,:};

eph_t = readtable('ephem');
[Xecef,Yecef,Zecef,bias] = sat_pos_ecef(eph_t,gnss);
xdiff(count) = gnss(22) - Xecef;
ydiff (count)= gnss(23) - Yecef;
zdiff(count) = gnss(24) - Zecef;
bdiff(count) = gnss(25) - bias;
end
figure
hold on
plot(xdiff)
plot(ydiff)
plot(zdiff)
plot(bdiff)
xlabel('sequance')
ylabel('difference between measured ECEF satellite and calculated (m)')
legend('xecef difference','yecef difference','zecef difference','bias difference')
title('check for the difference between calculated and provided  XYZB')
%% Problem 2
eph_t = readtable('ephem');
gnss_t = readtable('gnss_log');
Data = gnss_t{:,:};
Xecef = zeros(size(Data,1),1);
Yecef = zeros(size(Data,1),1);
Zecef = zeros(size(Data,1),1);
Bias = zeros(size(Data,1),1);
for i = 1:size(Data,1)
    rownum = i;
    gnss = gnss_t{rownum,:};
    [xecef,yecef,zecef,bias] = sat_pos_ecef(eph_t,gnss);
    Xecef(i) = xecef;
    Yecef(i) = yecef;
    Zecef(i) = zecef;
    Bias(i) = bias;
end
T = readtable('gnss_log');
Data = T{:,:};
Week_s = 604800;
TimeNanos = Data(:,3);
FullBiasNanos = Data(:,4);
TimeOffsetNanos = Data(:,6);
BiasNanos = Data(:,19);
ReceivedSvTimeNanos = Data(:,8);
N_w = floor(-FullBiasNanos/604800e9);
t_Rx_hw = TimeNanos+ TimeOffsetNanos;
b_hw = FullBiasNanos+ BiasNanos;
t_Rx_GPS = t_Rx_hw - b_hw;
t_Rx_w = t_Rx_GPS - N_w.*(604800e9);
rho_ns = t_Rx_w - ReceivedSvTimeNanos;

rho_m = rho_ns .* 299792458/1e9;

SV_ID = Data(:,5);

info = [Data(:,2) Data(:,5) Xecef Yecef Zecef Bias rho_m];
%% Problem II
%Implement a function that takes satellite positions (X_k) and 
%current position (X0) estimate asinputs and outputs the expected pseudorange.


% list all the time entries
time_inst = info(:,1);
time_inst = unique(time_inst);
time_inst = sort(time_inst);

% set up location array
X_Location_new = zeros(length(time_inst),1);
Y_Location_new = zeros(length(time_inst),1);
Z_Location_new = zeros(length(time_inst),1);
Radius_Check = zeros(length(time_inst),1);
count = 0;
time_inst = time_inst';
for time = time_inst
    count = count + 1;
    
    
    % initialize the user position and user time bias
    X0 = [0,0,0];
    b0 = 0;
    % find the measurements that correspond to the corresponding time
    % instance
    ind = find(info(:,1) == time);
    %slice the matrix for this time
    measure = info(ind,:);
    % check for svid uniqueness
    svid_check = unique(measure(:,2));
    x_location = nan;
    y_location = nan;
    z_location = nan;
    radius_check = nan;
    if length(svid_check)>=4
       % disp('check')
        % satellite position
        % [x1 y1 z1;x2 y2 z2;x3 y3 z3]
        X_k = measure(:,3:5);
        % satellite bias
        % [B1;B2;B3]
        B_k = measure(:,6);
        % [rho1;rho2;rho3]
        rho = measure(:,7);
        del_X = inf;
        del_b = inf;
        while (abs(del_X(1)) > 0.00001) || (abs(del_X(2)) > 0.00001)||(abs(del_X(2)) > 0.00001)||(abs(del_b) > 0.00001)
            %norm(del_X)>0.01 || del_b>0.01
            G = geometric_matrix(X_k,X0);
            rho_0 = expected_psedorange(X_k, X0, b0, B_k);
            %del_rho = rho - rho_0;
            del_rho = -(rho_0 - rho);
            %del_u = [del_x;del_bu]
            del_u = inv(G'*G)*G'*del_rho;

            del_X = del_u(1:3)';
            del_b = del_u(4);
            X0 = X0 + del_X; b0 = b0 + del_b;
        end
        x_location = X0(1);
        y_location = X0(2);
        z_location = X0(3);
        radius_check = norm([x_location,y_location,z_location]);
    end
    X_Location_new(count) = x_location;
    Y_Location_new(count) = y_location;
    Z_Location_new(count) = z_location;
    Radius_Check(count) = radius_check;
    
end

info = [Data(:,2) Data(:,5) Data(:,22) Data(:,23) Data(:,24) Data(:,25) rho_m];
%% Problem II
%Implement a function that takes satellite positions (X_k) and 
%current position (X0) estimate asinputs and outputs the expected pseudorange.


% list all the time entries
time_inst = info(:,1);
time_inst = unique(time_inst);
time_inst = sort(time_inst);

% set up location array
X_Location = zeros(length(time_inst),1);
Y_Location = zeros(length(time_inst),1);
Z_Location = zeros(length(time_inst),1);
Radius_Check = zeros(length(time_inst),1);
count = 0;
time_inst = time_inst';
for time = time_inst
    count = count + 1;
    
    
    % initialize the user position and user time bias
    X0 = [0,0,0];
    b0 = 0;
    % find the measurements that correspond to the corresponding time
    % instance
    ind = find(info(:,1) == time);
    %slice the matrix for this time
    measure = info(ind,:);
    % check for svid uniqueness
    svid_check = unique(measure(:,2));
    x_location = nan;
    y_location = nan;
    z_location = nan;
    radius_check = nan;
    if length(svid_check)>=4
       % disp('check')
        % satellite position
        % [x1 y1 z1;x2 y2 z2;x3 y3 z3]
        X_k = measure(:,3:5);
        % satellite bias
        % [B1;B2;B3]
        B_k = measure(:,6);
        % [rho1;rho2;rho3]
        rho = measure(:,7);
        del_X = inf;
        del_b = inf;
        while (abs(del_X(1)) > 0.00001) || (abs(del_X(2)) > 0.00001)||(abs(del_X(2)) > 0.00001)||(abs(del_b) > 0.00001)
            %norm(del_X)>0.01 || del_b>0.01
            G = geometric_matrix(X_k,X0);
            rho_0 = expected_psedorange(X_k, X0, b0, B_k);
            %del_rho = rho - rho_0;
            del_rho = -(rho_0 - rho);
            %del_u = [del_x;del_bu]
            del_u = inv(G'*G)*G'*del_rho;

            del_X = del_u(1:3)';
            del_b = del_u(4);
            X0 = X0 + del_X; b0 = b0 + del_b;
        end
        x_location = X0(1);
        y_location = X0(2);
        z_location = X0(3);
        radius_check = norm([x_location,y_location,z_location]);
    end
    X_Location(count) = x_location;
    Y_Location(count) = y_location;
    Z_Location(count) = z_location;
    Radius_Check(count) = radius_check;
    
end



figure('Renderer','opengl')
hold on
plot3(X_Location_new,Y_Location_new,Z_Location_new,'-')
plot3(X_Location,Y_Location,Z_Location,'-')
view(45,45)
legend('new solution','old solution')
xlabel('X_Location (m)')
ylabel('Y_Location (m)')
zlabel('Z_Location (m)')
title('Trajectory of the Solution')

%% Problem 3
gnss = gnss_t{end,:};
svid = gnss(5);
FullBiasNanos = gnss(4);
Measurement_ReceivedSvTimeNanos = gnss(8);
eph = table2cell(eph_t);
Ephem_GPSWeek = eph{svid,23};
Ephem_Toe = eph{svid,13};
% a

receiver_position = [X_Location_new(end),Y_Location_new(end),Z_Location_new(end)];

visiable_list = [];
for sat = 1:size(eph_t,1)
    
    [Xecef,Yecef,Zecef,~] = sat_pos_ecef_sv(eph_t,sat,FullBiasNanos,Measurement_ReceivedSvTimeNanos,0);
    
%     gnss = gnss_t{end-2,:};
%      [Xecef,Yecef,Zecef,bias] = sat_pos_ecef(eph_t,gnss)
    satellite_position = [Xecef,Yecef,Zecef];
   
     diff = (satellite_position - receiver_position)/norm(satellite_position - receiver_position);

    receiver_dot_vec = -receiver_position/norm(receiver_position);

   if dot(receiver_dot_vec,diff)<0
        visiable_list = [visiable_list sat];
    end
   
end
visiable_satellites = visiable_list
%b


% receiver position
visiable_list_1hr = []; 
for sat = 1:size(eph_t,1)
    
    [Xecef,Yecef,Zecef,~] = sat_pos_ecef_sv(eph_t,sat,FullBiasNanos,Measurement_ReceivedSvTimeNanos,3600);
    
%     gnss = gnss_t{end-2,:};
%      [Xecef,Yecef,Zecef,bias] = sat_pos_ecef(eph_t,gnss)
    satellite_position = [Xecef,Yecef,Zecef];
   
     diff = (satellite_position - receiver_position)/norm(satellite_position - receiver_position);

    receiver_dot_vec = -receiver_position/norm(receiver_position);

   if dot(receiver_dot_vec,diff)<0
        visiable_list_1hr = [visiable_list_1hr sat];
    end
   
end
visiable_satellites_1hr = visiable_list_1hr
%%

function r_enu = ecef2enu(lam_gs , phi_gs , r_ecef)
% ecef2enu  Converts  ECEF  position  to ENU  position
%
%    Note: This  function  assumes  an  elliptical  Earth
%
% Inputs:
%    lam_gs  - Ground  station  Geodetic  longitude [rad]
%    phi_gs  - Ground  station  Geodetic  latitude [rad]
%    r_ecef  - Satellite  ECEF  position  vector [km]
%
% Outputs:
%     r_enu  - Satellite  ENU  position  vector [km]
% Earth  parameters
rE = 6378.1;eE = 0.0818;
% Compute  the  satellite  position  wrt  the  ground  station  inECEF
N = rE/sqrt(1 - eE^2*( sin(phi_gs)^2));
r_GS_ecef = [N*cos(phi_gs)*cos(lam_gs);
    N*cos(phi_gs)*sin(lam_gs);
    N*(1-eE ^2)*sin(phi_gs)];
r_sat_GS = r_ecef  - r_GS_ecef;
% Develop  the  ECEF to ENU  rotation  matrix
E = [-sin(lam_gs);cos(lam_gs);0];
N = [-sin(phi_gs)*cos(lam_gs);
    -sin(phi_gs)*sin(lam_gs);
    cos(phi_gs)];
U = [cos(phi_gs)*cos(lam_gs);
    cos(phi_gs)*sin(lam_gs);
    sin(phi_gs)];
Recef2enu = [E, N, U].';
% Rotate  the  position  vector  into  ENU  coordinates
r_enu = Recef2enu*r_sat_GS;
end


function [lambda , phi , h] = ecef2geod(x, y, z, tol)
% ecef2geod  Converts  ECEF  coordinates  to  Geodetic  latitudeand  longitude
%
% Inputs:
%          x - ECEF x coordinate [km]
%          y - ECEF y coordinate [km]%          z - ECEF z coordinate [km]
%       tol - tolerance  of  coordinate  result [deg]
%
% Outputs:
%    lambda  - Geodetic  longitude  of  orbit [deg]
%       phi - Geodetic  latitude  of  orbit [deg]
%          h - Geodetic  altitude  of  orbit [km]
rE = 6378.1;
eE = 0.0818;
r = sqrt(x.^2 + y.^2 + z.^2);
rxy = sqrt(x.^2 + y.^2);
lambda = atan2d(y, x);
phi0 = asind(z./r);
N = rE/sqrt(1 - eE^2*( sind(phi0)^2));
diff = 2*tol;
while  abs(diff) > tol
    phi1 = atan2d ((z + N*eE^2* sind(phi0)), rxy);
    diff = phi1 - phi0;
    phi0 = phi1;
    N = rE/sqrt(1 - eE^2*( sind(phi0)^2));
end
phi = phi0;
h = rxy/cosd(phi) - N;
end


function rho_0 = expected_psedorange(X_k, X0, b0, B_k)

% X_k = [x1 y1 z1;x2 y2 z2] for satellite
% X0 = [x_user y_user z_user]
% B_k = [B1;B2;B3]
% rho_0 will be for all satellites in a vertical vector form

rho_0 = vecnorm(X_k - X0,2,2) + b0 - B_k;

end

%Implement a function that takes satellite positions (X_k)
%and the current position estimateas (X0) input and outputs the geometry matrix.  
%Remember that the first 3 columns of thegeometry matrix are a unit vector.

function G = geometric_matrix(X_k,X0)

num = X_k - X0;
deo = vecnorm(X_k - X0,2,2);
LOS = num./deo;
rownum = size(X_k,1);
One = ones(rownum,1);
G = [-LOS One];
end


function [Xecef,Yecef,Zecef,bias] = sat_pos_ecef(eph_t,gnss)
    
    svid = gnss(5);
    eph = table2cell(eph_t);
    bias = eph{svid,2}*299792458;
    
    sqrtA = eph{svid,12};

    %1
    a = sqrtA^2;
    mu = 3.986005e14;
    Omega_dot_E = 7.2921151467e-5;
    %2
    del_n  = eph{svid,7};
    n = sqrt(mu/(a^3)) + del_n;
    

    FullBiasNanos = gnss(4);
    N_w = floor(-FullBiasNanos/604800e9);
    Measurement_ReceivedSvTimeNanos = gnss(8);
    Measurement_Week = N_w;
    Ephem_GPSWeek = eph{svid,23};
    Ephem_Toe = eph{svid,13};
    %3
    tk = (Measurement_Week - Ephem_GPSWeek) *604800+...
        (Measurement_ReceivedSvTimeNanos/1e9 - Ephem_Toe);
    %4
    M0 = eph{svid,8};
    Mk = M0 + n*tk;

    
    %5
    e = eph{svid,10};
    tol = 1e-10;
    
     Ek = M2E (Mk , e , tol );
     Ek = wrapTo2Pi(Ek);
%      Mk = Mk
%      test = Ek - e*sin(Ek)

    %6
    sinnuk = sqrt(1 - e^2)*sin(Ek) / (1 - e*cos(Ek));
    cosnuk = (cos(Ek) - e) / (1 - e*cos(Ek));
    nuk = atan2(sinnuk, cosnuk);
    nuk = wrapTo2Pi(nuk);
    %7
    omega = eph{svid,19};
    phik = nuk + omega;
    phik = wrapTo2Pi(phik);
    
    %8
    Cus = eph{svid,11};
    Cuc = eph{svid,9};
    
    phik_temp = phik;
    for j = 1:100
        del_phik = Cus*sin(2*phik_temp) + Cuc* cos(2*phik_temp);

        %9
        uk  = phik + del_phik;
        phik_temp = uk;
    end
    
    %10
    Crs = eph{svid,6};
    
    Crc = eph{svid,18};
    
    del_rk = Crs * sin(2*phik) + Crc * cos(2*phik);
    
    %11
    Cis = eph{svid,16};
    
    Cic = eph{svid,14};
    
    del_ik = Cis * sin(2*phik) + Cic * cos(2*phik);
    
    %12
    Omega_0 = eph{svid,15};
    Omega_dot = eph{svid,20};
    t = (Measurement_Week - Ephem_GPSWeek) *604800+...
        Measurement_ReceivedSvTimeNanos/1e9;
    Omega_k = Omega_0 - Omega_dot_E * t + Omega_dot * tk;
    
    %13
    rk  = a *(1 - e*cos(Ek)) + del_rk;
    %14
    i0 = eph{svid,17};
    
    i_dot = eph{svid,21};
    
    ik = i0 + i_dot * tk + del_ik;
    
    %15
    xp = rk * cos(uk);
    
    %16
    
    yp = rk * sin(uk);
    
    %17
    
    Xecef = xp * cos(Omega_k) - yp *cos(ik) *sin(Omega_k);
    Yecef = xp * sin(Omega_k) + yp *cos(ik) *cos(Omega_k);
    Zecef = yp *sin(ik);
    
    
    
    function E = M2E (M , e , tol )
        % M2E solves Kepler 's equation for eccentric anomaly
        %
        % Note : This function uses a Newton - Raphson method to
        % numerically compute
        % the correct value for E
        % %
        % Inputs :
        % M - mean anomaly [ rad ]
        % e - eccentricity of orbit
        % tol - tolerance for Newton - Raphson iterator
        %
        % Outputs :
        % E - eccentric anomaly [ rad ]
        if M == 0 || M == pi
            Return the known solutions ( trivial )
            E = M ;
        else
            E0 = M;
            d = -( E0 - e * sin ( E0 ) - M ) /(1 - e * cos ( E0 ) ) ;
            
            % Loop until the solution converges
            while abs ( d ) > tol
                E1 = E0 + d ;
                d = -( E1 - e * sin ( E1 ) - M ) /(1 - e * cos ( E1 ) ) ;
                E0 = E1 ;
            end
            E = E0 ;
         end
    end
end
    
    
    
    
    
    
    



function [Xecef,Yecef,Zecef,bias] = sat_pos_ecef_sv(eph_t,svid,FullBiasNanos,Measurement_ReceivedSvTimeNanos,add_time)
    
   
    eph = table2cell(eph_t);
    bias = eph{svid,2}*299792458;
    
    sqrtA = eph{svid,12};

    %1
    a = sqrtA^2;
    mu = 3.986005e14;
    Omega_dot_E = 7.2921151467e-5;
    %2
    del_n  = eph{svid,7};
    n = sqrt(mu/(a^3)) + del_n;
    

    
    N_w = floor(-FullBiasNanos/604800e9);
    
    Measurement_Week = N_w;
    Ephem_GPSWeek = eph{svid,23};
    Ephem_Toe = eph{svid,13};
    %3
    tk = (Measurement_Week - Ephem_GPSWeek) *604800+...
        (Measurement_ReceivedSvTimeNanos/1e9 - Ephem_Toe) + add_time;
    %4
    M0 = eph{svid,8};
    Mk = M0 + n*tk;

    
    %5
    e = eph{svid,10};
    tol = 1e-10;
    
     Ek = M2E (Mk , e , tol );
     Ek = wrapTo2Pi(Ek);
%      Mk = Mk
%      test = Ek - e*sin(Ek)

    %6
    sinnuk = sqrt(1 - e^2)*sin(Ek) / (1 - e*cos(Ek));
    cosnuk = (cos(Ek) - e) / (1 - e*cos(Ek));
    nuk = atan2(sinnuk, cosnuk);
    nuk = wrapTo2Pi(nuk);
    %7
    omega = eph{svid,19};
    phik = nuk + omega;
    phik = wrapTo2Pi(phik);
    
    %8
    Cus = eph{svid,11};
    Cuc = eph{svid,9};
    
    phik_temp = phik;
    for j = 1:100
        del_phik = Cus*sin(2*phik_temp) + Cuc* cos(2*phik_temp);

        %9
        uk  = phik + del_phik;
        phik_temp = uk;
    end
    
    %10
    Crs = eph{svid,6};
    
    Crc = eph{svid,18};
    
    del_rk = Crs * sin(2*phik) + Crc * cos(2*phik);
    
    %11
    Cis = eph{svid,16};
    
    Cic = eph{svid,14};
    
    del_ik = Cis * sin(2*phik) + Cic * cos(2*phik);
    
    %12
    Omega_0 = eph{svid,15};
    Omega_dot = eph{svid,20};
    t = (Measurement_Week - Ephem_GPSWeek) *604800+...
        Measurement_ReceivedSvTimeNanos/1e9 + add_time;
    Omega_k = Omega_0 - Omega_dot_E * t + Omega_dot * tk;
    
    %13
    rk  = a *(1 - e*cos(Ek)) + del_rk;
    %14
    i0 = eph{svid,17};
    
    i_dot = eph{svid,21};
    
    ik = i0 + i_dot * tk + del_ik;
    
    %15
    xp = rk * cos(uk);
    
    %16
    
    yp = rk * sin(uk);
    
    %17
    
    Xecef = xp * cos(Omega_k) - yp *cos(ik) *sin(Omega_k);
    Yecef = xp * sin(Omega_k) + yp *cos(ik) *cos(Omega_k);
    Zecef = yp *sin(ik);
    
    
    
    function E = M2E (M , e , tol )
        % M2E solves Kepler 's equation for eccentric anomaly
        %
        % Note : This function uses a Newton - Raphson method to
        % numerically compute
        % the correct value for E
        % %
        % Inputs :
        % M - mean anomaly [ rad ]
        % e - eccentricity of orbit
        % tol - tolerance for Newton - Raphson iterator
        %
        % Outputs :
        % E - eccentric anomaly [ rad ]
        if M == 0 || M == pi
            Return the known solutions ( trivial )
            E = M ;
        else
            E0 = M;
            d = -( E0 - e * sin ( E0 ) - M ) /(1 - e * cos ( E0 ) ) ;
            
            % Loop until the solution converges
            while abs ( d ) > tol
                E1 = E0 + d ;
                d = -( E1 - e * sin ( E1 ) - M ) /(1 - e * cos ( E1 ) ) ;
                E0 = E1 ;
            end
            E = E0 ;
         end
    end
end