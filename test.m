clc
svid = 3;
gnss = gnss_log(1);
[Xecef,Yecef,Zecef,bias] = sat_pos_ecef_sv(eph,gnss);
function [Xecef,Yecef,Zecef,bias] = sat_pos_ecef_sv(eph,gnss)
 % eph = table2cell(eph_t);
 svid = gnss.Svid;
 eph_id = find(eph.PRN==svid);
 eph_id = eph_id(end);
  bias = eph.clock_bias(eph_id)*299792458;
  sqrtA = eph.sqrtA(eph_id);
  %1
  a = sqrtA^2;
  mu = 3.986005e14;
  Omega_dot_E = 7.2921151467e-5;
  %2
  del_n = eph.Delta_n(eph_id);
  n = sqrt(mu/(a^3)) + del_n;
  FullBiasNanos = gnss.FullBiasNanos;
  N_w = floor(-FullBiasNanos/604800e9);
  Measurement_Week = N_w;
  Ephem_GPSWeek = eph.GPS_weekday(eph_id);
  Ephem_Toe = eph.Toe(eph_id);
  %3
  tk = (Measurement_Week - Ephem_GPSWeek) *604800+...
    (Measurement_ReceivedSvTimeNanos/1e9 - Ephem_Toe) + add_time;
  %4
  M0 = eph.M0(eph_id);
  Mk = M0 + n*tk;
  %5
  e = eph.e(eph_id);
  tol = 1e-10;
   Ek = M2E (Mk , e , tol );
   Ek = wrapTo2Pi(Ek);
%   Mk = Mk
%   test = Ek - e*sin(Ek)
  %6
  sinnuk = sqrt(1 - e^2)*sin(Ek) / (1 - e*cos(Ek));
  cosnuk = (cos(Ek) - e) / (1 - e*cos(Ek));
  nuk = atan2(sinnuk, cosnuk);
  nuk = wrapTo2Pi(nuk);
  %7
  omega = eph.omega(eph_id);
  phik = nuk + omega;
  phik = wrapTo2Pi(phik);
  %8
  Cus = eph.Cus(eph_id);
  Cuc = eph.Cus(eph_id);
  phik_temp = phik;
  for j = 1:100
    del_phik = Cus*sin(2*phik_temp) + Cuc* cos(2*phik_temp);
    %9
    uk = phik + del_phik;
    phik_temp = uk;
  end
  %10
  Crs = eph.Crs(eph_id);
  Crc = eph.Crc(eph_id);
  del_rk = Crs * sin(2*phik) + Crc * cos(2*phik);
  %11
  Cis = eph.Cis(eph_id);
  Cic = eph.Cic(eph_id);
  del_ik = Cis * sin(2*phik) + Cic * cos(2*phik);
  %12
  Omega_0 = eph.OMEGA(eph_id);
  Omega_dot = eph.OMEGA_DOT(eph_id);
  t = (Measurement_Week - Ephem_GPSWeek) *604800+...
    Measurement_ReceivedSvTimeNanos/1e9 + add_time;
  Omega_k = Omega_0 - Omega_dot_E * t + Omega_dot * tk;
  %13
  rk = a *(1 - e*cos(Ek)) + del_rk;
  %14
  i0 = eph.i0(eph_id);
  i_dot = eph.I_DOT(eph_id);
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
    % M2E solves Kepler ‘s equation for eccentric anomaly
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
