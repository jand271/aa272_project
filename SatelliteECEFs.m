classdef SatelliteECEFs
    %SATELLITEECEFS Class pertaining to Satellite ECEF Positions
    
    properties(Constant, Access=private)
        mu = 3.986004418e14;
        c = 299792458;
        omega_dot_e = 7.2921151467e-5;
    end
    
    methods(Static, Access=public)
        function gnsslogdata = append_satellite_positions(ephemeris, gnsslogdata)
            
            ephemeris = SatelliteECEFs.fix_ephemeris_headers(ephemeris);
            
            n = size(gnsslogdata,1);
            Xs = zeros(n,1);
            Ys = zeros(n,1);
            Zs = zeros(n,1);
            Bs = zeros(n,1);
            
            for i = 1:n
                measurement_week = floor(-gnsslogdata.FullBiasNanos(i)/604800e9);
                measurement_time = gnsslogdata.ReceivedSvTimeNanos(i)/1e9;
                svid = gnsslogdata.Svid(i);
                switch gnsslogdata.ConstellationType(i)
                    case {1, 5, 6}
                        switch gnsslogdata.ConstellationType(i)
                            case 1
                                constellation_ephemeris = ephemeris.gps;
                            case 5
                                constellation_ephemeris = ephemeris.bds;
                            case 6
                                constellation_ephemeris = ephemeris.gal;
                        end
                        
                        toes = constellation_ephemeris.Toe;
                        toes(constellation_ephemeris.PRN ~= gnsslogdata.Svid(i)) = nan;
                        [~, time_index] = min(abs(toes - measurement_time));
                        
                        if isempty(time_index)
                            % when the hw-based unit test data doesn't
                            % match the RINEX format, fix to pass unit
                            % tests
                            time_index = gnsslogdata.Svid(i);
                        end
                        
                        [X,Y,Z] = SatelliteECEFs.eph_to_ecef(...
                            constellation_ephemeris.sqrtA(time_index), ...
                            constellation_ephemeris.e(time_index), ...
                            constellation_ephemeris.i0(time_index), ...
                            constellation_ephemeris.omega(time_index), ...
                            constellation_ephemeris.Delta_n(time_index), ...
                            constellation_ephemeris.M0(time_index), ...
                            constellation_ephemeris.GPS_week_num(time_index), ...
                            constellation_ephemeris.Cuc(time_index), ...
                            constellation_ephemeris.Cus(time_index), ...
                            constellation_ephemeris.Cic(time_index), ...
                            constellation_ephemeris.Cis(time_index), ...
                            constellation_ephemeris.Crc(time_index), ...
                            constellation_ephemeris.Crs(time_index), ...
                            constellation_ephemeris.OMEGA(time_index), ...
                            constellation_ephemeris.OMEGA_DOT(time_index), ...
                            constellation_ephemeris.I_DOT(time_index), ...
                            constellation_ephemeris.Toe(time_index), ...
                            measurement_week, ...
                            measurement_time);
                        B = constellation_ephemeris.clock_bias(time_index) * SatelliteECEFs.c;
                    case 3
                        X = nan;
                        Y = nan;
                        Z = nan;
                        B = nan;
                    otherwise
                        error('Not Implemented Constellation Time');
                end
                
                Xs(i) = X;
                Ys(i) = Y;
                Zs(i) = Z;
                Bs(i) = B;
            end
            gnsslogdata.X = Xs;
            gnsslogdata.Y = Ys;
            gnsslogdata.Z = Zs;
            gnsslogdata.B = Bs;
        end
        function [X,Y,Z] = hw4_query(ephemeris_struct, svid, query_gps_time)
            
            data_index = find(ephemeris_struct.TOE < query_gps_time & ephemeris_struct.prm == svid, 1, 'last');
            
            [X,Y,Z] = SatelliteECEFs.eph_to_ecef(...
                ephemeris_struct.sqrta(data_index), ...
                ephemeris_struct.e(data_index), ...
                ephemeris_struct.i0(data_index), ...
                ephemeris_struct.omega(data_index), ...
                ephemeris_struct.dn(data_index), ...
                ephemeris_struct.M0(data_index), ...
                ephemeris_struct.gpsweek(data_index), ...
                ephemeris_struct.Cuc(data_index), ...
                ephemeris_struct.Cus(data_index), ...
                ephemeris_struct.Cic(data_index), ...
                ephemeris_struct.Cis(data_index), ...
                ephemeris_struct.Crc(data_index), ...
                ephemeris_struct.Crs(data_index), ...
                ephemeris_struct.Omega0(data_index), ...
                ephemeris_struct.Omega_dot(data_index), ...
                ephemeris_struct.idot(data_index), ...
                ephemeris_struct.TOE(data_index), ...
                ephemeris_struct.gpsweek(data_index), ...
                query_gps_time);
        end
    end
    
    methods(Static, Access=private)
        function [X,Y,Z] = eph_to_ecef(...
                sqrta, ...
                e, ...
                i, ...
                omega, ...
                dn, ...
                M0, ...
                epoch_week, ...
                cuc, ...
                cus, ...
                cic, ...
                cis, ...
                crc, ...
                crs, ...
                Omega_0, ...
                Omega_dot, ...
                idot, ...
                time_of_ephemeris, ...
                measurement_week, ...
                time_measurement)
            
            %1
            a = sqrta^2;
            
            %2
            n = sqrt(SatelliteECEFs.mu/a^3) + dn;
            
            %3
            tk = (measurement_week - epoch_week) * 604800 + time_measurement - time_of_ephemeris;
            
            %4
            Mk = M0 + n*tk;
            
            %5
            Ek = SatelliteECEFs.kepler_equation(Mk, e, 1e-10);
            
            %6
            nuk = atan2(sqrt(1 - e^2) * sin(Ek), cos(Ek) - e);
            
            %7
            phik = nuk + omega;
            
            %8
            del_phik = cus*sin(2*phik) + cuc* cos(2*phik);
            
            %9
            uk  = phik + del_phik;
            
            %10
            del_rk = crs * sin(2*phik) + crc * cos(2*phik);
            
            %11
            del_ik = cis * sin(2*phik) + cic * cos(2*phik);
            
            %12
            t = (measurement_week - epoch_week) * 604800 + time_measurement;
            Omega_k = Omega_0 - SatelliteECEFs.omega_dot_e * t + Omega_dot * tk;
            
            %13
            rk  = a *(1 - e*cos(Ek)) + del_rk;
            
            %14
            ik = i + idot * tk + del_ik;
            
            %15
            xp = rk * cos(uk);
            
            %16
            yp = rk * sin(uk);
            
            %17
            X = xp * cos(Omega_k) - yp *cos(ik) *sin(Omega_k);
            Y = xp * sin(Omega_k) + yp *cos(ik) *cos(Omega_k);
            Z = yp *sin(ik);
        end
        
        function E = kepler_equation(M, e, tol)
            if M == 0 || M == pi
                E = M;
                return;
            else
                E0 = M;
                d = -( E0 - e * sin ( E0 ) - M ) /(1 - e * cos ( E0 ) );
                
                while abs ( d ) > tol
                    E1 = E0 + d;
                    d = -(E1 - e*sin(E1) - M)/(1 - e*cos(E1));
                    E0 = E1;
                end
                E = E0 ;
            end
        end
        function ephemeris = fix_ephemeris_headers(ephemeris)
            
            if istable(ephemeris.gps)
                ephemeris.gps = table2struct(ephemeris.gps,'ToScalar',true);
            end
            if ~isfield(ephemeris.gps,'e')
                ephemeris.gps.e = ephemeris.gps.Eccentricity;
            end
            if ~isfield(ephemeris.gps,'i0')
                ephemeris.gps.i0 = ephemeris.gps.Io;
            end
            if ~isfield(ephemeris.gps,'Delta_n')
                ephemeris.gps.Delta_n = ephemeris.gps.DeltaN;
            end
            if ~isfield(ephemeris.gps,'I_DOT')
                ephemeris.gps.I_DOT = ephemeris.gps.IDOT;
            end
            if ~isfield(ephemeris.gps,'GPS_week_num')
                ephemeris.gps.GPS_week_num = ephemeris.gps.GPSWeek;
            end
            if ~isfield(ephemeris.gps,'OMEGA') % could be wrong
                ephemeris.gps.OMEGA = ephemeris.gps.Omega0;
            end
            if ~isfield(ephemeris.gps,'OMEGA_DOT')
                ephemeris.gps.OMEGA_DOT = ephemeris.gps.OmegaDot;
            end
            if ~isfield(ephemeris.gps,'clock_bias')
                ephemeris.gps.clock_bias = ephemeris.gps.SVclockBias;
            end
            if ~isfield(ephemeris.gps,'PRN')
                ephemeris.gps.PRN = ephemeris.gps.Svid;
            end
        end
    end
end

