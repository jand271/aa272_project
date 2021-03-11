classdef PsuedoRangeGroupGNSSLog < PsuedoRangeGroup
    %PSUEDORANGEGROUPGNSSLOG Methods pertaining to PsuedoRangeGroups from
    %GNSS Log data
    
    properties(Access = protected)
        gnss_log_data
    end
    
    methods
        function obj = PsuedoRangeGroupGNSSLog(gnss_log_data, correct_ionosphere)
            
            mrhos = PsuedoRangeGroupGNSSLog.data_pseudoranges(gnss_log_data);
            xsats = gnss_log_data.X;
            ysats = gnss_log_data.Y;
            zsats = gnss_log_data.Z;
            bsats = gnss_log_data.B;
            
            if correct_ionosphere
                % correct ionosphere
                svids = gnss_log_data.Svid;
                bands = PsuedoRangeGroupGNSSLog.determine_bands(gnss_log_data);
                [mrhos, xsats, ysats, zsats, bsats] = PsuedoRangeGroupGNSSLog. correct_ionosphere(...
                    svids, bands, mrhos, xsats, ysats, zsats, bsats);
            end
            
            obj@PsuedoRangeGroup(mrhos,xsats,ysats,zsats,bsats);
            obj.gnss_log_data = gnss_log_data;
        end
    end
    
    methods(Static, Access=protected)
        function rho_m = data_pseudoranges(gnss_log_data)
                        
            NANO_PER_WEEK = 604800e9;
            NANO_PER_DAY = 86400e9;
            
            t_Tx = gnss_log_data.ReceivedSvTimeNanos;
            t_Rx_GNSS = gnss_log_data.TimeNanos + gnss_log_data.TimeOffsetNanos - ...
                (gnss_log_data.FullBiasNanos(1)+gnss_log_data.BiasNanos(1));
            
            n_measurements = size(gnss_log_data, 1);
            t_Rx = zeros(n_measurements, 1);
            for i = 1:n_measurements
                switch gnss_log_data.ConstellationType(i)
                    case {1,6}
                        weekNumberNanos = floor(-gnss_log_data.FullBiasNanos(1)/NANO_PER_WEEK) * NANO_PER_WEEK;
                        t_Rx(i) = t_Rx_GNSS(i) - weekNumberNanos;
                    case 5
                        weekNumberNanos = floor(-gnss_log_data.FullBiasNanos(1)/NANO_PER_WEEK) * NANO_PER_WEEK;
                        t_Rx(i) = t_Rx_GNSS(i) - weekNumberNanos - 14e9;
                    case 3
                        DayNumberNanos = floor(-gnss_log_data.FullBiasNanos(1)/NANO_PER_DAY) * NANO_PER_DAY;
                        t_Rx(i) = t_Rx_GNSS(i) - DayNumberNanos + 3 * 3600e9 - (37-18)*1e9; % 18 needs to not be static!!
                    otherwise
                        error('Constellation Not implemmented');
                end
            end
            
            rho_m = (t_Rx - t_Tx) / 1e9 * 299792458;
        end   
        function bands = determine_bands(gnss_log_data)
            n_measurements = size(gnss_log_data,1);
            bands = zeros(n_measurements, 1);
            fs = gnss_log_data.CarrierFrequencyHz;
            
            for i = 1:n_measurements
                if abs(fs(i) - 1575.42e6) < 1e6
                    bands(i) = 1;
                elseif abs(fs(i) - 1176.45e6) < 1e6
                    bands(i) = 5;
                else
                    error('Bad Carrier Frequency Received');
                end   
            end
        end
        function [mrhos_out, xsats_out, ysats_out, zsats_out, bsats_out] = correct_ionosphere(...
                svids, bands, mrhos, xsats, ysats, zsats, bsats)
            
            % very inefficient wrt data structure complexity
            
            svs_unique = unique(svids);
            n_svs = length(svs_unique);
            
            mrhos_out = zeros(n_svs, 1);
            xsats_out = zeros(n_svs, 1);
            ysats_out = zeros(n_svs, 1);
            zsats_out = zeros(n_svs, 1);
            bsats_out = zeros(n_svs, 1);
            
            for i = 1:n_svs
                n_svid = sum(svids == svs_unique(i));
                
                switch n_svid
                    case 1
                        index = find(svids == svs_unique(i));
                        mrhos_out(i) = mrhos(index);
                        xsats_out(i) = xsats(index);
                        ysats_out(i) = ysats(index);
                        zsats_out(i) = zsats(index);
                        bsats_out(i) = bsats(index);
                    case 2
                        indices = svids == svs_unique(i);
                        L1_index = 1 == bands .* indices;
                        L5_index = 5 == bands .* indices;
                        
                        f_L1 = 1575.42e6;
                        f_L5 = 1176.45e6;
                        L1_weight = f_L1^2 / (f_L1^2 - f_L5^2);
                        L5_weight = f_L5^2 / (f_L1^2 - f_L5^2);
                        
                        mrhos_out(i) = L1_weight * mrhos(L1_index) - L5_weight * mrhos(L5_index);
                        xsats_out(i) = mean(xsats(indices));
                        ysats_out(i) = mean(ysats(indices));
                        zsats_out(i) = mean(zsats(indices));
                        bsats_out(i) = mean(bsats(indices));
                    otherwise
                        error('Received bad number of band measurements')
                end
            end
        end
    end
end

