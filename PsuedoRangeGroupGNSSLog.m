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
            N_w = floor(-gnss_log_data.FullBiasNanos/604800e9);
            trx_hw = gnss_log_data.TimeNanos + gnss_log_data.TimeOffsetNanos;
            b_hw = gnss_log_data.FullBiasNanos + gnss_log_data.BiasNanos;
            trx_gps = trx_hw - b_hw;
            trx_w = trx_gps - N_w * 604800e9;
            rho_ns = trx_w - gnss_log_data.ReceivedSvTimeNanos;
            rho_m = rho_ns * 299792458/1e9;   
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

