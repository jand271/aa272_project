classdef PsuedoRangeGroupGNSSLog < PsuedoRangeGroup
    %PSUEDORANGEGROUPGNSSLOG Methods pertaining to PsuedoRangeGroups from
    %GNSS Log data
    
    properties(Access = protected)
        gnss_log_data
    end
    
    methods
        function obj = PsuedoRangeGroupGNSSLog(gnss_log_data)
            
            mrhos = PsuedoRangeGroupGNSSLog.data_pseudoranges(gnss_log_data);
            xsats = gnss_log_data.X;
            ysats = gnss_log_data.Y;
            zsats = gnss_log_data.Z;
            bsats = gnss_log_data.B;
            
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
    end
end

