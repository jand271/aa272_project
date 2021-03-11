classdef TestPseudoRangeGroup < matlab.unittest.TestCase
    %TESTPSEUDORANGEGROUP Class to test PsuedoRangeGroup
    
    methods(Test)
        function test_against_hw2(testCase)
            load(fullfile(fileparts(mfilename('fullpath')),'hw2_first_gnss_solution_data.mat'));
            x = -2.704071790285519e+06;
            y = -4.263427501694307e+06;
            z = 3.885149158026692e+06;
            b = -2.169197695509361e+03;
            xr_correct = [x;y;z;b];
            
            prg = PsuedoRangeGroupGNSSLog(gnsslogdata, false);
            xr_computed = prg.solve_newton_raphson();
            
            testCase.verifyEqual(xr_computed,xr_correct,'abstol',1e-7);
        end
        
        function test_against_hw3(testCase)
            
            % get correct satellite ECEF positions provided in hw3
            load(fullfile(fileparts(mfilename('fullpath')),'hw3_solution_data.mat'));
            xsats_correct = gnsslogdata.x_ECEFs_c;
            ysats_correct = gnsslogdata.y_ECEFs_c;
            zsats_correct = gnsslogdata.z_ECEFs_c;
            bsats_correct = gnsslogdata.b_ECEFs_c;
            
            % remove the satellite ECEF positions
            gnsslogdata_woXYZB = gnsslogdata(:,1:end-4);
            
            % load ephemeris file
            load(fullfile(fileparts(mfilename('fullpath')),'hw3_ephemeris_data.mat')); 
            
            % regenerate satellite ECEF positions
            gnsslogdata_wXYZB = SatelliteECEFs.append_satellite_positions(ephemeris, gnsslogdata_woXYZB);
            
            % get computed satellite ECEF positions
            xsats_computed = gnsslogdata_wXYZB.X;
            ysats_computed = gnsslogdata_wXYZB.Y;
            zsats_computed = gnsslogdata_wXYZB.Z;
            bsats_computed = gnsslogdata_wXYZB.B;
            
            % assert that computed and correct values are close to each
            % other
            testCase.verifyEqual(xsats_computed,xsats_correct,'abstol',5);
            testCase.verifyEqual(ysats_computed,ysats_correct,'abstol',5);
            testCase.verifyEqual(zsats_computed,zsats_correct,'abstol',5);
            testCase.verifyEqual(bsats_computed,bsats_correct,'abstol',5);
        end
        
        function test_multi_constelation_friday_oval_data(testCase)
            
            return;
            
            % load Madison's Friday Oval Data
            load(fullfile(fileparts(mfilename('fullpath')),'MM_Friday_Oval_data.mat'));
            
            % build tables that seperate the data into constellations
            gnsslogdata_ALL = gnsslogdata(gnsslogdata.ConstellationType ~= 4,:);
            gnsslogdata_GPS = gnsslogdata(gnsslogdata.ConstellationType == 1,:);
            gnsslogdata_GLO = gnsslogdata(gnsslogdata.ConstellationType == 3,:);
            gnsslogdata_BEI = gnsslogdata(gnsslogdata.ConstellationType == 5,:);
            gnsslogdata_GAL = gnsslogdata(gnsslogdata.ConstellationType == 6,:);
            
            % appends satellite ECEF positions to the data in data fields
            % X, Y, Z, B
            gnsslogdata_ALL = SatelliteECEFs.append_satellite_positions(eph, gnsslogdata_ALL);
            gnsslogdata_BEI = SatelliteECEFs.append_satellite_positions(eph, gnsslogdata_BEI);
            gnsslogdata_GPS = SatelliteECEFs.append_satellite_positions(eph, gnsslogdata_GPS);
            gnsslogdata_GLO = SatelliteECEFs.append_satellite_positions(eph, gnsslogdata_GLO);
            gnsslogdata_GAL = SatelliteECEFs.append_satellite_positions(eph, gnsslogdata_GAL);
            
            % Construct PseudoRangeGroups with seperated psuedorange data
            prg_ALL = PsuedoRangeGroupGNSSLog(gnsslogdata_ALL, false);
            prg_GPS = PsuedoRangeGroupGNSSLog(gnsslogdata_GPS, false);
            prg_GLO = PsuedoRangeGroupGNSSLog(gnsslogdata_GLO, false);
            prg_BEI = PsuedoRangeGroupGNSSLog(gnsslogdata_BEI, false);
            prg_GAL = PsuedoRangeGroupGNSSLog(gnsslogdata_GAL, false);
            
            % compute positions with each seperated  psuedorange group
            xr_ALL = prg_ALL.solve_newton_raphson();
            xr_GPS = prg_GPS.solve_newton_raphson();
            xr_GLO = prg_GLO.solve_newton_raphson();
            xr_BEI = prg_BEI.solve_newton_raphson();
            xr_GAL = prg_GAL.solve_newton_raphson();
            
            % assert that the positions are within 20m of each other
            testCase.verifyLessThan(norm(xr_GPS-xr_ALL), 20);
            testCase.verifyLessThan(norm(xr_GLO-xr_ALL), 20);
            testCase.verifyLessThan(norm(xr_BEI-xr_ALL), 20);
            testCase.verifyLessThan(norm(xr_GAL-xr_ALL), 20);
        end  
        function test_against_hw4(testCase)
            load(fullfile(fileparts(mfilename('fullpath')),'hw4_first_gnss_solution_data.mat'));
            x = -2703950.14549877;
            y = -4263087.17552819;
            z = 3885060.28066444;
            b = -1300.23569491921;
            xr_correct = [x;y;z;b];
            
            prg = PsuedoRangeGroupGNSSLog(gnsslogdata, true);
            xr_computed = prg.solve_newton_raphson();
            
            testCase.verifyEqual(xr_computed,xr_correct,'abstol',1e-7);
        end
        function test_multi_constelation_friday_oval_data_GPS_BDS_GLA(testCase)
            
            % load Madison's Friday Oval Data
            load(fullfile(fileparts(mfilename('fullpath')),'MM_Friday_Oval_data.mat'));
             
            % build tables that seperate the data into constellations
            gnsslogdata_ALL = gnsslogdata(...
                bitand(gnsslogdata.ConstellationType ~= 4, gnsslogdata.ConstellationType ~= 3),:);
            gnsslogdata_GPS = gnsslogdata(gnsslogdata.ConstellationType == 1,:);
            gnsslogdata_BEI = gnsslogdata(gnsslogdata.ConstellationType == 5,:);
            gnsslogdata_GAL = gnsslogdata(gnsslogdata.ConstellationType == 6,:);
            
            % appends satellite ECEF positions to the data in data fields
            % X, Y, Z, B
            gnsslogdata_ALL = SatelliteECEFs.append_satellite_positions(ephemeris_data, gnsslogdata_ALL);
            gnsslogdata_BEI = SatelliteECEFs.append_satellite_positions(ephemeris_data, gnsslogdata_BEI);
            gnsslogdata_GPS = SatelliteECEFs.append_satellite_positions(ephemeris_data, gnsslogdata_GPS);
            gnsslogdata_GAL = SatelliteECEFs.append_satellite_positions(ephemeris_data, gnsslogdata_GAL);
            
            % Construct PseudoRangeGroups with seperated psuedorange data
            prg_ALL = PsuedoRangeGroupGNSSLog(gnsslogdata_ALL, false);
            prg_GPS = PsuedoRangeGroupGNSSLog(gnsslogdata_GPS, false);
            prg_BEI = PsuedoRangeGroupGNSSLog(gnsslogdata_BEI, false);
            prg_GAL = PsuedoRangeGroupGNSSLog(gnsslogdata_GAL, false);
            
            % compute positions with each seperated  psuedorange group
            xr_ALL = prg_ALL.solve_newton_raphson();
            xr_GPS = prg_GPS.solve_newton_raphson();
            xr_BEI = prg_BEI.solve_newton_raphson();
            xr_GAL = prg_GAL.solve_newton_raphson();
            
            % assert that the positions are within 20m RMS of each other
            testCase.verifyLessThan(norm(xr_GPS-xr_ALL), 20);
            testCase.verifyLessThan(norm(xr_BEI-xr_ALL), 20);
            testCase.verifyLessThan(norm(xr_GAL-xr_ALL), 20);
        end  
    end
end

