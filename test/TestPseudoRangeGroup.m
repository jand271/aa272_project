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
            
            prg = PsuedoRangeGroupGNSSLog(gnsslogdata);
            xr_computed = prg.solve_newton_raphson();
            
            testCase.verifyEqual(xr_computed,xr_correct,'abstol',1e-7);
        end
        
        function test_against_hw3(testCase)
            
            return; % REMOVE THIS STATEMENT WHEN append_satellite_positions IS READY
            
            % get correct satellite ECEF positions provided in hw2
            load(fullfile(fileparts(mfilename('fullpath')),'hw2_first_gnss_solution_data.mat'));
            xsats_correct = gnsslogdata.X;
            ysats_correct = gnsslogdata.Y;
            zsats_correct = gnsslogdata.Z;
            bsats_correct = gnsslogdata.B;
            
            % remove the satellite ECEF positions
            gnsslogdata_woXYZB = gnsslogdata(:,1:end-4);
            
            % load ephemeris file
            load(fullfile(fileparts(mfilename('fullpath')),'hw3_ephemeris_data.mat')); 
            
            % regenerate satellite ECEF positions
            gnsslogdata_wXYZB = append_satellite_positions(eph, gnsslogdata_woXYZB);
            
            % get computed satellite ECEF positions
            xsats_computed = gnsslogdata_wXYZB.X;
            ysats_computed = gnsslogdata_wXYZB.Y;
            zsats_computed = gnsslogdata_wXYZB.Z;
            bsats_computed = gnsslogdata_wXYZB.B;
            
            % assert that computed and correct values are close to each
            % other
            testCase.verifyEqual(xsats_computed,xsats_correct,'abstol',1e-7);
            testCase.verifyEqual(ysats_computed,ysats_correct,'abstol',1e-7);
            testCase.verifyEqual(zsats_computed,zsats_correct,'abstol',1e-7);
            testCase.verifyEqual(bsats_computed,bsats_correct,'abstol',1e-7);
        end
        
        function test_multi_constelation_10_min_data(testCase)
            
            return; % REMOVE THIS STATEMENT WHEN append_satellite_positions IS READY
            
            % load Madison's 10-min Data set
            load(fullfile(fileparts(mfilename('fullpath')),'MM_10_min_data.mat'));
            
            % build tables that seperate the data into constellations
            gnsslogdata_ALL = gnsslogdata;
            gnsslogdata_GPS = gnsslogdata(gnsslogdata.ConstellationType == 1,:);
            gnsslogdata_GLO = gnsslogdata(gnsslogdata.ConstellationType == 3,:);
            gnsslogdata_BEI = gnsslogdata(gnsslogdata.ConstellationType == 5,:);
            gnsslogdata_GAL = gnsslogdata(gnsslogdata.ConstellationType == 6,:);
            
            % appends satellite ECEF positions to the data in data fields
            % X, Y, Z, B
            gnsslogdata_ALL = append_satellite_positions(eph, gnsslogdata_ALL);
            gnsslogdata_GPS = append_satellite_positions(eph, gnsslogdata_GPS);
            gnsslogdata_GLO = append_satellite_positions(eph, gnsslogdata_GLO);
            gnsslogdata_BEI = append_satellite_positions(eph, gnsslogdata_BEI);
            gnsslogdata_GAL = append_satellite_positions(eph, gnsslogdata_GAL);
           
            % Construct PseudoRangeGroups with seperated psuedorange data
            prg_ALL = PsuedoRangeGroupGNSSLog(gnsslogdata_ALL);
            prg_GPS = PsuedoRangeGroupGNSSLog(gnsslogdata_GPS);
            prg_GLO = PsuedoRangeGroupGNSSLog(gnsslogdata_GLO);
            prg_BEI = PsuedoRangeGroupGNSSLog(gnsslogdata_BEI);
            prg_GAL = PsuedoRangeGroupGNSSLog(gnsslogdata_GAL);
            
            % compute positions with each seperated  psuedorange group
            xr_ALL = prg_ALL.solve_newton_raphson();
            xr_GPS = prg_GPS.solve_newton_raphson();
            xr_GLO = prg_GLO.solve_newton_raphson();
            xr_BEI = prg_BEI.solve_newton_raphson();
            xr_GAL = prg_GAL.solve_newton_raphson();
            
            % assert that the positions are within 1m RMS of each other
            testCase.verifyLessThan(norm(xr_GPS-xr_ALL), 1);
            testCase.verifyLessThan(norm(xr_GLO-xr_ALL), 1);
            testCase.verifyLessThan(norm(xr_BEI-xr_ALL), 1);
            testCase.verifyLessThan(norm(xr_GAL-xr_ALL), 1);
        end  
    end
end

