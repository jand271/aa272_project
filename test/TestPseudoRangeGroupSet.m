classdef TestPseudoRangeGroupSet < matlab.unittest.TestCase
    %TESTPSEUDORANGEGROUP Class to test PsuedoRangeGroupSet
    
    methods(Test)
        function test_repeating_group(testCase)
            load(fullfile(fileparts(mfilename('fullpath')),'hw2_first_gnss_solution_data.mat'));
            
            prg = PsuedoRangeGroupGNSSLog(gnsslogdata, false);
            xr_correct = prg.solve_newton_raphson();
            
            prgs = PsuedoRangeGroupSet([prg,prg,prg,prg,prg,prg]);
            xr_computed = prgs.mean_newton_raphson_solution();
            xr_cov_computed = prgs.covariance_newton_raphson_solution();
            
            testCase.verifyEqual(xr_computed,xr_correct,'abstol',1e-9);
            testCase.verifyEqual(xr_cov_computed,zeros(4,4),'abstol',1e-9);
        end
        
        function test_bootstrap_group(testCase)
            rng(0);
            load(fullfile(fileparts(mfilename('fullpath')),'hw2_first_gnss_solution_data.mat'));
            
            prg = PsuedoRangeGroupGNSSLog(gnsslogdata, false);
            xr_correct = prg.solve_newton_raphson();
            
            prgs = BootstrapPsuedoRangeGroupSet(prg, 100, 5);
        end
    end
end

