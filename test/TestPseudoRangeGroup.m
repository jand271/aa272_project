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
    end
end

