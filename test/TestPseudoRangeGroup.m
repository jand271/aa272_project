classdef TestPseudoRangeGroup < matlab.unittest.TestCase
    %TESTPSEUDORANGEGROUP Class to test PsuedoRangeGroup
    
    methods(Test)
        function test_against_hw2(testCase)
            load(fullfile(fileparts(mfilename('fullpath')),'hw2_first_gnss_solution_data.mat'));
            x = -2.703952629951560e+06;
            y = -4.263104282040041e+06;
            z = 3.885071813390337e+06;
            b = -1.281768152170775e+03;
            xr_correct = [x;y;z;b];
            
            prg = PsuedoRangeGroupGNSSLog(gnsslogdata);
            xr_computed = prg.solve_newton_raphson();
            
            testCase.verifyEqual(xr_computed,xr_correct,'abstol',1e-9);
        end
    end
end

