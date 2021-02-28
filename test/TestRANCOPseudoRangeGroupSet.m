classdef TestRANCOPseudoRangeGroupSet < matlab.unittest.TestCase
    %TESTRANCOPSEUDORANGEGROUP Class to test RANCOPsuedoRangeGroupSet
    
    methods(Test)
        function test_ranco_group_size(testCase)
            load(fullfile(fileparts(mfilename('fullpath')),'hw2_first_gnss_solution_data.mat'));
            
            prg = PsuedoRangeGroupGNSSLog(gnsslogdata);
            
            rprgs = RANCOPsuedoRangeGroupSet(prg, 4);
            
            number_of_ranco_sets = nchoosek(prg.n_measurements, 4);
            
            testCase.verifyEqual(length(rprgs.psuedoRangeGroupSet), number_of_ranco_sets);
        end
    end
end
