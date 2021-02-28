%% Create RANCO Error Plots

clear all;
close all;

%% Load some nice data from the AA272 HW

% load mat file from repo
load(fullfile(fileparts(mfilename('fullpath')),'../test','hw2_first_gnss_solution_data.mat'));

% create PsuedoRangeGroup object
prg = PsuedoRangeGroupGNSSLog(gnsslogdata);

%% Compute RANCO Errors

% generate the permuation of all combinations of 4
rprgs = RANCOPsuedoRangeGroupSet(prg, 4);

% compute the errors of all satellites on each subset of 4
sat_errors = rprgs.satellite_errors();

%% Plot

figure;
plot(sat_errors, 'bo');
xlabel('Subset Combitation Number');
ylabel('Error Equation 6 from RANCO Paper');
title('RANCO Errors For Each Satellite Against Each Subset');
