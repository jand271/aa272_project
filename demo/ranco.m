%% Create RANCO Error Plots

clear all;
close all;

%% Import data

this_file_path = fileparts(mfilename('fullpath'));

ephemeris = struct();
ephemeris.gps = loadRINEXNavigation('G',fullfile(this_file_path,'../data/oval30_1'),'BRDC00WRD_R_20210710000_01D_GN.rnx');
ephemeris.bds = loadRINEXNavigation('C',fullfile(this_file_path,'../data/oval30_1'),'BRDC00WRD_R_20210710000_01D_CN.rnx');

gnsslogdata = readtable(fullfile(this_file_path,'../data/oval30_1/raw.csv'));

gnsslogdata_gps_bds = gnsslogdata(bitor(gnsslogdata.ConstellationType == 1, gnsslogdata.ConstellationType == 5),:);
gnsslogdata_gps_bds = SatelliteECEFs.append_satellite_positions(ephemeris, gnsslogdata_gps_bds);

gnsslogdata_gps_bds = gnsslogdata_gps_bds(~isnan(gnsslogdata_gps_bds.X),:);

gnsslogdata_gps_bds = gnsslogdata_gps_bds(gnsslogdata_gps_bds.CarrierFrequencyHz > 1.2e9, :);


%% 

times = unique(gnsslogdata_gps_bds.TimeNanos);
n_times = size(times,1);

times = times(50:n_times);
n_times = size(times,1);

xrs = zeros(4, n_times);


satellite_exclusion_solution = zeros(4, n_times);
solution_ranco = zeros(4, n_times);
solution_raw = zeros(4, n_times);
parfor t = 1:n_times
   gnsslogdata_instance = gnsslogdata_gps_bds(gnsslogdata_gps_bds.TimeNanos == times(t),:);
   solution_raw(:,t) = PsuedoRangeGroupGNSSLog(gnsslogdata_instance,false).solve_newton_raphson();
   [solution_ranco(:,t), satellite_exclusion_solution(:,t)] = ranco_satellite(gnsslogdata_instance);
end

%%
figure
hold on
plot3(solution_ranco(1,:),solution_ranco(2,:),solution_ranco(3,:),'*')
plot3(satellite_exclusion_solution(1,:),satellite_exclusion_solution(2,:),satellite_exclusion_solution(3,:),'*')
plot3(solution_raw(1,:),solution_raw(2,:),solution_raw(3,:),'*');
legend('ranco','satellite exclusion','raw')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('comparision solution for ranco vs satellite separation')