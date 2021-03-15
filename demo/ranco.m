%% Create RANCO Error Plots

clear all;
close all;

%% Import data

this_file_path = fileparts(mfilename('fullpath'));

ephemeris = struct();
ephemeris.gps = loadRINEXNavigation('G',fullfile(this_file_path,'../data/oval30_1'),'BRDC00WRD_R_20210710000_01D_GN.rnx');
ephemeris.bds = loadRINEXNavigation('C',fullfile(this_file_path,'../data/oval30_1'),'BRDC00WRD_R_20210710000_01D_CN.rnx');
% ephemeris.gps = loadRINEXNavigation('G',fullfile(this_file_path,'../data/gsb30_1'),'BRDC00WRD_R_20210710000_01D_GN.rnx');
% ephemeris.bds = loadRINEXNavigation('C',fullfile(this_file_path,'../data/gsb30_1'),'BRDC00WRD_R_20210710000_01D_CN.rnx');

%gnsslogdata = readtable(fullfile(this_file_path,'../data/gsb30_1/raw.csv'));

gnsslogdata = readtable(fullfile(this_file_path,'../data/oval30_1/raw.csv'));


gnsslogdata_gps_bds = gnsslogdata(bitor(gnsslogdata.ConstellationType == 1, gnsslogdata.ConstellationType == 5),:);
gnsslogdata_gps_bds = SatelliteECEFs.append_satellite_positions(ephemeris, gnsslogdata_gps_bds);

gnsslogdata_gps_bds = gnsslogdata_gps_bds(~isnan(gnsslogdata_gps_bds.X),:);

gnsslogdata_gps_bds = gnsslogdata_gps_bds(gnsslogdata_gps_bds.CarrierFrequencyHz > 1.2e9, :);


%% 

times = unique(gnsslogdata_gps_bds.TimeNanos);
n_times = size(times,1);

times = times(55:n_times);
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
plot(solution_ranco(1,:),solution_ranco(2,:),'c*')
PlottingUtils.plot_percentile_ellipse(solution_ranco(1:2,:), 0.05, 0.90,'b');
% plot(satellite_exclusion_solution(1,:),satellite_exclusion_solution(2,:),'r*')
% PlottingUtils.plot_percentile_ellipse(satellite_exclusion_solution(1:2,:), 0.05, 0.90,'r');
plot(solution_raw(1,:),solution_raw(2,:),'g*');
PlottingUtils.plot_percentile_ellipse(solution_raw, 0.05, 0.90,'k');
legend('ranco','ranco 90%','raw','raw 90%','location','best')
xlabel('X (m)')
ylabel('Y (m)')
title('Compare ECEF solution of RACO against raw data')
hold off



figure
hold on
% plot(solution_ranco(1,:),solution_ranco(2,:),'b*')
% PlottingUtils.plot_percentile_ellipse(solution_ranco(1:2,:), 0.05, 0.90,'b');
plot(satellite_exclusion_solution(1,:),satellite_exclusion_solution(2,:),'c*')
PlottingUtils.plot_percentile_ellipse(satellite_exclusion_solution(1:2,:), 0.05, 0.90,'b');
plot(solution_raw(1,:),solution_raw(2,:),'g*');
PlottingUtils.plot_percentile_ellipse(solution_raw, 0.05, 0.90,'k');
legend('satellite exclusion','satellite exclusion 90%','raw','raw 90%','location','best')
hold off
xlabel('X (m)')
ylabel('Y (m)')
title('Compare ECEF solution of satellite exclusion method against raw data')

% legend('ranco','ranco 90%','satellite exclusion','satellite exclusion 90%','raw','raw 90%','location','best')
%legend('ranco','satellite exclusion','raw','location','best')
solution_ranco1 = solution_ranco(1,:);
solution_ranco1 = solution_ranco1(~isnan(solution_ranco1));
solution_ranco2 = solution_ranco(2,:);
solution_ranco2 = solution_ranco2(~isnan(solution_ranco2));

ranco_mean_solution = [mean(solution_ranco1),mean(solution_ranco2)]

solution_raw1 = solution_raw(1,:);
solution_raw1 = solution_raw1(~isnan(solution_raw1));
solution_raw2 = solution_raw(2,:);
solution_raw2 = solution_raw2(~isnan(solution_raw2));

raw_mean_solution = [mean(solution_raw1),mean(solution_raw2)]

% 37.428153513491736, -122.16179918380078 23.4696  gsb

%37.429713, -122.169527 21.336ovel 

Truth = lla2ecef([37.429713 -122.169527 21.336]);
Truth = Truth(1:2);

ranco_diff = norm(ranco_mean_solution - Truth)
raw_diff = norm(raw_mean_solution - Truth)

% xlabel('X (m)')
% ylabel('Y (m)')
% title('Compare ECEF solution of RANCO against raw data')