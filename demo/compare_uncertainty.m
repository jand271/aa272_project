
%% Import data

this_file_path = fileparts(mfilename('fullpath'));

ephemeris = struct();
ephemeris.gps = loadRINEXNavigation('G',fullfile(this_file_path,'../data/oval30_1'),'BRDC00WRD_R_20210710000_01D_GN.rnx');
ephemeris.bds = loadRINEXNavigation('C',fullfile(this_file_path,'../data/oval30_1'),'BRDC00WRD_R_20210710000_01D_CN.rnx');

gnsslogdata = readtable(fullfile(this_file_path,'../data/oval30_1/raw.csv'));

gnsslogdata_gps_bds = gnsslogdata(bitor(gnsslogdata.ConstellationType == 1, gnsslogdata.ConstellationType == 5),:);
gnsslogdata_gps_bds = SatelliteECEFs.append_satellite_positions(ephemeris, gnsslogdata_gps_bds);

gnsslogdata_gps_bds = gnsslogdata_gps_bds(~isnan(gnsslogdata_gps_bds.X),:);


%% 

times = unique(gnsslogdata_gps_bds.TimeNanos);
n_times = size(times,1);

times = times(50:n_times);
n_times = size(times,1);

xrs = zeros(4, n_times);

for t = 1:n_times
   gnsslogdata_instance = gnsslogdata_gps_bds(gnsslogdata_gps_bds.TimeNanos == times(t),:);
   xrs(:,t) = PsuedoRangeGroupGNSSLog(gnsslogdata_instance, false).solve_newton_raphson();
end

static_cov = cov(xrs(1:2,:)');
static_mean = mean(xrs(1:2,:)')';

%%

figure;
hold on;
plot(xrs(1,:),xrs(2,:),'bo');
e = error_ellipse(static_cov, static_mean,'conf',0.95);
e.LineWidth = 4;

%% BootStrap

number_of_samples = 1000;
number_per_group = 13;
 
t = floor(n_times/2);
gnsslogdata_instance = gnsslogdata_gps_bds(gnsslogdata_gps_bds.TimeNanos == times(t),:);
prg = PsuedoRangeGroupGNSSLog(gnsslogdata_instance, false);
brpg = BootstrapPsuedoRangeGroupSet(prg, number_of_samples, number_per_group);
xrs_bs = brpg.solve_each_newton_raphson();

plot(xrs_bs(1,:),xrs_bs(2,:),'ro');



mean_bs = brpg.mean_newton_raphson_solution();
cov_bs = brpg.covariance_newton_raphson_solution();

mean_bs = mean_bs(1:2);
cov_bs = cov_bs(1:2,1:2);

e = error_ellipse(cov_bs, mean_bs,'conf',0.95);
e.LineWidth = 4;

[bs_lowers, bs_uppers] = brpg.bs_ranges(0.975,0.025);
bs_means = (bs_uppers + bs_lowers)/2;
bs_radii = (bs_uppers - bs_lowers)/2;

PlottingUtils.plot_ellipse_2d(mean_bs(1:2), bs_radii(1:2))








