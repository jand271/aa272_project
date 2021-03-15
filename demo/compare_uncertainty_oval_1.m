
close all;
clear all;

%% Import data

this_file_path = fileparts(mfilename('fullpath'));

ephemeris = struct();
ephemeris.gps = loadRINEXNavigation('G',fullfile(this_file_path,'../data/oval30_1'),'BRDC00WRD_R_20210710000_01D_GN.rnx');
ephemeris.bds = loadRINEXNavigation('C',fullfile(this_file_path,'../data/oval30_1'),'BRDC00WRD_R_20210710000_01D_CN.rnx');
gnsslogdata = readtable(fullfile(this_file_path,'../data/oval30_1/raw.csv'));

gps_bds = bitor(gnsslogdata.ConstellationType == 5, gnsslogdata.ConstellationType == 1);
gnsslogdata = gnsslogdata(gps_bds,:);
gnsslogdata = SatelliteECEFs.append_satellite_positions(ephemeris, gnsslogdata);

gnsslogdata_nonan = gnsslogdata(~isnan(gnsslogdata.X),:);


%% 

times = unique(gnsslogdata_nonan.TimeNanos);
n_times = size(times,1);

times = times(50:n_times);
n_times = size(times,1);

xrs = zeros(4, n_times);

for t = 1:n_times
   gnsslogdata_instance = gnsslogdata_nonan(gnsslogdata_nonan.TimeNanos == times(t),:);
   try
       xrs(:,t) = PsuedoRangeGroupGNSSLog(gnsslogdata_instance, false).solve_newton_raphson();
   catch
       xrs(:,t) = nan(4,1);
   end
end

static_mean = mean(xrs,2);
static_cov = cov(xrs');

%%

number_of_samples = 1000;
number_per_group = 7;
 
t = floor(n_times/2);
gnsslogdata_instance = gnsslogdata_nonan(gnsslogdata_nonan.TimeNanos == times(t),:);
prg = PsuedoRangeGroupGNSSLog(gnsslogdata_instance, false);

figure;
hold on;
plot(xrs(1,t),xrs(2,t),'b*');

for n = 7:12
    brpg = BootstrapPsuedoRangeGroupSet(prg, number_of_samples, n);
    xrs_bs = brpg.solve_each_newton_raphson();
    brpg.plot_percentile_ellipse(0.05,0.95);   
    
    if n == 7   
        plot(xrs_bs(1,:),xrs_bs(2,:),'ro');
    end
end


[xr_instance, n, ~] = prg.solve_newton_raphson();
[S,PDOP,TDOP,GDOP] = prg.DOPcalcs(xr_instance,  sqrt(mean(n.^2)));
e = error_ellipse(S(1:2,1:2), xr_instance(1:2),'conf',0.90);
e.LineWidth = 8;

e = error_ellipse(static_cov(1:2,1:2),static_mean(1:2),0.90);
e.LineWidth = 4;

xlim(-2699975 + [-150,150]);
ylim(-4292640 + [-150,150]);
axis equal
title('Comparison of Uncertainty Strategies Oval');
xlabel('x ECEF');
ylabel('y ECEF');
legend('raw','BS7','BS7 samples','BS8','BS9','BS10','BS11','BS12','DOP','static','location','best');

saveas(gcf,'Oval_1_BS.png');
