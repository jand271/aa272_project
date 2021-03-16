function [solution_ranco, solution_satellite_exclusion,multipath_check_se,multipath_check_ranco] = ranco_satellite(gnsslogdata)

% create PsuedoRangeGroup object
warning('off')
prg = PsuedoRangeGroupGNSSLog(gnsslogdata,false);

%% Compute RANCO Errors

% base line paramenter: 
% n_selection = 8
% upper = 80
% lower = 20

% generate the permuation of all combinations of 4
n_selected = 8;
rprgs = RANCOPsuedoRangeGroupSet(prg, n_selected);

% compute the errors of all satellites on each subset of 4
sat_errors = rprgs.satellite_errors();

% nchoosek(1:psuedoRangeGroup.n_measurements, number_per_group);

% compute the thresholds of the same
 thresh = rprgs.threshold();

% array with only errors that are outliers
outliers_indices = abs(sat_errors)>thresh;

% upper = prctile(sat_errors(:),80).*ones(size(sat_errors,1),1);
% lower = prctile(sat_errors(:),20).*ones(size(sat_errors,1),1);
% outliers_indices = (sat_errors > upper) | ((sat_errors < lower));

outlier_errors = sat_errors;
outlier_errors(~outliers_indices) = nan;
% 
% figure
% hold on
% plot(sat_errors,'o')
% plot(upper,'k-','LineWidth',2)
% plot(lower,'k-','LineWidth',2)
% % legend('Satellite Errors','Upper Bound','Lower Bound')
% ylabel('satellite errors (m)')
% xlabel('combinations')
% title('Satellite Errors with thresholds for all combinations for the first time instant of Oval')
% hold off

% figure
% hold on
% plot(sat_errors,'o')
% plot(-thresh,'r-','LineWidth',2)
% plot(thresh,'r-','LineWidth',2)
% % legend('Satellite Errors','Upper Bound','Lower Bound')
% ylabel('satellite errors (m)')
% xlabel('combinations')
% title('Satellite Errors with Ranco thresholds for all combinations for the first time instant of Oval')
% hold off


%% Ranco

Ranco_outlier_frequency = sum(outliers_indices');

all_ranco_solutions = rprgs.solve_each_newton_raphson();


all_ranco_solutions = all_ranco_solutions(:,Ranco_outlier_frequency == min(Ranco_outlier_frequency));

solution_ranco = [mean(all_ranco_solutions(1,:)) , mean(all_ranco_solutions(2,:)), mean(all_ranco_solutions(3,:)), mean(all_ranco_solutions(4,:))];

% ranco multipath detector
combination_index = rprgs.combo_ind;

combination_index = combination_index';

ranco_inlier_index = combination_index(:,Ranco_outlier_frequency == min(Ranco_outlier_frequency));
ranco_outlier_index = combination_index(:,Ranco_outlier_frequency ~= min(Ranco_outlier_frequency));

ranco_inlier_index = ranco_inlier_index(:);
inlier_gnss = gnsslogdata(ranco_inlier_index,:);
percentage_inlier_multipath = sum(inlier_gnss.MultipathIndicator)/length(ranco_inlier_index);

ranco_outlier_index = ranco_outlier_index(:);
outlier_gnss = gnsslogdata(ranco_outlier_index,:);
percentage_outlier_multipath = sum(outlier_gnss.MultipathIndicator)/length(ranco_outlier_index);


% if the multipath detect 1's percentage for inlier is larger than outlier,
% we indicate the check fails, and output 0
multipath_check_ranco = 1;
if percentage_outlier_multipath < percentage_inlier_multipath
    multipath_check_ranco = 0;
end
%% Compute which satellites are most frequently outliers
satellite_outlier_frequency = sum(outliers_indices);

%% Best Performing Combination index


[~, index] = sort(satellite_outlier_frequency);
best_satellite_combo_residual      = satellite_outlier_frequency(sort(index(1:n_selected)));

% Extract the elements of a at those indexes.
best_satellite_combo_ind = find(ismember(satellite_outlier_frequency, best_satellite_combo_residual));

excluded_satellite_combo_ind = find(~ismember(satellite_outlier_frequency, best_satellite_combo_residual));

best_satellite_combo = gnsslogdata(best_satellite_combo_ind,:);

multi_path_check_best_se = gnsslogdata(best_satellite_combo_ind,:).MultipathIndicator;

multi_path_check_worst_se = gnsslogdata(excluded_satellite_combo_ind,:).MultipathIndicator;

multipath_check_se = 1; % one means checks, 0 means fails

if sum(multi_path_check_best_se == 1)>sum(multi_path_check_worst_se == 1)
    rmultipath_check_se = 0;
end
prg_satellite = PsuedoRangeGroupGNSSLog(best_satellite_combo,false);
solution_satellite_exclusion = prg_satellite.solve_newton_raphson();
end