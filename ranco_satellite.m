function [solution_ranco, satellite_exclusion_solution] = ranco_satellite(gnsslogdata)

% create PsuedoRangeGroup object
prg = PsuedoRangeGroupGNSSLog(gnsslogdata,false);

%% Compute RANCO Errors

% generate the permuation of all combinations of 4
n_selected = 10;
rprgs = RANCOPsuedoRangeGroupSet(prg, n_selected);

% compute the errors of all satellites on each subset of 4
sat_errors = rprgs.satellite_errors();

% compute the thresholds of the same
thresh = rprgs.threshold();

% array with only errors that are outliers
outliers_indices = abs(sat_errors)>thresh;
outlier_errors = sat_errors;
outlier_errors(~outliers_indices) = nan;

%% Ranco

Ranco_outlier_frequency = sum(outliers_indices');
all_ranco_solutions = rprgs.solve_each_newton_raphson();
all_solution_ranco = all_ranco_solutions(:,Ranco_outlier_frequency == min(Ranco_outlier_frequency));
% all_x = rmoutliers(all_solution_ranco(1,:));
% all_y = rmoutliers(all_solution_ranco(2,:));
% all_z = rmoutliers(all_solution_ranco(3,:));
% all_b = rmoutliers(all_solution_ranco(4,:));
solution_ranco = [mean([all_solution_ranco(1,:)]) , mean([all_solution_ranco(2,:)]), mean([all_solution_ranco(3,:)]), mean([all_solution_ranco(4,:)])];
%% Compute which satellites are most frequently outliers
satellite_outlier_frequency = sum(outliers_indices);

%% Best Performing Combination index


[~, index] = sort(satellite_outlier_frequency);
best_satellite_combo_residual      = satellite_outlier_frequency(sort(index(1:n_selected)));

% Extract the elements of a at those indexes.
best_satellite_combo_ind = find(ismember(satellite_outlier_frequency, best_satellite_combo_residual));
best_satellite_combo = gnsslogdata(best_satellite_combo_ind,:);
prg_satellite = PsuedoRangeGroupGNSSLog(best_satellite_combo,false);
satellite_exclusion_solution = prg_satellite.solve_newton_raphson();
end