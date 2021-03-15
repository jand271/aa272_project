function [solution_ranco, solution_satellite_exclusion] = ranco_satellite(gnsslogdata)

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

% compute the thresholds of the same
% thresh = rprgs.threshold();

% array with only errors that are outliers
%outliers_indices = abs(sat_errors)>thresh;

upper = prctile(sat_errors(:),80).*ones(size(sat_errors,1),1);
lower = prctile(sat_errors(:),20).*ones(size(sat_errors,1),1);
outliers_indices = (sat_errors > upper) | ((sat_errors < lower));

outlier_errors = sat_errors;
outlier_errors(~outliers_indices) = nan;

% figure
% hold on
% plot(sat_errors,'o')
% plot(upper,'-','LineWidth',2)
% plot(lower,'-','LineWidth',2)
% hold off

%% Ranco

Ranco_outlier_frequency = sum(outliers_indices');
all_ranco_solutions = rprgs.solve_each_newton_raphson();
all_ranco_solutions = all_ranco_solutions(:,Ranco_outlier_frequency == min(Ranco_outlier_frequency));
solution_ranco = [mean(all_ranco_solutions(1,:)) , mean(all_ranco_solutions(2,:)), mean(all_ranco_solutions(3,:)), mean(all_ranco_solutions(4,:))];
%% Compute which satellites are most frequently outliers
satellite_outlier_frequency = sum(outliers_indices);

%% Best Performing Combination index


[~, index] = sort(satellite_outlier_frequency);
best_satellite_combo_residual      = satellite_outlier_frequency(sort(index(1:n_selected)));

% Extract the elements of a at those indexes.
best_satellite_combo_ind = find(ismember(satellite_outlier_frequency, best_satellite_combo_residual));
best_satellite_combo = gnsslogdata(best_satellite_combo_ind,:);
prg_satellite = PsuedoRangeGroupGNSSLog(best_satellite_combo,false);
solution_satellite_exclusion = prg_satellite.solve_newton_raphson();
end