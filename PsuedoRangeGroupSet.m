classdef PsuedoRangeGroupSet < handle
    %BOOTSTRAPPSUEDORANGEGROUPSET Class pertaining to sets of
    %PseudoRangeGroups
    
    properties
        psuedoRangeGroupSet;
    end
    
    methods(Access = public)
        function obj = PsuedoRangeGroupSet(psuedoRangeGroupArray)
            obj.psuedoRangeGroupSet = psuedoRangeGroupArray;
        end
        function [xrs, ns, Gs] = solve_each_newton_raphson(obj)
            N = length(obj.psuedoRangeGroupSet);
            xrs = zeros(4,N);
            ns = zeros(obj.psuedoRangeGroupSet(1).n_measurements,N);
            [~,~,G0] = obj.psuedoRangeGroupSet(1).solve_newton_raphson();
            Gs = zeros(N, size(G0,1), size(G0,2));
            for i = 1:N
                try
                    [xr, n, G] = obj.psuedoRangeGroupSet(i).solve_newton_raphson();
                catch
                    xr = nan(4,1);
                    n = nan(obj.psuedoRangeGroupSet(1).n_measurements,1);
                    G = nan(size(G0,1), size(G0,2));
                end
                xrs(:,i) = xr;
                ns(:,i) = n;
                Gs(i,:,:) = G;
            end
        end
        function xr_mean = mean_newton_raphson_solution(obj)
            xrs = obj.solve_each_newton_raphson();
            xr_mean = mean(xrs,2);
        end
        function xr_cov = covariance_newton_raphson_solution(obj)
            xrs = obj.solve_each_newton_raphson();
            xr_cov = cov(xrs');
        end
        function plot_percentile_ellipse(obj, lower_bound, upper_bound)
            xrs = obj.solve_each_newton_raphson();
            
            xrs = xrs(:,~isnan(xrs(1,:)));
            
            n = size(xrs,2);
            
            [U,~,~] = svd(xrs(1:2,:));
            
            dxrs = xrs - mean(xrs,2);
            dxrs = dxrs(1:2,:);
            
            p1_projection = sum(dxrs.*U(:,1),1);
            p1_sorted = sort(p1_projection);
            p1_lower = p1_sorted(floor(lower_bound*n));
            p1_upper = p1_sorted(ceil(upper_bound*n));
            p1_radius = (p1_upper - p1_lower)/2;
            v1 = p1_radius * U(:,1);
            
            p2_projection = sum(dxrs.*U(:,2),1);
            p2_sorted = sort(p2_projection);
            p2_lower = p2_sorted(floor(lower_bound*n));
            p2_upper = p2_sorted(ceil(upper_bound*n));
            p2_radius = (p2_upper - p2_lower)/2;
            v2 = p2_radius * U(:,2);
            
            t = linspace(0,2*pi);
            p = v1 * cos(t) + v2*sin(t) + mean(xrs(1:2,:),2);
            
            plot(p(1,:),p(2,:),'LineWidth',4);
        end
        function Gs = geometry_matrices(obj)
            number_per_group = size(obj.psuedoRangeGroupSet(1).geometry_matrix([0;0;0;0]),1);
            n = length(obj.psuedoRangeGroupSet);
            Gs = zeros(n, number_per_group, 4);
            for i = 1:n
                x0 = obj.psuedoRangeGroupSet(i).solve_newton_raphson();
                Gs(i, :, :) = obj.psuedoRangeGroupSet(i).geometry_matrix(x0);
            end
        end
    end
end

