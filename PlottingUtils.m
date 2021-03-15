classdef PlottingUtils
    methods(Static)     
        function plot_ellipse_2d(mean, radii)
            t = linspace(0,2*pi);
            xs = mean(1) + radii(1)*cos(t);
            ys = mean(2) + radii(2)*sin(t);
            plot(xs,ys,'LineWidth',4);
        end
        function plot_percentile_ellipse(xrs, lower_bound, upper_bound, LineSpec)
%            xrs = xrs(~isnan(xrs));
            
            
            %[U,~,~] = svd(xrs(1:2,:));
            xrs1 = xrs(1,:);
            xrs2 = xrs(2,:);
            xrs1 = xrs1(~isnan(xrs1));
            xrs2 = xrs2(~isnan(xrs2));
            [U,~,~] = svd([xrs1;xrs2]);
           
            xrs = [xrs1;xrs2];
            n = size(xrs,2);
            
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
            
            plot(p(1,:),p(2,:), LineSpec,'LineWidth',4);
        end
    end
end

