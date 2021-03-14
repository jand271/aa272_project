classdef PlottingUtils
    methods(Static)     
        function plot_ellipse_2d(mean, radii)
            t = linspace(0,2*pi);
            xs = mean(1) + radii(1)*cos(t);
            ys = mean(2) + radii(2)*sin(t);
            plot(xs,ys,'LineWidth',4);
        end
    end
end

