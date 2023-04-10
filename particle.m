classdef particle < handle
    properties
        sim
        R %Diameter in micron
        Kb = 1.38064852e-23*1e6*1e15; % m2 kg s-2 K-1 to m N K-1 to micron fN K-1
        T = 300; %Temperature K
        m = 11; %pg
        vis = 0.001*1e+15/((1e+6)^2); %Ns/m^2  to fN s /micronm^2
        gamma
        np = 1.52;
        nm = 1.33;
        n 
        D  %diff constant
        Time = [];
        x = [] % x Pasitions
        y = []% y positions
        fx = []%Force components
        fy = []%Force components
        color = 'r'
    end
    methods
        function obj = particle(sim, R, n1, n2, color)
            obj.sim = sim;
            obj.R = R;
            obj.nm = n1;
            obj.np = n2;
            obj.gamma = 6*pi*obj.vis*obj.R;
            obj.D = obj.Kb*obj.T/obj.gamma;
            obj.n = obj.np/obj.nm;
            obj.color = color;
        end
        function NextPos(obj, t)
            
            x = obj.x(end) + obj.fx(end)*obj.sim.Dt/obj.gamma + sqrt(2*obj.D*obj.sim.Dt)*normrnd(0,1);
            y = obj.y(end) + obj.fy(end)*obj.sim.Dt/obj.gamma + sqrt(2*obj.D*obj.sim.Dt)*normrnd(0,1);
            
            if x > obj.sim.Lx - obj.R/2
                x = obj.sim.Lx - obj.R/2;
            end
            if x < - obj.sim.Lx + obj.R/2
                x = - obj.sim.Lx + obj.R/2;
            end
            if y > obj.sim.Ly - obj.R/2
                y = obj.sim.Ly - obj.R/2;
            end
            if y < - obj.sim.Ly + obj.R/2
                y = - obj.sim.Ly + obj.R/2;
            end
            
            obj.x = [obj.x; x];
            obj.y = [obj.y; y];
            obj.Time = [obj.Time; t];

        end
        function initpos(obj, t, x0, y0)
            obj.Time = t;
            obj.x = x0;
            obj.y = y0;
            
        end
        function plotpos(obj, log)
            if log
                plot(obj.x, obj.y , obj.color)
            end
            
            th = 0:pi/50:2*pi;
            xunit = obj.R/2 * cos(th) + obj.x(end);
            yunit = obj.R/2 * sin(th) + obj.y(end);
            plot(xunit, yunit, obj.color, 'LineWidth', 2);
        end
        function plotTrajectories(obj)
            figure, plot(obj.Time, obj.x)
            figure, plot(obj.Time, obj.y)
        end
        function plotMSD(obj)
            num = length(obj.x);
            dt = obj.Time(2) - obj.Time(1);
            
            for k = 1:num/2
                MSDx(k) =mean((obj.x(k+1:1:end) - obj.x(1:1:end-k)).^2);
                MSDy(k) =mean((obj.y(k+1:1:end) - obj.y(1:1:end-k)).^2); 
                Zn(k)=k*dt;
            end
            
            figure,
            loglog(Zn(1:1:end),MSDx(1:1:end)/1e-12,'sr')
            hold on
            loglog(Zn(1:1:end),MSDy(1:1:end)/1e-12,'db')
            
        end
        function [fx, fy] = ForceCalc(obj, I)
            
            rr = obj.R/2;
            fx = 0; fy = 0;
            dr = rr/50;
            r = 0 + dr/2 : dr : obj.R/2 - dr/2;
            
            dw = 2*pi/360;
            w = 0 + dw/2 : dw : 2*pi - dw/2;
            
            data = [];
            for i = 1 : length(r)
                for j = 1 : length(w)
                    data = [data; r(i)*cos(w(j)), r(i)*sin(w(j)), r(i), w(j)];
                end
            end
            

            Dr = (data(:,3)/rr).*sqrt(1-(data(:,3)/obj.n/rr).^2)- ...
                (data(:,3)/obj.n/rr).*sqrt(1-(data(:,3)/rr).^2);
            Dr = real(Dr);
            
            II = interp2(obj.sim.Fx, obj.sim.Fy, I, obj.x(end) + data(:,1), obj.y(end) + data(:,2));%, 'spline'
            
            Dx = II.*Dr.*cos(data(:,4)).*data(:,3).*dr.*dw;
            Dy = II.*Dr.*sin(data(:,4)).*data(:,3).*dr.*dw;

%             a = sum(data(:,3).*dr.*dw)
%             aa = (pi*(obj.R/2)^2)

            fx = obj.sim.kapha*sum(Dx);
            fy = obj.sim.kapha*sum(Dy);
            
            obj.fx = [obj.fx; fx];
            obj.fy = [obj.fy; fy];
        end
        function f = F_Harmonics(obj, k, x)
            f = - k*x;
        end
        function f = F_Kramers(obj, a, b, c, x)
            f = - a*(x.^3) + b*x + c;
        end
    end
end