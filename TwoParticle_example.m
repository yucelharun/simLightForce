close all;
clearvars;




%Base parameters for simulation
fps = 100; %Simulation frame per second
psi = [12, -12]; %pattern motion velocity and direction
per = [6, 8];%Pattern parameter
Dia = [3.0, 4.0];%Particle diameter
Kapha = 20;%Kapha parameter corresponding to light power.

TotalSimulationTime = 20; % in second

%Create a simulation base
% = SimBase(Simulation size in pixels, Simulation size in micron, Frame Per
% Second, Total Simulation Time, Kapha parameter)
simulation = SimBase([500, 1000], [10, 20], fps, TotalSimulationTime*fps, Kapha);
    
    
%create pattern rate
exposuretime = 0.50; %second
exposurerate = [0.66, 0.33];
Imax = [255, 255];

%Desription of particles
% = particle(SimBase object, Diameter of the particle, medium refractive
% index, particle refractive index, Color for graphics
p = {}; ParticleNumber = 8;
for k = 1 : 2 : ParticleNumber
    p{k} = particle(simulation, Dia(1), 1.33, 1.59, 'r');
    p{k + 1} = particle(simulation, Dia(2), 1.33, 1.42, 'y');
end


timeseries = [0, 2.5, 5, 7.5, 10, 12.5, 15 , 17.5, 19.99];%Times to take a picture of thr simulation
%disp([num2str(psi) ' ' num2str(per) ' ' num2str(Dia)]);

%Recording information
path = 'ExampleSim';
FileName = [path '\p1_' num2str(Dia(1)) '_p2_' num2str(Dia(2)) '_fps' num2str(fps)];
vidrecord = false;
if vidrecord
    vid = VideoWriter(FileName, 'MPEG-4');
    vid.FrameRate = simulation.fps;
    open(vid);
end

%Create visualization figure that the results present.
f = figure('Units', 'centimeters', 'position', [5, 5, simulation.Lx, simulation.Ly]);

%Main simulation loop, calculation and visualization
for i = 1 : simulation.FrameNumber
    
    %Calculation
    if i == 1
        p = initpositions(p);
        I = PatternCombination(simulation, Imax, per, simulation.Time(i), psi, exposuretime, exposurerate);
        p = ForceCalcs(p, I); 
    else
        p = nextpos(p, simulation.Time(i));
        I = PatternCombination(simulation, Imax, per, simulation.Time(i), psi, exposuretime, exposurerate); %Lansdspace Field
        p = ForceCalcs(p, I);
    end
    
    %Visualization
    clf('reset')
    imagesc('XData', simulation.Sx,'YData', simulation.Sy,'CData', I);
    colormap(gray(256)); caxis([0 255])
    set(gca, 'Box', 'on')
    axis equal
    hold on
    for k = 1 : 1 : length(p)
        p{k}.plotpos(0);
    end
    xlabel('x Direction in micron')
    ylabel('y Direction in micron')
    text(p{1}.sim.Lx*0.70, p{1}.sim.Ly*0.8, [num2str(simulation.Time(i), '%.2f') ' sn'], 'Color', 'black', 'FontSize', 14, 'BackgroundColor', 'white')
    hold off
    
    drawnow();
    
    if vidrecord
        frame = getframe(gcf);
        im = frame2im(frame);
        writeVideo(vid,im)
    end

    if ~isequal(sum(timeseries == simulation.Time(i)),0)
        set(gcf, 'PaperPositionMode', 'auto');
        print([path '\fig_' num2str(simulation.Time(i)) '.eps'], '-depsc2')
    end
end


%Close and delete the regordings object.
if vidrecord
    close(vid)
    for k = 1 : 1 : ParticleNumber
        data = [p{k}.Time, p{k}.x, p{k}.y, p{k}.fx, p{k}.fy];
        xlswrite([FileName '_' num2str(k) 'th particle.xlsx'], data);
    end
end
delete(f);


function I = PatternCombination(sim, Imax, per, t, psi, exposuretime, exposurerate)
    
    I1 = Imax(1)*(cos(2*pi*(sim.Fx + (per(1)*psi(1)/2)*(t))/per(1)).^2);
    I2 = Imax(2)*(cos(2*pi*(sim.Fx + (per(2)*psi(2)/2)*(t))/per(2)).^2);
    
%     I = I2;
%     I = (exposurerate(1)/(exposurerate(1) + exposurerate(2)))*I1 + (exposurerate(2)/(exposurerate(1) + exposurerate(2)))*I2;

    if mod(t, exposuretime) < exposurerate(1)*exposuretime
        I = I1;
    else
        I = I2;
    end
end
function p = initpositions(p)
    for k = 1 : length(p)
        yenileme = 0;
        while yenileme == 0
            p{k}.initpos(0, randi([-15 15], 1, 1), randi([-7 7], 1, 1));
            yenileme = 1;
            if k > 1
                for j = 1 : k - 1
                    d = sqrt((p{k}.x(end) - p{j}.x(end)).^2 + (p{k}.y(end) - p{j}.y(end)).^2);
                    if d < p{k}.R/2 + p{j}.R/2 + 0.1
                        yenileme = 0;
                    end
                end
            end
        end
    end
end
function p = ForceCalcs(p, I)
    for k = 1 : length(p)
        p{k}.ForceCalc(I);
    end
end
function p = nextpos(p, t)
    for k = 1 : length(p)
        p{k}.NextPos(t);

    end
end



