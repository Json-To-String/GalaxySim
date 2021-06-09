%% Galaxy Project
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultAxesFontSize',18)
set(0,'defaultfigureposition',[0 0 700 350]); 
format compact; 
clear all; close all; clc; 

%% Tracer Stars

pC = physicsConstants();
yr2sec = 31557600; % s
mSun = 1.989e30; % kg
mGal = (10e11)*mSun; % kg
kpc = 3.0857e19; % kpc to m

initCond0 = [0 0 0 0]; % kpc kpc m/s m/s 
initCond1 = [25 175 -100e3 -300e3]; % kpc kpc m/s m/s
initCond2 = [-25 -175 100e3 300e3]; % kpc kpc m/s m/s

[xStars0, yStars0, vX0, vY0] = galaxies(initCond0, mGal, -1); % Galaxy on its lonesome
[xStars1, yStars1, vX1, vY1] = galaxies(initCond1, mGal, -1); % Parent gal w/ initial conditions
[xStars2, yStars2, vX2, vY2] = galaxies(initCond2, mGal, 1); % Intruder gal w/ vel opp parent
[xStars3, yStars3, vX3, vY3] = galaxies(initCond2, mGal, -1); % Intruder gal w/ vel along parent
[xStars4, yStars4, vX4, vY4] = galaxies(initCond2, .9*mGal, -1); % Light intruder gal w/ vel along parent
[xStars5, yStars5, vX5, vY5] = galaxies(initCond2, 1.1*mGal, -1); % Heavy intruder gal w/ vel along parent

%%
figure(10)
hold on
    plot(xStars0/kpc, yStars0/kpc, 'r*')
    plot(0,0,'k+','markersize',25)
    quiver(xStars0/kpc, yStars0/kpc, vX0/kpc, vY0/kpc,'b')
    title('Galaxy on its own')
    xlabel('x[kpc]')
    ylabel('y[kpc]')
    axis([-25 25 -25 25])

%%
figure(11)
hold on
    plot(xStars1/kpc, yStars1/kpc, 'r*')
    plot(25,175,'k+','markersize',25)
    quiver(xStars1/kpc, yStars1/kpc, vX1/kpc, vY1/kpc, 'b')
    title('Galaxy with initial conditions')
    xlabel('x[kpc]')
    ylabel('y[kpc]')
    axis([0 50 150 200])

%%
timeSpan = [0 1e9*yr2sec]; % Seconds
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
initGal1 = [xStars1 yStars1 vX1 vY1 ...
		initCond1(1)*kpc initCond1(2)*kpc initCond1(3) initCond1(4) ...
		xStars2 yStars2 vX2 vY2 ...
		initCond2(1)*kpc initCond2(2)*kpc initCond2(3) initCond2(4)];

galaxyODE1 = @(t,vals) [vals(101:150); vals(151:200); ...
			galAccel(t,vals(1:204),mGal,pC.G,vals(205:end)); ...
			vals(305:354); vals(355:404); ...
			galAccel(t,vals(205:end),mGal,pC.G,vals(1:204))];

[t1,data1] = ode45(galaxyODE1,timeSpan,initGal1,opts);

initGal2 = [xStars1 yStars1 vX1 vY1 ...
		initCond1(1)*kpc initCond1(2)*kpc initCond1(3) initCond1(4) ...
		xStars3 yStars3 vX3 vY3 ...
		initCond2(1)*kpc initCond2(2)*kpc initCond2(3) initCond2(4)];

galaxyODE2 = @(t,vals) [vals(101:150); vals(151:200); ...
			galAccel(t,vals(1:204),mGal,pC.G,vals(205:end)); ...
			vals(305:354); vals(355:404); ...
			galAccel(t,vals(205:end),mGal,pC.G,vals(1:204))];

[t2,data2] = ode45(galaxyODE2,timeSpan,initGal2,opts);

initGal3 = [xStars1 yStars1 vX1 vY1 ...
		initCond1(1)*kpc initCond1(2)*kpc initCond1(3) initCond1(4) ...
		xStars4 yStars4 vX4 vY4 ...
		initCond2(1)*kpc initCond2(2)*kpc initCond2(3) initCond2(4)];

galaxyODE3 = @(t,vals) [vals(101:150); vals(151:200); ...
			galAccel(t,vals(1:204),mGal,pC.G,vals(205:end)); ...
			vals(305:354); vals(355:404); ...
			galAccel(t,vals(205:end),mGal,pC.G,vals(1:204))];

[t3,data3] = ode45(galaxyODE3,timeSpan,initGal3,opts);


initGal4 = [xStars1 yStars1 vX1 vY1 ...
		initCond1(1)*kpc initCond1(2)*kpc initCond1(3) initCond1(4) ...
		xStars5 yStars5 vX5 vY5 ...
		initCond2(1)*kpc initCond2(2)*kpc initCond2(3) initCond2(4)];

galaxyODE4 = @(t,vals) [vals(101:150); vals(151:200); ...
			galAccel(t,vals(1:204),mGal,pC.G,vals(205:end)); ...
			vals(305:354); vals(355:404); ...
			galAccel(t,vals(205:end),mGal,pC.G,vals(1:204))];

[t4,data4] = ode45(galaxyODE4,timeSpan,initGal4,opts);

figure(12)
hold on
axis([-150 150 -200 200])
title('Equal Mass Retrograde Passage')
xlabel('x[kpc]')
ylabel('y[kpc]')
    half = 204; 
    p1 = plot(data1(1,1:50)/kpc, data1(1,51:100)/kpc, 'r.');
    p2 = plot(data1(1,half-3)/kpc, data1(1,half-2)/kpc, 'k+');
    p3 = plot(data1(1,205:254)/kpc, data1(1,255:304), 'b.');
    p4 = plot(data1(1,end-3)/kpc, data1(1,end-2)/kpc, 'k+'); 
    plot(data1(:,half-3)/kpc, data1(:,half-2)/kpc,'r')
    plot(data1(:,end-3)/kpc, data1(:,end-2)/kpc,'b')

        for N=2:length(t1)
            set(p1,'XData',data1(N,1:50)/kpc,'YData',data1(N,51:100)/kpc);
            set(p2,'XData',data1(N,half-3)/kpc,'YData',data1(N,half-2)/kpc);
            set(p3,'XData',data1(N,205:254)/kpc,'YData',data1(N,255:304)/kpc);
            set(p4,'XData',data1(N,end-3)/kpc,'YData',data1(N,end-2)/kpc);
            pause(0.0001)
        end

figure(13)
hold on
axis([-150 150 -200 200])
title('Equal Mass Direct Passage')
xlabel('x[kpc]')
ylabel('y[kpc]')
    half = 204; 
    p1 = plot(data2(1,1:50)/kpc, data2(1,51:100)/kpc, 'r.');
    p2 = plot(data2(1,half-3)/kpc, data2(1,half-2)/kpc, 'k+');
    p3 = plot(data2(1,205:254)/kpc, data2(1,255:304), 'b.');
    p4 = plot(data2(1,end-3)/kpc, data2(1,end-2)/kpc, 'k+'); 
    plot(data2(:,half-3)/kpc, data2(:,half-2)/kpc,'r')
    plot(data2(:,end-3)/kpc, data2(:,end-2)/kpc,'b')

        for N=2:length(t2)
            set(p1,'XData',data2(N,1:50)/kpc,'YData',data2(N,51:100)/kpc);
            set(p2,'XData',data2(N,half-3)/kpc,'YData',data2(N,half-2)/kpc);
            set(p3,'XData',data2(N,205:254)/kpc,'YData',data2(N,255:304)/kpc);
            set(p4,'XData',data2(N,end-3)/kpc,'YData',data2(N,end-2)/kpc);
            pause(0.0001)
        end

figure(14)
hold on
axis([-150 150 -200 200])
title('Light Intruder Direct Passage')
xlabel('x[kpc]')
ylabel('y[kpc]')
    half = 204; 
    p1 = plot(data3(1,1:50)/kpc, data3(1,51:100)/kpc, 'r.');
    p2 = plot(data3(1,half-3)/kpc, data3(1,half-2)/kpc, 'k+');
    p3 = plot(data3(1,205:254)/kpc, data3(1,255:304), 'b.');
    p4 = plot(data3(1,end-3)/kpc, data3(1,end-2)/kpc, 'k+'); 
    plot(data3(:,half-3)/kpc, data3(:,half-2)/kpc,'r')
    plot(data3(:,end-3)/kpc, data3(:,end-2)/kpc,'b')

        for N=2:length(t3)
            set(p1,'XData',data3(N,1:50)/kpc,'YData',data3(N,51:100)/kpc);
            set(p2,'XData',data3(N,half-3)/kpc,'YData',data3(N,half-2)/kpc);
            set(p3,'XData',data3(N,205:254)/kpc,'YData',data3(N,255:304)/kpc);
            set(p4,'XData',data3(N,end-3)/kpc,'YData',data3(N,end-2)/kpc);
            pause(0.0001)
        end        

figure(15)
hold on
axis([-150 150 -200 200])
title('Heavy Intruder Direct Passage')
xlabel('x[kpc]')
ylabel('y[kpc]')
    half = 204; 
    p1 = plot(data4(1,1:50)/kpc, data4(1,51:100)/kpc, 'r.');
    p2 = plot(data4(1,half-3)/kpc, data4(1,half-2)/kpc, 'k+');
    p3 = plot(data4(1,205:254)/kpc, data4(1,255:304), 'b.');
    p4 = plot(data4(1,end-3)/kpc, data4(1,end-2)/kpc, 'k+'); 
    plot(data4(:,half-3)/kpc, data4(:,half-2)/kpc,'r')
    plot(data4(:,end-3)/kpc, data4(:,end-2)/kpc,'b')

        for N=2:length(t3)
            set(p1,'XData',data4(N,1:50)/kpc,'YData',data4(N,51:100)/kpc);
            set(p2,'XData',data4(N,half-3)/kpc,'YData',data4(N,half-2)/kpc);
            set(p3,'XData',data4(N,205:254)/kpc,'YData',data4(N,255:304)/kpc);
            set(p4,'XData',data4(N,end-3)/kpc,'YData',data4(N,end-2)/kpc);
            pause(0.0001)
            if N > 4449
                break
            end
        end        
% figure(16)
% hold on
%     plot(data1(:,half-3)/kpc, data1(:,half-2)/kpc)
%     plot(data1(:,end-3)/kpc, data1(:,end-2)/kpc)
%     xlabel('x[kpc]')
%     ylabel('y[kpc]')
%     axis equal
   
%% Local Functions

function [xStars, yStars, vX, vY] = galaxies(initCond, mass, dir)
	
	pC = physicsConstants();
	kpc = 3.0857e19; % kpc to m
	r = [5 10 15 20]; 
	rAll = kpc.*r; % m

	theta1 = linspace(0, 2*pi*(r(1) / (r(1)+1)), r(1)); 
	theta2 = linspace(0, 2*pi*(r(2) / (r(2)+1)), r(2)); 
	theta3 = linspace(0, 2*pi*(r(3) / (r(3)+1)), r(3)); 
	theta4 = linspace(0, 2*pi*(r(4) / (r(4)+1)), r(4)); 

	xStars = initCond(1)*kpc + [rAll(1)*cos(theta1) rAll(2)*cos(theta2)...
        rAll(3)*cos(theta3) rAll(4)*cos(theta4)];
	yStars = initCond(2)*kpc + [rAll(1)*sin(theta1) rAll(2)*sin(theta2)...
        rAll(3)*sin(theta3) rAll(4)*sin(theta4)];

	v = @(r,theta) -theta.*sqrt((pC.G.*mass)./r);
	vX = initCond(3) - dir.*[v(rAll(1), sin(theta1)) v(rAll(2),...
        sin(theta2)) v(rAll(3), sin(theta3)) v(rAll(4), sin(theta4))];
	vY = initCond(4) + dir.*[v(rAll(1), cos(theta1)) v(rAll(2),...
        cos(theta2)) v(rAll(3), cos(theta3)) v(rAll(4), cos(theta4))];

end

function accel = starAccel(t,vals,mass,G,xGal,yGal)
	
	x = vals(1:50) - xGal;
	y = vals(51:100) - yGal;
	r = sqrt((x).^2 + (y).^2);	
	a = (-G*mass)./(r.^2);
	ax = (a.*(x))./r;
	ay = (a.*(y))./r;
	accel = [ax; ay];
end

function accel = galAccel(t,vals,mass,G,intruder)
	
	x1 = vals(1:50) - vals(end-3);
    x2 = vals(1:50) - intruder(end-3);  
	y1 = vals(51:100) - vals(end-2);
    y2 = vals(51:100) - intruder(end-2);   
	r1 = sqrt((x1).^2 + (y1).^2);
   	r2 = sqrt((x2).^2 + (y2).^2);	

	a1 = (-G*mass)./(r1.^2);
    a2 = (-G*mass)./(r2.^2);

	ax1 = (a1.*(x1))./r1;
   	ax2 = (a2.*(x2))./r2;
	ay1 = (a1.*(y1))./r1;
	ay2 = (a2.*(y2))./r2;

	ax_tot = ax1 + ax2;
	ay_tot = ay1 + ay2;

	vXGal = vals(end-1);
	vYGal = vals(end);
	rGal = sqrt((vals(end-3)-intruder(end-3))^2 + (vals(end-2)-intruder(end-2))^2);
	aGal = (-G*mass)/(rGal^2);
	aXGal = aGal*((vals(end-3)-intruder(end-3))/rGal);
	aYGal = aGal*((vals(end-2)-intruder(end-2))/rGal);
	accel = [ax_tot; ay_tot; vXGal; vYGal; aXGal; aYGal;];
end
