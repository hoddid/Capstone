%% Description
%This function is used to calculate the disturbance torques caused by
%various phenomona on the spacecraft, including gravity-gadient and solar
%pressure. It also calculates the total amount of fuel that would be
%required to conteract that torque (including de-saturating the reaction
%wheels), as well as design a magnetorquer to counteract control the craft
%during jovian orbit.
%% Disturbance torque calculator
clear;close all;clc;

%% General Parameters
N = 10000; %number of samples to run
muSun = 1.3271244e11;
muEarth = 3.9860045e5;
muJupiter = 1.26686545e9;
muMars = 4.282837e4;
esc = 0.601470573561797;

ef= 0.655452978476002;

theta_intersect1 = 1.485638084937971;
theta_intersect2 = 1.287134158073098;
hsc =5.638717998683439e+09;
hf = 5.983015618347387e+09;
Ixx = 37648020.309663706;
Iyy = 10004089.107927233;
Izz = 38509783.135031499;
I = [Ixx,Iyy,Izz]./1e4;
I1 = I(1);
I2 = I(2);
I3 = I(3);

%% Find solar rad vs time
Theta1 = linspace(0,theta_intersect1,N);
Theta2 = linspace(theta_intersect2,pi,N);
Rt = hsc^2/muSun *1./(1 + esc*cos(Theta1));
Rf = hf^2/muSun *1./(1 + ef*cos(Theta2)); 
T1 = [];
T2 = [];
for i = 1:N
    T1(i) = TravelTime(hsc,esc,0,Theta1(i),muSun);
    T2(i) = TravelTime(hf,ef,theta_intersect2,Theta2(i),muSun);
end

T3 = T2+T1+8.058301e6;
time = [T1,T3];%,T3(end)+220752000];
dist = [Rt,Rf];

%% Find Jupiter rad vs time 
hfinal = 3.920475810234752e+06;
efinal = 0.969669803104752;
Theta1 = linspace(0,pi,N);
RJup = hfinal^2/muJupiter *1./(1 + efinal*cos(Theta1));
TJup = [];
for i = 1:N
    TJup(i) = TravelTime(hfinal,efinal,0,Theta1(i),muJupiter);
end
TJextend = [];
adds = zeros(1,length(TJup));
for i = 1:length(TJup)-1
    adds(i+1) = TJup(i+1)-TJup(i);
end
RJ = [RJup,flip(RJup)];
TJ = [TJup,flip(flip(TJup)+adds)];%TJextend];
figure(1)
plot(TJup,RJup,'.');
title('Distance from Jupiter V Time in orbit')
xlabel('Time (sec)')
ylabel('Distance (m)')

%% Solar Torque transfer
distPlanets = [1;1.52;3;5.2];
distrs = distPlanets.*1.496e8;
press2 = @(r) 1.03e11/(r^2);
dishA = .392975; % m^2
thrusterA = .60544459;% m^2
A= thrusterA-dishA;
Torques_solar = [];%zeros(LOM);
t = 0;
for t = 1:length(time)
    Fabs = press2(dist(t))*A*cosd(0);
    Fref = 2*press2(dist(t))*A*0.5*cosd(0);
    T = (Fref+Fabs)*2.3/2;
    Torques_solar = [Torques_solar,T];
end
figure(2)
plot(dist,Torques_solar)
title('Torque from Solar Pressure')
xlabel('Mission Time (Seconds)')
ylabel('Torque (N-M)')
trans_solar_T = trapz(time,Torques_solar);
fprintf('The total Solar angular momentum imparted during the transfer is %3f N-M-S \n',trans_solar_T);
%% Solar torque in jupiters orbit
T3 = 8.058301e6;
time2 = [T3(end),T3(end)+220752000];
dist2 = [Rf(end),Rf(end)];
dishA = .392975; % m^2
thrusterA = .60544459;% m^2
A= thrusterA-dishA;
Torques_solar_juip = [];%zeros(LOM);
for t = 1:length(time2)
    Fabs = .5*press2(dist2(t))*A*cosd(0);
    Fref = 2*press2(dist2(t))*A*0.5*cosd(0);
    T = (Fref+Fabs)*2.3/2;
    Torques_solar_juip = [Torques_solar_juip,T];
end
solar_t_mission = trapz(time2,Torques_solar_juip);
fprintf('The total Solar angular momentum imparted during the Jupiter orbits is %3f N-M-S\n',solar_t_mission);

%% Gravity gradient torques
Marscount = 0;
Jupcount = 0;
Earthcount =0;
Suncount = 0;
Torques_grav_Sun = [];
Torques_grav_Jup = [];
TimeJup = [];
Jup_dist = [];
Torques_grav_Mars = [];
TimeMars = [];
Mars_dist = [];
Torques_grav_Earth = [];
TimeEarth = [];
Earth_dist = [];
for t = 1:length(time)
    %if in jupiter SOI
    if (dist(t) < 777920000+48.2e6) &&  (dist(t) > 777920000-48.2e6)
        Mu = muJupiter*1e9;
        dist_local = (dist(t)-777920000)*1e3;
        if dist_local > 0
            Jup_dist = [Jup_dist,dist_local];
            T_grav = (3*Mu)/(2*(dist_local)^3)*abs(I(3)-I(2));
            %Ensure that we never accidientally calculate from
            %inside of jupiters radius
            if dist_local < 69.91*1e6
                if max(Torques_grav_Jup)>0
                    T_grav = Torques_grav_Jup(end);
                else
                    T_grav = 0;
                end
            end
            Torques_grav_Jup = [Torques_grav_Jup,T_grav];
            TimeJup = [TimeJup,time(t)];
        end
        %if in mars SOI
    elseif (dist(t) < 227392000+.576e6) &&  (dist(t) > 227392000-.576e6)
        Mu = muMars*1e9;
        dist_local = abs(dist(t)-227392000)*1e3;
        if dist_local > 0
            T_grav = max((3*Mu)/(2*(dist_local)^3)*abs(I(3)-I(2)));
            if dist_local < 3.389*1e6
                T_grav = Torques_grav_Mars(end);
            end
            Torques_grav_Mars = [Torques_grav_Mars,T_grav];
            TimeMars = [TimeMars,time(t)];
            Mars_dist = [Mars_dist,dist_local];
        end 

    %if in earth SOI
    elseif (dist(t) < 149600000+.924e6) &&  (dist(t) > 149600000-.924e6)
        Mu = muEarth*1e9;
        dist_local = abs(dist(t)-149600000)*1e3;
        if dist_local > 0
            T_grav = max((3*Mu)/(2*(dist_local)^3)*abs(I(3)-I(2)));
            if dist_local < 6.378*1e6
                if t<5
                    T_grav = 0;
                else
                    T_grav = Torques_grav_Earth(end);
                end
            end
            Torques_grav_Earth = [Torques_grav_Earth,T_grav];
            TimeEarth = [TimeEarth,time(t)];
            Earth_dist = [Earth_dist,dist_local];
        end
    end
    % Else must be in sun SOI
    muSun_m = muSun*1e9;
    T_grav_sun = max((3*muSun_m)/(2*(dist(t)*1e3)^3)*abs(I(3)-I(2)));
    Torques_grav_Sun = [Torques_grav_Sun,T_grav_sun];
    %Mu = muSun*1e9;
    %T_grav = max(((3*Mu)/(2*(dist(t)*1e3)^3))*abs(I(2)-I(1))*sind(theta));

end
%TG = smooth(Torques_grav)
figure(3)
plot(dist,Torques_grav_Sun)
title('     Torque from Solar Gravity Gradient (Transfer Orbit)')
xlabel('Distance from Sun (Km)')
ylabel('Torque (N-M)')
xlim([1.496e8 8e8])
total_grav_sun = trapz(time,Torques_grav_Sun);

figure(4)
plot(Jup_dist./1000,Torques_grav_Jup)
title('   Torque from Jovian Gravity Gradient (Transfer Orbit)')
xlabel('Distance from jupiter (Km)')
ylabel('Torque (N-M)')
xlim([9.55e4 5e5])

total_grav_jup = 2*trapz(TimeJup,Torques_grav_Jup./100);

figure(5)
plot(Mars_dist./1000,Torques_grav_Mars)
title('Gravity Gradient from Mars (along approach)')
xlabel('Mars Distance (Km)')
ylabel('Torque (N-M)')
xlim([0 1e5])
total_grav_mars = 2*trapz(TimeMars,Torques_grav_Mars);

figure(6)
plot(Earth_dist./1000,Torques_grav_Earth)
title('Torque from Earth Gravity Gradient (Transfer Orbit)')
xlabel('Earth Distance (Km)')
ylabel('Torque (N-M)')
xlim([6500 .5e5])
total_grav_earth = 2*trapz(TimeEarth,Torques_grav_Earth./100);

total_transfer_ang_mom = total_grav_earth + total_grav_jup+total_grav_earth+total_grav_sun+trans_solar_T;

fprintf('Total Gravity Gradient from the Sun is %3f N-M-S\n',total_grav_sun);
fprintf('Total Gravity Gradient from Jupiter (along approach) is %3f N-M-S\n',total_grav_jup);
fprintf('Total Gravity Gradient from the Mars (along approach) is %3f N-M-S\n',total_grav_mars);
fprintf('Total Gravity Gradient from the Earth is %3f N-M-S \n\n',total_grav_earth);
fprintf('Total angular momentum imparted during the transfer is %3f N-M-S\n',total_transfer_ang_mom);

%% Orbiting Jupiter Gravity Gradient
Mu = muJupiter*1e9;
Jup_dist = RJup.*1000 +69.91*1e6;
Torques_Jup_orbit = [];%zeros(length(Jup_dist),1);
for i = 1:length(Jup_dist)
    %The spacecraft will spend approximatly a 15% of its time in the 
    %worst orientation for grav-grad
    T_grav1 = max((3*Mu)/(2*(Jup_dist(i))^3)*abs(I(3)-I(2)))*.15;
    T_grav2 = max((3*Mu)/(2*(Jup_dist(i))^3)*abs(I(3)-I(1)))*.85;
    T_grav = T_grav1+T_grav2;
    %Ensure that we never accidientally calculate from
    %inside of jupiters radius
    if Jup_dist(i) < 70*1e6
        if max(Torques_Jup_orbit)>0
            T_grav = Torques_Jup_orbit(end);
        else
            T_grav = 0;
        end
    end
    Torques_Jup_orbit = [Torques_Jup_orbit,T_grav];
end
figure(7)
plot(TJup,Torques_Jup_orbit)
title('Torque from Jovian Gravity Gradient (Per half orbit)')
xlabel('Distance from Jupiter (Km)')
ylabel('Torque (N-M)')
xlim([0 4000])
Grav_torque_per_orbit = 2*trapz(TJup,Torques_Jup_orbit);
ang_mom_per_orbit = (solar_t_mission/2.542336677777778e+04) + Grav_torque_per_orbit;
fprintf('Total angular momentum imparted per jovian orbit is %3f N_M_S per orbit \n',ang_mom_per_orbit)

%% Fuel Calculations
Isp = 235; % isp in Sec
mdot = 10.4/1000; %mdot max n kg/sec
g0 = 9.81;
T = Isp*mdot*g0; 
L_arm = 2.3;
T_motor = T*L_arm;
Time_burn_trans = total_transfer_ang_mom/T_motor;
Time_burn_orbit = ang_mom_per_orbit/T_motor;
M_trans = Time_burn_trans*mdot*1000;
M_orbit = Time_burn_orbit*mdot*1000*142*2;
fprintf('Total Fuel to desaturate for the tranfer is: %3f \n',M_trans)
fprintf('Total Fuel to desaturate during the duration of the orbiting is: %3f \n',M_orbit)

%% Desaturate with Mag Torque Tubes
B0 = 4.4e-4; %teslass
Turns = 1000;
current = 1;
r_tube = 7.25/100;
A = pi*r_tube^2;
R0 = 69.91*1e6; %Jupiter Radius in m
L = deg2rad(45);%orbit arount 45 degrees offset from upiter equator
BJup = @(r) ((B0*R0^3)/r^3)*sqrt(3*sin(L)+1);
Jup_dist = RJup.*1000 +69*1e6;
Torques_tube = [];%zeros(length(Jup_dist),1);
Blocals = [];
for i = 1:length(Jup_dist)
    Blocal = BJup(Jup_dist(i));
    Blocals = [Blocals,Blocal];
    T = Turns*current*A*Blocal*sind(45);
    Torques_tube = [Torques_tube,T];
end
figure(8)
plot(Jup_dist./1000,Blocals)
title('Jovian Magnetic Field Strength')
xlabel('Distance from Jupiter (Km)')
ylabel('Magnetic field (Teslas)')
xlim([7.4e4 5.7e5])
figure(9)
plot(RJup,Torques_tube)
title('Torque from Magnetic Torque Tube')
xlabel('Altitude (Km from Jupiter)')
ylabel('Avaliable Torque (N-M))')
Tot_mag_ang = 2*trapz(TJup,Torques_tube);
fprintf('Total angular momentum dissapated by torque tubes per orbit is %3f \n',Tot_mag_ang)

%% Calc power req by Torque Tubes;
revistivity = .05e-8; %revistivity Coppper at 20K
Rwire = .1/1000; %radius wire in M
Awire = pi*Rwire^2;
Lcoil = (Turns*2*pi*r_tube);
Rtube = revistivity*Lcoil/Awire;
V = current*Rtube;
P = V*current;
fprintf('This would require %3f Watts of power when fully on\n',P)
Power = P*1000/60; %Near Jupiter for approximatly 1000 seconds
P_desat = Power*(1.001177/Tot_mag_ang);
fprintf('or %3f KiloWatts-hours per orbit to desaturate \n',Power/1000)
vol = Awire*Lcoil;
mass_coil = vol*8960;
fprintf('This would require a coil of mass %3f Kg\n',mass_coil)


