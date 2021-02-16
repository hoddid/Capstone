close all;clear;clc;
%% Run Control Model:
%First, Declare the initial conditions:
%Controller can handle 5 rad/sec before it has to use 
%propellant

%Declare Angular Velocities (Rad/sec)
ominit1 = 2;
ominit2 = 3;
ominit3 = -4;
%Declare Angular Position (rad)
EAinit1 = .5;
EAinit2 = 1;
EAinit3 = -1.3;


%declare desired angles and speeds
EAd1 = 0;
EAd2 = 0;
EAd3 = 0;
omd1 = 0;
omd2 = 0;
omd3 = 0;

% Input Moments of inertia

Ixx = 37648020.309663706;
Iyy = 10004089.107927233;
Izz = 38509783.135031499;
I = [Ixx,Iyy,Izz]./1e4;
I1 = I(1);
I2 = I(2);
I3 = I(3);



%Calculate the maximum torue the system can use:
T_wheels = 65; %Wheels can do 65 N-M on their own.
Thrust = 22.4; % thusters can produce 22.4 N each
boom_length = 2.3; % Boom is 2.3 meters from center 
thruster_torque = 2*Thrust*boom_length + sind(45)*Thrust;
Max_torque = thruster_torque + T_wheels;

%Determine how long the simulation will run for
Tsim = 5;
fprintf('Running Simulation, this may take a few moments.\n')
% Run the simulation
sim('sliding_lin_working.slx')

%plot things
figure(1)
hold on
plot(time,om1)
plot(time,om2)
plot(time,om3)
title('Angular Velocity Response')
ylabel('Radians Per second')
xlabel('Time (seconds)')
legend('Theta_1','Theta_2','Theta_3')
hold off


figure(2)
hold on
plot(time,ea1)
plot(time,ea2)
plot(time,ea3)
title('Angular Position Response')
ylabel('Radians')
xlabel('Time (seconds)')
legend('Theta_1','Theta_2','Theta_3')
hold off

figure(3)
hold on
plot(time,M1)
plot(time,M2)
plot(time,M3)
title('Control Moments Applied')
legend('Theta_1','Theta_2','Theta_3')
ylabel('Moment (N-M)')
xlabel('Time (seconds)')