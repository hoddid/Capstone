function t = TravelTime(h,e,theta1,theta2,mu)
%time it took to get to set location
%NOTE: E is in radians
T = 2*pi*(h/sqrt(1 - e^2))^3/mu^2; %period of orbit in seconds
x = sqrt((1 - e)/(1 + e))*tan(theta1/2); %intermediate calculation
E = 2*atan(x); %the big E in kepler's equation
Me = E - e*sin(E); %mean anomaly
t1 = Me*T/(2*pi); %travel time in seconds to theta1 in sec

x = sqrt((1 - e)/(1 + e))*tan(theta2/2); %intermediate calculation
E = 2*atan(x); %the big E in kepler's equation
Me = E - e*sin(E); %mean anomaly
t2 = Me*T/(2*pi); %travel time in seconds to second location

t = t2 - t1; %time between 2 locations
end