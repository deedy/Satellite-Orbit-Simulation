%ORBIT ANALYZER
%Violet Attitude Control Subsystem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Converts the State Vectors -
%1. Position Vector [Rx Ry Rz]
%2. Velocity Vector [Vx Vy Vz]
% into Kepler's Orbital Elements -
%1. Length of Semi-major Axis (a)
%2. Eccentricity (e)
%3. Inclination (inc)
%4. Right Ascension of Ascending Node (RAAN)
%5. Argument of Perigee (w)
%6. Initial Mean Anomaly (M0)
% after which it visually displays the orbit of the diagram in 3D space.
% The Blue Sphere denotes the Earth
% The Green Sphere denotes the initial position of the satellite
% The Red Spheres denote the positions of the satellites at fixed time
% intervals
% The dotted green line denotes the initial velocity vector of the
% satellite
% The three orthogonal black lines denote the 3 positive ECI axes
% 
% At the end, the altitude ranges are displayed, as well as the orbit type
% - LEO or MEO
% The simulation only takes into account the the influence of the Earth
%CODED BY DEBARGHYA DAS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CUSTOMIZABLE STATE VECTORS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SAMPLE VALUES:
%MOLNIYA ORBIT
%r= [0 -7000 -14000];
%v= [2 1 -0.5];
%vmag= 6175;
%runspeed= 300;
%GEOSTATIONARY ORBIT
%r= [42241 0 0]; 
%v= [0 1 0];
%vmag= 3072;
%runspeed=1000;
%POLAR ORBIT COVERING ALL LONGITUDES IN ONE REV
%r= [42241 0 0];
%v= [0 0 1];
%vmag=3072;
%runspeed=1000;
%LOW ECCENTRICITY LOW EARTH ORBIT
%r=[7000 0 0];
%v=[0 1 1];
%vmag=7500;
%runspeed=100;
%NORMAL ORBIT
%r= [7000 0 0]; 
%v= [0 1 0.5];
%vmag=7500;
%runspeed=100;

%INITIAL POSITION VECTOR (in km)
r= [0 -7000 -14000];
%INITIAL VELOCITY DIRECTION VECTOR
v= [2 1 -0.5];
%INITIAL VELOCITY MAGNITUDE (in m/s)
vmag=6175;
v= vmag*(v/sqrt(dot(v, v))); %Velocity Vector
%RUN SPEED OF SIMULATION (Increasing this value increases the time interval
%between the plotting of the satellite)
runspeed=300;
%NUMBER OF REVOLUTIONS 
rev_no=5;
%Number of Days since Jan 1, 2000
J2000_days=103752/24; % = 4321 on 30th October, 2011 http://www.timeanddate.com/counters/year2000.html


% Constant parameters
mu = 398.6004418e12;  % Planetary gravitational constant for Earth, (mu = GMearth) (m^3/s^2)
earth_rad=6371000; % in m



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONVERTING STATE VECTORS INTO ORBITAL ELEMENTS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V=v;
r=r*1000; %Position vector conversion to meters

rmag = sqrt(dot(r, r)); %Position Magnitude
vmag = sqrt(dot(v, v)); %Velocity Magnitude

rhat = r/rmag; %Position Unit Vector
vhat = v/vmag; %Velocity Unit Vector

hv = cross(r, v); %Angular Momentum Vector
hmag = sqrt(dot(hv, hv)); %Angular Momentum Magnitude
hhat = hv/hmag; %Angular Momentum Unit Vector

%Eccentricity Vector
vtmp = v / mu;
ecc = cross(vtmp, hv);
ecc = ecc - rhat;

%SEMIMAJOR AXIS (a)
a = 1 / (2 / rmag - vmag * vmag / mu);

p = hhat(1) / (1 + hhat(3));
q = -hhat(2) / (1 + hhat(3));
const1 = 1 / (1 + p * p + q * q);
fhat(1) = const1 * (1 - p * p + q * q);
fhat(2) = const1 * 2 * p * q;
fhat(3) = -const1 * 2 * p;
ghat(1) = const1 * 2 * p * q;
ghat(2) = const1 * (1 + p * p - q * q);
ghat(3) = const1 * 2 * q;
h = dot(ecc, ghat);
xk = dot(ecc, fhat);
x1 = dot(r, fhat);
y1 = dot(r, ghat);

%ECCENTRICITY (e) %0<e<1
e = sqrt(h * h + xk * xk)

%INCLINATION (inc) %in rad
inc = 2 * atan(sqrt(p * p + q * q));

xlambdat = atan3(y1,x1);

%RIGHT ASCENSION OF ASCENDING NODE (RAAN) %in rad
if (inc > 0.00000001)
    RAAN = atan3(p,q);
else
   RAAN = 0;
end

%ARGUMENT OF PERIGEE (w) %in rad
if (e > 0.00000001)
   w = atan3(h,xk)-RAAN; 
else
   w = 0;
end

%True Anomaly %in rad
v = xlambdat - RAAN - w;

%MEAN ANOMALY (M0)
M0 = 2*atan(sqrt((1-e)/(1+e))*tan(v/2)) - e*sqrt(1-e^2)*sin(v)/(1+e*cos(v)); %in rad



%ERROR if Launch position is inside the Earth
if (sqrt (r(1)*r(1)+r(2)*r(2)+r(3)*r(3)) <= 6731000)
    blast (r(1), r(2), r(3), 2000000);
    ErrorMsg='Launch Position Inside Earth'
    return;
end


%Final Adjustments to RAAN
RAAN=pi/2-RAAN;
if RAAN < 0
    RAAN = RAAN + 2*pi;
end


%Final Adjustments to Argument of Perigee
w=2*pi-w;



%Final Adjustments to Initial Mean Anomaly
if M0<0
    M0=-M0;
end
    E=M0;
    for i=1:5
        E = E + (M0 + e*sin(E) - E)/(1 - e*cos(E));
    end
    v= 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    R = a*(1-e*cos(E));
    Xeci = R*(cos(w + v)*cos(RAAN) - sin(w+v)*sin(RAAN)*cos(inc));
    Yeci = R*(cos(w + v)*sin(RAAN) + sin(w+v)*cos(RAAN)*cos(inc));
    Zeci = R*(sin(w + v)*sin(inc));
    c=0;
    while abs(r(1)-Xeci)>100 && abs(r(2)-Yeci)>100 && abs(r(3)-Zeci)>100 && c<15
        if c~=0
            if (c<3)
                M0=2*pi-M0;
            end
            if (c<5 && c>=3)
                M0=-M0;
            end
            if (c<10 && c>=5)
                M0=M0+(pi/2);
                if (c==9)
                    M0=iniM0;
                end
            end
            if (c>=10 && c<15)
                M0=M0-(pi/2);
                if (c==15)
                    M0=iniM0;
                end
            end 
        end
        E=M0;
        for i=1:5
            E = E + (M0 + e*sin(E) - E)/(1 - e*cos(E));
        end
        v= 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
        R = a*(1-e*cos(E));
        b  = a*sqrt(1-e^2);
        Xeci = R*(cos(w + v)*cos(RAAN) - sin(w+v)*sin(RAAN)*cos(inc));
        Yeci = R*(cos(w + v)*sin(RAAN) + sin(w+v)*cos(RAAN)*cos(inc));
        Zeci = R*(sin(w + v)*sin(inc));
        c=c+1;
end
if (M0>2*pi)
    M0=M0-2*pi;
end
if (M0<0)
    M0=M0+2*pi;     
end
  
%Orbital Period 
orbital_period = sqrt((a*a*a*4*pi*pi)/mu)
%DISPLAYING ORBITAL ELEMENTS
a_km= a/1000 %Semi-major axis (in m)
b_km  = a*sqrt(1-e^2)/1000
e %Eccentricity
inc_deg= inc*180/pi %Inclination (in degrees)
RAAN_deg= RAAN*180/pi %Right Ascension of Ascending Node (in degrees)
w_deg= w*180/pi %Argument of Perigee (in degrees)
M0_deg= M0*180/pi %Mean Anomaly (in degrees)

%ERROR if eccentricity is greater than supported
if (e<0)
    if (e>0.95)
        Errormsg='Eccentricity greater than supported'
    end
    if (e>1)
        Errormsg='Eccentricity invalid'
    end
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DRAWING THE STATIC VISUALIZATION COMPONENTS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initializing the Drawing Space
close all
set(gcf,'Menubar','default','Name','Orbit Visualization', ... 
    'NumberTitle','off','Position',[10,350,750,750], ... 
    'Color',[0.38 0.26 0.67]); 
lim=(1+e)*a;%Setting the limits of the graph
clf
axis([-lim, lim, -lim, lim, -lim, lim])	
view(150,15) 
axis equal
shg
hold on
grid on
title('Orbital Visualization');

%Plotting the Earth
equat_rad=6378137.00;
    polar_rad=6356752.3142;
    [xx yy zz]=ellipsoid (0,0,0,equat_rad, equat_rad, polar_rad);
    load('topo.mat','topo','topomap1');
    topo2 = [topo(:,181:360) topo(:,1:180)];
    pro.FaceColor= 'texture';
    pro.EdgeColor = 'none';
    pro.FaceLighting = 'phong';
    pro.Cdata = topo2;
   earth= surface(xx,yy,zz,pro)
    colormap(topomap1)
omega_earth = 7.292115855377074e-005; % (rad/sec)  
Go = 1.727564365843028; % (rad) http://www.amsat.org/amsat/articles/g3ruh/106.html 
GMST = Go + omega_earth*86400*(J2000_days + 0.5);
GMST = GMST - 2*pi*floor(GMST/(2*pi));
GMST_deg=GMST*(180/pi)
rotate (earth, [0 0 1], GMST_deg);
Xaxis= line([0 lim],[0 0],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
Yaxis= line([0 0],[0 lim],[0 0],'Color', 'red', 'Marker','.','LineStyle','-');
rotate (Xaxis, [0 0 1], GMST_deg);
rotate (Yaxis, [0 0 1], GMST_deg);
%Sun=light('Position',[1 0 0],'Style','infinite');


%Plotting the ECI Axes
line([0 lim],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
line([0 0],[0 lim],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
line([0 0],[0 0],[0 lim],'Color', 'black', 'Marker','.','LineStyle','-')

%Plotting Initial Velocity Vector
line([r(1) r(1)+1000*V(1)],[r(2) r(2)+1000*V(2)],[r(3) r(3)+1000*V(3)],'Color', 'green','Marker','.','LineWidth', 2, 'MarkerSize', 8,'LineStyle','-');

%Plotting the initial poisition of the satellite
plot3 (r(1), r(2), r(3),'o', 'MarkerEdgeColor', 'black','MarkerFaceColor','green','MarkerSize', 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DRAWING THE DYNAMIC VISUALIZATION COMPONENTS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=0;
long=1:rev_no*ceil(orbital_period/runspeed);
%Plotting the movement of the satellite
Xcoord(1)=r(1);
Ycoord(1)=r(2);
Zcoord(1)=r(3);
%Sun=light('Position',[1 0 0],'Style','infinite', 'Visible', 'on');

for time = 1:rev_no*ceil(orbital_period/runspeed)
    
    k=k+1;
    %Computing Eccentric Anomaly
    E=M0;
    for i=1:5
        E = E + (M0 + e*sin(E) - E)/(1 - e*cos(E));
    end
   
    %Computing the True Anomaly
    v= 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    
    %Computing 'r' in polar coordinates
    r = a*(1-e*cos(E));
    GMST=GMST+(runspeed/86400)*2*pi;
    %Computes the Cartesian Co-ordinates in ECI frame from 'r' and orbital
    %elements
    Xeci = r*(cos(w + v)*cos(RAAN) - sin(w+v)*sin(RAAN)*cos(inc));
    Yeci = r*(cos(w + v)*sin(RAAN) + sin(w+v)*cos(RAAN)*cos(inc));
    Zeci = r*(sin(w + v)*sin(inc));
    rotate (earth, [0 0 1], (runspeed)*(360/86400))
    
    rotate (Xaxis, [0 0 1], (runspeed)*(360/86400))
    rotate (Yaxis, [0 0 1], (runspeed)*(360/86400))

     
    
    %Drawing the red sphere
    array(k)=plot3 (Xeci, Yeci, Zeci,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 6);
    position(k)=line([0 Xeci],[0 Yeci], [0 Zeci],'Color', 'yellow', 'LineWidth', 2);
    
    if (k~=1)
    set (array(k-1), 'Visible', 'off');
    set (position(k-1), 'Visible', 'off');
    end
    if (time~=1 && time<=ceil(orbital_period/runspeed)+1)
        Xcoord(k)=Xeci;
        Ycoord(k)=Yeci;
        Zcoord(k)=Zeci;
        line([Xcoord(k-1) Xcoord(k)],[Ycoord(k-1) Ycoord(k)], [Zcoord(k-1) Zcoord(k)],'Color', 'black', 'LineWidth', 2);
    end
    
    
  
    if (GMST>2*pi)
        GMST=GMST-2*pi
    end
    
    lat(k)=atan(Zeci/sqrt(Xeci*Xeci+Yeci*Yeci))*(180/pi);
    ECIX=[cos(GMST) sin(GMST) 0];
    Pos=[Xeci Yeci 0];
    %Posmag = sqrt(dot(Pos, Pos));
    cvec = cross(ECIX,Pos);
    angleyz = mod(sign(dot([0 0 1],cvec))*atan2(norm(cvec),dot(ECIX,Pos)),2*pi);
    long(k) =(180/pi)* angleyz;
    
        
    %Pause
    pause (0.01);
    
    %Blast condition
    if (sqrt (Xeci*Xeci+Yeci*Yeci+Zeci*Zeci) <= 6731000)
        blast (Xeci, Yeci, Zeci, 2000000);
        ErrorMsg='Blast Effected'
        break;
    end
    
    M0=M0+sqrt(mu/(a*a*a))*runspeed; %Updating Mean Anomaly for next iteration
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DISPLAYING ALTITUDE RANGE AND ORBIT TYPE%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

perigee=a*(1-e);
apogee=a*(1+e);
altitude_low=(perigee-earth_rad)/1000
altitude_high=(apogee-earth_rad)/1000
if altitude_low>200 && altitude_high<2000
    'Low Earth Orbit'
end
if altitude_low>2000 && altitude_high<35786
    'Middle Earth Orbit'
end
figure (2);
set(gcf,'Menubar','none','Name','Earth Track', ... 
    'NumberTitle','off','Position',[10,350,1000,500], ... 
    'Color',[0.4 0.3 0.7]); 
hold on
image([0 360],[-90 90],topo,'CDataMapping', 'scaled');
colormap(topomap1);
axis equal
axis ([0 360 -90 90]);
plot (167.717,8.717,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 10);
plot (360-76.496, 42.440,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 10);
for i=1:k
    plot (long(i),lat(i),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 6);
    if (i~=1 && abs(long(i-1)-long(i))<100)
        line([long(i-1) long(i)],[lat(i-1) lat(i)],'Color', 'red', 'LineWidth', 2);
    end
    
    pause (0.001);
end
