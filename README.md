Satellite-Orbit-Simulation
==========================

Simulates a satellite's Low Earth Orbit trajectory and additionally outputs it's Ground Track for a custom launch position and speed.

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
