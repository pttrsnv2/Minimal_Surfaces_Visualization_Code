clear, clc, close all

% Nasiha Muna and Albert E. Patterson
% Code to accompany 
% "Simple 3-D Visualization of Some Common Mathematical Minimal Surfaces using MATLAB"
% Technical Report, Alabama A&M University and University of Illinois at
% Urbana-Champaign

% Minimal Surface Plot Generation Code
% This code can be used and modified to produce 3-D plots of some common
% minimal surfaces. This code can be used, remixed, and modified as long as
% proper attribution is given to the original authors.

%% Gyroid minimal surface
% The gyroid is usually described by the following equation in terms of three
% variables: 

f_gy = @(x,y,z) sin(x).*cos(y) + sin(y).*cos(z) + sin(z).*cos(x);

% The limits for a 3-D plot should be defined here. For this case, they are
% defined as:
x_gy_l = -2*pi;
x_gy_u = 2*pi;
y_gy_l = -2*pi;
y_gy_u = 2*pi;
z_gy_l = -2*pi;
z_gy_u = 2*pi;

% The limits can be easily changed here to any values in the domain of the 
% function decribed above 

% Several factors are controllable here, but the most important for this
% visualization are the mesh density and the face color. Line color is
% assumed to be the default color. If controlling this parameter is
% desired, simply add 'edgecolor','color' to the solution command below. 

% Mesh density
MD_gy = 50; % Change as needed

% Face color
C_gy = 'yellow'; % Change as needed

% Since the function is defined in terms of the three variables, the best 
% way to plot the surface is to use the 3-D implicit function in Matlab, 
% which is set up as:  

figure    % figure command retains all of the plots instead of over-writing each one

fimplicit3(f_gy, [x_gy_l x_gy_u y_gy_l y_gy_u z_gy_l z_gy_u], 'meshdensity',MD_gy,'facecolor',C_gy)
title('Gyroid Plot')
xlabel('X')
ylabel('Y')
zlabel('Z')

%% Lidinoid minimal surface
% The lidinoid surface is usually described by the following equation in terms of three
% variables: 

f_lid = @(x,y,z) 0.5*(sin(2*x).*cos(y).*sin(z) + sin(2*y).*cos(z).*sin(x)...
    + sin(2*z).*cos(x).*sin(y)) - 0.5*(cos(2*x).*cos(2*y) + cos(2*y).*cos(2*z)...
    + cos(2*z).*cos(2*x)) + 0.15;

% The limits for a 3-D plot should be defined here. For this case, they are
% defined as:
x_lid_l = -pi;
x_lid_u = pi;
y_lid_l = -pi;
y_lid_u = pi;
z_lid_l = -0.5*pi; 
z_lid_u = 1.5*pi; % z limits are non-uniform in order to ensure that no 
% disconnected geometry is present in the plot

% The limits can be easily changed here to any values in the domain of the 
% function decribed above 

% Several factors are controllable here, but the most important for this
% visualization are the mesh density and the face color. Line color is
% assumed to be the default color. If controlling this parameter is
% desired, simply add 'edgecolor','color' to the solution command below. 

% Mesh density
MD_lid = 50; % Change as needed

% Face color
C_lid = 'blue'; % Change as needed

% Since the function is defined in terms of the three variables, the best 
% way to plot the surface is to use the 3-D implicit function in Matlab, 
% which is set up as:  
 
figure    % figure command retains all of the plots instead of over-writing each one

fimplicit3(f_lid, [x_lid_l x_lid_u y_lid_l y_lid_u z_lid_l z_lid_u], 'meshdensity',MD_lid,'facecolor',C_lid)
title('Lidinoid Plot')
xlabel('X')
ylabel('Y')
zlabel('Z')

%% Schwarz P minimal surface
% The Schwarz P surface is usually described by the following equation in terms of three
% variables: 

f_sp = @(x,y,z) cos(x) + cos(y) + cos(z);

% The limits for a 3-D plot should be defined here. For this case, they are
% defined as:
x_sp_l = -3*pi;
x_sp_u = 3*pi;
y_sp_l = -3*pi;
y_sp_u = 3*pi;
z_sp_l = -3*pi; 
z_sp_u = 3*pi;

% The limits can be easily changed here to any values in the domain of the 
% function decribed above 

% Several factors are controllable here, but the most important for this
% visualization are the mesh density and the face color. Line color is
% assumed to be the default color. If controlling this parameter is
% desired, simply add 'edgecolor','color' to the solution command below. 

% Mesh density
MD_sp = 50; % Change as needed

% Face color
C_sp = 'green'; % Change as needed

% Since the function is defined in terms of the three variables, the best 
%way to plot the surface is to use the 3-D implicit function in Matlab, 
%which is set up as:  

figure    % figure command retains all of the plots instead of over-writing each one

fimplicit3(f_sp, [x_sp_l x_sp_u y_sp_l y_sp_u z_sp_l z_sp_u], 'meshdensity',MD_sp,'facecolor',C_sp)
title('Schwarz P Plot')
xlabel('X')
ylabel('Y')
zlabel('Z')

%% Schwarz D minimal surface
% The Schwarz D surface is usually described by the following equation in terms of three
% variables: 

f_sd = @(x,y,z) sin(x).*sin(y).*sin(z) + sin(x).*cos(y).*cos(z)+cos(x).*sin(y).*cos(z)...
    + cos(x).*cos(y).*sin(z);

% The limits for a 3-D plot should be defined here. For this case, they are
% defined as:
x_sd_l = -2*pi;
x_sd_u = 2*pi;
y_sd_l = -2*pi;
y_sd_u = 2*pi;
z_sd_l = -2*pi; 
z_sd_u = 2*pi;

% The limits can be easily changed here to any values in the domain of the 
% function decribed above 

% Several factors are controllable here, but the most important for this
% visualization are the mesh density and the face color. Line color is
% assumed to be the default color. If controlling this parameter is
% desired, simply add 'edgecolor','color' to the solution command below. 

% Mesh density
MD_sd = 50; % Change as needed

% Face color
C_sd = 'red'; % Change as needed

% Since the function is defined in terms of the three variables, the best 
%way to plot the surface is to use the 3-D implicit function in Matlab, 
%which is set up as:  

figure    % figure command retains all of the plots instead of over-writing each one

fimplicit3(f_sd, [x_sd_l x_sd_u y_sd_l y_sd_u z_sd_l z_sd_u], 'meshdensity',MD_sd,'facecolor',C_sd)
title('Schwarz D Plot')
xlabel('X')
ylabel('Y')
zlabel('Z')

%% Scherk minimal surface
% The Scherk surface is usually described by the following equation in terms of three
% variables: 

f_sc = @(x,y,z) sinh(x).*sinh(y) - sin(z);

% The limits for a 3-D plot should be defined here. For this case, they are
% defined as:
x_sc_l = -pi;
x_sc_u = pi;
y_sc_l = -pi;
y_sc_u = pi;
z_sc_l = -(5/2)*pi; 
z_sc_u = (5/2)*pi;

% The limits can be easily changed here to any values in the domain of the 
% function decribed above 

% Several factors of the plot are controllable here, but the most important 
% for this visualization are the mesh density and the face color. Line color 
% is assumed to be the default color. If controlling this parameter is
% desired, simply add 'edgecolor','color' to the solution command below. 

% Mesh density
MD_sc = 50; % Change as needed

% Face color
C_sc = 'cyan'; % Change as needed

% Since the function is defined in terms of the three variables, the best 
% way to plot the surface is to use the 3-D implicit function in Matlab, 
% which is set up as:  

figure    % figure command retains all of the plots instead of over-writing each one

fimplicit3(f_sc, [x_sc_l x_sc_u y_sc_l y_sc_u z_sc_l z_sc_u], 'meshdensity',MD_sc,'facecolor',C_sc)
title('Scherk Plot')
xlabel('X')
ylabel('Y')
zlabel('Z')

%% Helicoid minimal surface
% The helicoid is one of the parametric minimal surfaces and is usually
% defined in terms of three variables and two parameters u and t such that

syms  u t
x = u.*cos(t);
y = u.*sin(t);
z = (2/3)*t;

% Since the variables are functions of the two parameters, it is only 
% necessary to place limits on the parameters. For this visualization, 
% the following were selected as limits:  
u_h_l = 0;
u_h_u = 2*pi;
t_h_l = -1;
t_h_u = 1;

% The limits can be easily changed here to any values in the domain of the 
% function decribed above 

% Several factors of the plot are controllable here, but the most important 
% for this visualization are the mesh density and the face color. Line color 
% is assumed to be the default color. If controlling this parameter is
% desired, simply add 'edgecolor','color' to the solution command below. 

% Mesh density
MD_h = 50; % Change as needed

% Face color
C_h = 'white'; % Change as needed

% Since the function is defined in terms of just the two parameters here,
% it is not necessary to use the implict solver. Therefore, the standard
% fsurf can be used. 

figure    % figure command retains all of the plots instead of over-writing each one

fsurf(x,y,z,[u_h_l u_h_u t_h_l t_h_u],'meshdensity',MD_h,'facecolor',C_h)
title('Helicoid Plot')
xlabel('X')
ylabel('Y')
zlabel('Z')

%% Catenoid minimal surface
% The catenoid is one of the parametric minimal surfaces and is usually
% defined in terms of three variables and two parameters u and v such that 

syms u v
alpha = pi/2;
x = cos(alpha).*sinh(v).*sin(u) + sin(alpha).*cosh(v).*cos(u);
y = -cos(alpha).*sinh(u)*v*cos(u)+sin(alpha)*cosh(v)*sin(u);
z = u*cos(alpha)+v*sin(alpha);

% Since the variables are functions of the two parameters, it is only 
% necessary to place limits on the parameters. For this visualization, 
% the following were selected as limits:  
u_c_l = 0;
u_c_u = 2*pi;
t_c_l = -1;
t_c_u = 1;

% The limits can be easily changed here to any values in the domain of the 
% function decribed above 

% Several factors of the plot are controllable here, but the most important 
% for this visualization are the mesh density and the face color. Line color 
% is assumed to be the default color. If controlling this parameter is
% desired, simply add 'edgecolor','color' to the solution command below. 

% Mesh density
MD_c = 50; % Change as needed

% Face color
C_c = 'magenta'; % Change as needed

% Since the function is defined in terms of just the two parameters here,
% it is not necessary to use the implict solver. Therefore, the standard
% fsurf can be used. 

figure    % figure command retains all of the plots instead of over-writing each one

fsurf(x,y,z,[u_c_l u_c_u t_c_l t_c_u],'meshdensity',MD_c,'facecolor',C_c)
title('Catenoid Plot')
xlabel('X')
ylabel('Y')
zlabel('Z')

%% Mobius strip minimal surface
% The mobius strip is one of the parametric minimal surfaces and is usually
% defined in terms of three variables and two parameters u and v such that 

syms u v
x = cos(u)*(1 + (v/2)*cos(u/2));
y = sin(u)*(1 + (v/2)*cos(u/2));
z = (v/2)*sin(u/2);

% Since the variables are functions of the two parameters, it is only 
% necessary to place limits on the parameters. For this visualization, 
% the following were selected as limits:  
u_m_l = 0;
u_m_u = 2*pi;
t_m_l = -1;
t_m_u = 1;

% The limits can be easily changed here to any values in the domain of the 
% function decribed above 

% Several factors of the plot are controllable here, but the most important 
% for this visualization are the mesh density and the face color. Line color 
% is assumed to be the default color. If controlling this parameter is
% desired, simply add 'edgecolor','color' to the solution command below. 

% Mesh density
MD_m = 50; % Change as needed

% Face color
C_m = [0.3010, 0.7450, 0.9330]; % Change as needed

% Since the function is defined in terms of just the two parameters here,
% it is not necessary to use the implict solver. Therefore, the standard
% fsurf can be used. 

figure    % figure command retains all of the plots instead of over-writing each one

fsurf(x,y,z,[u_m_l u_m_u t_m_l t_m_u],'meshdensity',MD_m,'facecolor',C_m)
title('Mobius Strip Plot')
xlabel('X')
ylabel('Y')
zlabel('Z')

















