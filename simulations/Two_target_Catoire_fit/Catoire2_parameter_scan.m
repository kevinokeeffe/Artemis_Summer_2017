clear all
% compare experimental raster scans to those produced by 2 target, two
% trajectory model after calculating and comparing an error metric for all
% permutations of parameter values defined by the user
% square rasters only (equal number of T1 and T2 z-locations)
% import experimental data, see below for naming convention
A = importdata('2017-08-22_170826_Sort_focus_235uJ_raster_data.mat');

% define constants of experiment
Constants.order = 15; % order of raster scans analysed
Constants.wavelength = 0.78; % wavelength
Constants.Max_z = max(A.zT1); % max z position in raster scan
Constants.R = 0; % radial coordinate of raster scan
Constants.z_det = 0.69E6; %source-detector distance

% Variable names
Var_names = {'w0', 'del_phi', 'A_long', 'qeff', 'alpha_short', 'alpha_long'...
   , 'z_off', 'z_space'};

% if adding a new variable also add name to Var_names and additional index
% later to main calculation (not ideal)
% input variable ranges here.
Variables.w0 = 47:48:49; % spot size (microns)
Variables.del_phi = 0.3*pi; % relative phase shift of trajectories
Variables.A_long = 0.3; % relative amplitude of long trajectory
Variables.qeff = 6; % effective nonlinearity
Variables.alpha_short = 1; %\alpha_{short} x peak intensity
Variables.alpha_long = 18; %\alpha_{long} x peak intensity
Variables.z_off = 0*1E3; % distance between focal plane and T2
Variables.z_space = 1.2*1E3; % distance between T1 and T2 at "zero" separation

% get the number of parameter values for each variable
for ii = 1:length(Var_names)
    index_max(ii) = length(eval((['Variables.' Var_names{ii}])));
end

% total number of permutations (16,000 permutations takes ~ 3 hours to
% compute)
number_instances = prod(index_max);

% coordinates
dz  = mean(diff(A.zT1)); % longitudinal step size from experiment
z = 1E3*[(-max(A.zT1)-2):dz:(max(A.zT1)+2)]; % make the longitudinal coordinate
% assume step size was constant in experiment.

L_r = 0.5E3; % number of radial points
r = linspace(0,1E3,L_r); % radial coordinate
dr = mean(diff(r)); % radial step size step size


% import experimental raster scans and radial coordinates
% this part can be generalised a bit better

R_ras(:,1) = 0;
R_ras(:,2) = A.r_midaxis;
R_ras(:,3) = A.r_offaxis;

expt(:,:,1) = A.I_onaxis/max(max(A.I_onaxis));
expt(:,:,2) = A.I_midaxis/max(max(A.I_midaxis));
expt(:,:,3) = A.I_offaxis/max(max(A.I_offaxis));


% pre-allocate grid of possible indexs for subsequent arrayfun operation
[rd1,rd2,rd3,rd4,rd5,rd6,rd7,rd8] = ndgrid(1:index_max(1),1:index_max(2),...
    1:index_max(3),1:index_max(4),1:index_max(5),1:index_max(6),...
    1:index_max(7),1:index_max(8));

%
tic

phs_curvIO = 0; % Switch waveefront curvature on of off. 1 = one, 0 = off


% main calculation, run function Model_instance_n_plot on all possible permutations of parameter values 
allsums = arrayfun(@(r1,r2,r3,r4,r5,r6,r7,r8) ...
Model_instance_n_plot(expt,Constants,Variables,[r1,r2,r3,r4,r5,r6,r7,r8], ...
R_ras, r,z,length(A.zT1),phs_curvIO,'calc'),...
rd1, rd2, rd3, rd4, rd5, rd6, rd7, rd8, 'UniformOutput', false); % 
% Model_instance_n for comparing multiple raster scan radial positions
    
t=toc

% reshape output of main calculation back into correct form
chi_sq = reshape(vertcat(allsums{:}),index_max);

% extract indexes of minimum chi squared value
[min_chi_sq,min_ind] = min(chi_sq(:));
min_chi_val = chi_sq(min_ind);
% probably an easier way of coding this...
[chi_ind(1),chi_ind(2),chi_ind(3),chi_ind(4), chi_ind(5), chi_ind(6), ...
    chi_ind(7), chi_ind(8)] = ind2sub(size(chi_sq),min_ind);
%

% plot experimental raster scans and best fit scans from model
set(0,'defaultaxesfontsize',32)
Model_instance_n_plot(expt,Constants, Variables, chi_ind, ...
R_ras, r,z,length(A.zT1),phs_curvIO,'plot');


% add best fit values to array B
for ii = 1:length(Var_names)
    C = eval(['Variables.' Var_names{ii}]);
    B(ii) = C(chi_ind(ii));
end
