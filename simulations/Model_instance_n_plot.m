% Model_instance(expt. raster scan, expt. constants, expt_variables, indexes to evaluate
% off axis position, radial coordinate, z coordinate, expt. z length,
% driver curvature on/off)

% if plot_model_instance.m is changed, make the same changes here

function  yy = Model_instance_n_plot(expt, Constants,Variables, arr_ind, R_ras, r, z, L_expt, Phs_on, mode)

w = Variables.w0(arr_ind(1)); % IR driver waist
lam_0 = Constants.wavelength; % driver wavelength 
q = Constants.order; % order
lam_q = lam_0/q; % harmonic wavelength
k_vec = 2*pi/lam_q; % harmonic wavevector
k0 = 2*pi/lam_0; % driver wavevector
zr = pi*w^2/lam_0; % Rayleigh range

% 
del_phi_ls = Variables.del_phi(arr_ind(2)); % relative phase of trajectories
Long_amp = Variables.A_long(arr_ind(3)); % relative strength of long trajectory
% 
qeff = Variables.qeff(arr_ind(4)); % effective nonlinearity
alpha_short = Variables.alpha_short(arr_ind(5)); % long trajectory alpha (10^-14 cm^2/W)
alpha_long = Variables.alpha_long(arr_ind(6)); % short trajectory alpha
I_max = 1; % peak intensity  10^14 W/cm^-2
% 
z_off = Variables.z_off(arr_ind(7)); %  zero gas cell separation offset from focal plane
z_off = z_off+1E-13; % NaN killer
% 
% % longitudinal coordinate and quantities
L_z = length(z); % number of z locations
L_r = length(r);
% 
CEP = -atan((z-z_off)/zr); % Gouy phase

% extend 1-D coordinates to 2-D
rr = kron(r,ones(L_z,1));
zz = kron(z,ones(L_r,1));

% step sizes
dr = mean(diff(r));
dz = mean(diff(z));

R_z = (zz-z_off).*(1+(zr./(zz-z_off)).^2); % Radius of curvature (Gaussian beam)

if Phs_on == 1
    Phs_curv = k0*rr.^2./(2*R_z'); % Driver phase front
elseif Phs_on == 0
    Phs_curv = 0; % plane wave driver
else
    disp('Wrong phase activation command')
end

I_z = 1./(1+((zz-z_off)./zr).^2); % Fundamental intensity variation with z
% 
EIR = exp(-rr.^2/(2*w^2)).*sqrt(I_z'); % IR driver transverse FIELD amplitude
% % denominator factor 2 added from consistency with Rayleigh range def.


% Extend to 2-D. Gouy phase has no transverse dependence in this model
CEP = kron(CEP,ones(L_r,1));
CEP = CEP';


Z = Constants.z_det; % distance to propagate after second cell 0.6E6

R = 0:2*pi/(L_r*dr):2*pi/dr-2*pi/(L_r*dr); % far-field transverse plane
R = Z*R/k_vec; % scaling factor the same as a Fourier transform

 
% pre-allocate
Exit_field = zeros(L_z,L_z,L_r);
allsums{L_z,L_z} = [];

% Harmonic field as per simple Catoire et al model
E_HHG = @(alpha,CEP,del_phi_ls)...
    (abs(EIR)).^qeff.*exp(1i*(del_phi_ls+alpha*I_max.*abs(EIR).^2+q*(CEP+Phs_curv)));
% % now with driver curvature (Phs_curv) imparted to the harmonic field. %
% alpha plus or minus?
 
% add long and short contributions together (coherently)
m = Long_amp*E_HHG(alpha_long,CEP,del_phi_ls)+E_HHG(alpha_short,CEP,0);

[ridx1, ridx2] = ndgrid(1:size(m, 1)); % grids for arrayfun 
 
% add together harmonic field for all combinations of T1 and T2 positions
allsums = arrayfun(@(r1, r2) m(r1, :) + m(r2, :), ridx1, ridx2, 'UniformOutput', false);
Exit_field = reshape(vertcat(allsums{:}),L_z,L_z,L_r);

% % Propagation using Hankel transform: extended Marcel Leutenegger's
% % function from the file exchange to accept L_z^2 x L_r size functions
h = Exit_field(1:L_z,1:L_z,:);
h = reshape(h,size(h,2).^2,L_r);
h = h';
 
r_temp=0:L_r-1;
[r_n,~]=sort(r_temp(:).');
 
k=pi/L_r*(0:L_r-1);
 
r_temp=[(r_temp(2:end) + r_temp(1:end-1))/2 r_temp(end)];
I=2*pi./k(:)*r_temp.*besselj(1,k(:)*r_n); % Propagation "Kernal"
I(k == 0,:)=pi*r_temp.*r_temp;
I=I - [zeros(numel(k),1) I(:,1:end-1)];
 
RS_Far_field = I*h; % actual propagation step

RS_Far_field = reshape(RS_Far_field',L_z,L_z,L_r);
 
ofs = 1E-3*(Variables.z_space(arr_ind(8))); % convert to mm for convenience

max_z_expt = (L_expt-1)*dz*1E-3; % (L_expt-1)*dz = 8 = max z in expt

% indexes to select correct sub-set from model raster scan
[~,b0] = min(abs(1E-3*z));
[~,b0ofs] = min(abs(1E-3*z-ofs));
[~,b8] = min(abs(1E-3*z-max_z_expt-ofs)); % 
[~,b_8] = min(abs(1E-3*z+max_z_expt));

% discritizing bug fix
if length(b0ofs:b8) == L_expt-1
    b8 = b8+1;
end

% either perform the calculation of the total chi squared value for n
% raster scans considered, or plot the n raster scans

if strcmp(mode,'calc')

% calculate chi square values
    for ii = 1:size(R_ras,2)

        [~,b_r] = min(abs(R-R_ras(:,ii)));
        subset = abs(RS_Far_field(b_8:b0,b0ofs:b8,b_r)).^2; 
        model = flipud(subset/max(max(subset)));

        chi_sq(ii) = sum(sum(abs(model-[expt(:,:,ii)]').^2));

    end

% sum over n chi square values (n = number of raster scans)
yy = sum(chi_sq);

elseif strcmp(mode,'plot')
    for ii = 1:size(R_ras,2)
        [~,b_r] = min(abs(R-R_ras(ii))); % R_ras in microns; finds location off-axis

        figure
        subplot(2,1,1)
        imagesc(1E-3*z(b0ofs:b8),1E-3*z(b_8:b0),abs(RS_Far_field(b_8:b0,b0ofs:b8,b_r)).^2)
        view(2)
        axis tight
        daspect([1 1 1])
        title(['Model: R = ' num2str(round(1E-3*R_ras(ii),2)) 'mm'])
        subplot(2,1,2)
        z_expt = [0:dz*1E-3:max_z_expt];
        imagesc(z_expt,-z_expt,expt(:,:,ii)') %FIX
        view(2)
        axis tight
        daspect([1 1 1])
        title(['Expt: R = ' num2str(round(1E-3*R_ras(ii),2)) 'mm'])
    end
else
    disp('Wrong mode (either "calc" or "plot")')
end