clear all
% NORMAL INCIDENT FOR SIMPLICITY
% Current OCT sample is (from the front surface inwards): 5 nm Au, 2 nm
% Ti, 100 nm Si, 2 nm Ti, 5 nm Au, 2 nm Ti, SiO2 substrate
% method follows Fuchs, S. et al. Nanometer resolution optical coherence 
% tomography using broad bandwidth XUV and so?t x-ray radiation. Sci. Rep. 6, 
% 20658; doi: 10.1038/srep20658 (2016) supplementary material

if ~exist('smooth', 'builtin')
	smooth = @(y, w) smoothdata(y, 'movmean', w);
end

if ~exist('fittype', 'builtin')
	fittype = @(fun, varargin) fprintf('Not implemented');
end

t1 = 7E-9; % combine Au and Ti layers together for simplicity
t2 = 100E-9;% silicon separating layer
t3 = t1; % buried layer is same composition at capping layer

% Import data
dname = fullfile(pwd, '..', '..', 'data', 'OCT');
Si_R = importdata(fullfile(dname, 'Si_mirror_35-55eV_Reflect.dat'));	% raw Si reflectivity at normal incidence
Au_R = importdata(fullfile(dname, 'Au_mirror_35-55eV_Reflect.dat'));	% raw Au reflectivity at normal incidence
Au_abs_5 = importdata(fullfile(dname, 'Au_7nm_35-55eV_Abs.dat'));		% 7 nm thick Au transmission
Si_abs_140 = importdata(fullfile(dname, 'Si_100nm_35-55eV_Abs.dat'));	% 100 nm thick Si transmission


%
E = Si_abs_140(:,1); % Input Energy Spectrum 35 - 55 eV
L = length(E);
conv_E_omeg = 1.6E-19*2*pi/6.63E-34; %converts eV to seconds 
omeg = E*conv_E_omeg; % angular frequency
d_omeg = mean(diff(omeg)); % angular frequency step size
c = 3E8; % speed of light

% Expected axial resolution
E_span = E(end)-E(1); % span in energy
E_centre = mean(E); % centre photon energy
Res = 6.63E-34*c./(2*log(2)/pi*(E_centre)^2/E_span)/1.6E-19; % Assumes a Gaussian spectrum

%
Si_R = Si_R(:,2); % Reflectivity of an Si slab
Au_R = Au_R(:,2); % Reflectivity of an Au slab
Si_abs_140 = Si_abs_140(:,2); % 100 nm Si absorption
Au_abs_5 = Au_abs_5(:,2); % 7 nm Au absorption

%
Au_n = (1-sqrt(Au_R))./(sqrt(Au_R)+1); %refractive index of Au from rearranged Fresnel coefficient (dubiously small values?)
Si_n = (1-sqrt(Si_R))./(sqrt(Si_R)+1); % as above for Si (values more reasonable)

%
Si_Au_R = abs((Si_n-Au_n)./(Si_n+Au_n)).^2; % reflectivity at interface layer from Fresnel coefficient

%
Spec_in = [zeros(1,10) ones(1,L-20) zeros(1,10)];
Spec_in = smooth(Spec_in,45)'; % smooth input spectrum to make it a little more realistic

%
R1 = Au_R; % Front surface capping layer reflectivity
R2 = Si_Au_R.*(1-Au_R).^2.*Au_abs_5.^2; % Reflection from 1st capping Au-Si interface
R3 = Si_Au_R.*(1-Au_R).^2.*Au_abs_5.^2.*(1-Si_Au_R).^2.*Si_abs_140.^2; % Reflection from Si-rear buried Au interface
R4 = Si_Au_R.*(1-Au_R).^2.*Au_abs_5.^2.*(1-Si_Au_R).^2.*Si_abs_140.^2.*(1-Au_R).^2.*Au_abs_5.^2; % reflection from buried Au Si substrate interface 

% Cosine arguments
Arg_1 = (omeg*2/c).*Au_n*t1;
Arg_2 = (omeg*2/c).*Si_n*t2;
Arg_3 = (omeg*2/c).*Au_n*t1;

%
Term_1 = 2*sqrt(R1.*R2).*cos(Arg_1); % i = 2, j = 1 
Term_2 = 2*sqrt(R2.*R3).*cos(Arg_2)+2*sqrt(R1.*R3).*cos(Arg_1+Arg_2); %i=3, j=1,2
Term_3 = 2*sqrt(R4.*R3).*cos(Arg_3)+2*sqrt(R4.*R2).*cos(Arg_2+Arg_3)...
    +2*sqrt(R4.*R1).*cos(Arg_1+Arg_2+Arg_3);% i = 4, j=1,2,3
%
Spec_out = Spec_in'.*(R1+R2+R3+R4+Term_1+Term_2+Term_3);
figure(1)
plot(E,Spec_in*max(Spec_out),':')
hold on
plot(E,Spec_out)
hold off
xlabel('Photon Energy (eV)')
ylabel('Intensity (arb.)')
set(gca,'linewidth',4)
legend('Incident spectrum','Reflected spectrum','location','best')
title('Spectrum In - Spectrum Out')
%
% Reconstruction with padding and reflectivity
% trend correction.
F_type_exp = fittype('a.*exp(-x./b)','coefficients',{'a','b'});  
[a,b] = findpeaks(Spec_out,'npeaks',4); % find peaks in reflected spectrum

[f,g] = fit(E(b),Spec_out(b),F_type_exp,'Startpoint',[2 12]); % fit a decay to approximate reflectivity

omeg_pad =0:d_omeg:350*conv_E_omeg; % extend spectral range to 0-350 eV
L2 = length(omeg_pad);
[~,b_1] = min(abs(omeg_pad-E(1)*conv_E_omeg));
[~,b_end] = min(abs(omeg_pad-E(end)*conv_E_omeg));
Spec_out_pad = [zeros(1,b_1-1) Spec_out' zeros(1,length(omeg_pad)-b_end)];

Recon_pad = abs(fft(Spec_out_pad'./feval(f,omeg_pad/conv_E_omeg)));
%Recon_pad = abs(fft(Spec_out_pad));
%
Delay_pad = [1:L2]*2*pi/(L2*d_omeg); % padded delay axis
Depth = mean(Si_n)*1E9*Delay_pad*3E8/2;% Fudge factor 2 in here. Also not the correct way to correct for dispersion.
[~,b1] = min(abs(Depth*1E-9-(t1+t2))); % find where buried layer begins
[~,b2] = min(abs(Depth*1E-9-(t1+t2+t3))); % find where buried layer ends
figure(2)
plot(Depth, Recon_pad) 
hold on
plot(ones(1,10)*Depth(b1),linspace(0, max(Recon_pad),10),'k:','linewidth',3)
plot(ones(1,10)*Depth(b2),linspace(0, max(Recon_pad),10),'k:','linewidth',3)
hold off
xlim([0 200])
xlabel('Depth (nm)')
ylabel('|F[S(\omega)]|')
set(gca,'linewidth',4)
title('Reconstructed Depth Profile')
legend('Reconstruction','Buried layer location')