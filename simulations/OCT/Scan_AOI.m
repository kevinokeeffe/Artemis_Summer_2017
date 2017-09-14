clear all
% German OCT sample is (from the front surface inwards): 5 nm Au, 2 nm
% Ti, 100 nm Si, 2 nm Ti, 5 nm Au, 2 nm Ti, SiO2 substrate
% treat 5 nm Au + 2 nm Ti as 7 nm Au
% angle convention: normal = 0 degrees/radians.
% Dubiously low gold refractive index values inferred from CXRO normal
% incidence reflectivity data, this means that total reflection from
% capping surface occurs at smaller angles of incidence (AOI>20^o) than is
% likely in reality.
% reconstruction quality varies somewhat with AOI

dname = fullfile(pwd, '..', '..', 'data', 'OCT');

if ~exist('smooth', 'builtin')
	smooth = @(y, w) smoothdata(y, 'movmean', w);
end

if ~exist('fittype', 'builtin')
	fittype = @(fun, varargin) fprintf('Not implemented');
end

t1 = 7E-9; % combine Au and Ti together
t2 = 100E-9;
t3 = t1;

angles = linspace(0,pi/4.5,300);
for i =1:length(angles)
    AOI = angles(i);

%
    Si_R = importdata(fullfile(dname, 'Si_mirror_35-55eV_Reflect.dat'));
    Au_R = importdata(fullfile(dname, 'Au_mirror_35-55eV_Reflect.dat'));
    Au_abs_5 = importdata(fullfile(dname, 'Au_7nm_35-55eV_Abs.dat'));
    Si_abs_140 = importdata(fullfile(dname, 'Si_100nm_35-55eV_Abs.dat')); %ignore the 140, it's for 100 nm

    %
    E = Si_abs_140(:,1); % Input Energy Spectrum 35 - 55 eV
    L = length(E);
    conv_E_omeg = 1.6E-19*2*pi/6.63E-34; %converts eV to rad/s 
    omeg = E*conv_E_omeg; % angular frequency
    d_omeg = mean(diff(omeg)); % angular frequency step size
    c = 3E8; % speed of light

    % Expected axial resolution
    E_span = E(end)-E(1); % span in energy
    E_centre = mean(E); % centre photon energy
    Res = 6.63E-34*c./(2*log(2)/pi*(E_centre)^2/E_span)/1.6E-19; % Assumes a Gaussian spectrum

    %
    Si_R = Si_R(:,2); % Reflectivity of an Si slab at normal incidence
    Au_R = Au_R(:,2); % Reflectivity of an Au slab at normal incidence


    %
    Au_n = (1-sqrt(Au_R))./(sqrt(Au_R)+1); %refractive index of Au from rearranged Fresnel coefficient (dubiously small values?)
    Au_n = Au_n+0.1; % Au refractive index seemed low (0.58 - 0.77) so bumped it up
    Si_n = (1-sqrt(Si_R))./(sqrt(Si_R)+1); % as above for Si (values more reasonable)

    %
    Spec_in = [zeros(1,10) ones(1,L-20) zeros(1,10)];
    Spec_in = smooth(Spec_in,45)'; % make smoothed top-hat input spectrum




    % Refracted angles
    alpha_1 = asin(1.*sin(AOI)./Au_n); % inside capping layer
    alpha_2 = sin(Au_n.*sin(alpha_1)./Si_n); % inside silicon
    alpha_3 = asin(Si_n.*sin(alpha_2)./Au_n); % inside buried layer
    if isreal(alpha_1) == 0 
        disp('Total reflection')
        return
    end

    % modify reflectivities for non-normal AOI using Fresnel equations
    AOI_mod1 = sqrt(1-(1./Au_n.*sin(AOI)).^2); % modifying factor
    AOI_mod2 = sqrt(1-(Si_n/Au_n*sin(alpha_1)).^2); % modifying factor
    Si_Au_R = abs((Si_n.*AOI_mod2-Au_n.*cos(alpha_1))./(Si_n.*AOI_mod2+Au_n.*cos(alpha_1))).^2; % reflectivity at interface layer from Fresnel coefficient
    Au_R = abs((1.*AOI_mod1-Au_n.*cos(AOI))./(1.*AOI_mod1+Au_n.*cos(AOI))).^2;

    % absorption needs to be increased because of longer path length
    Au_abs_5 = Au_abs_5(:,2); % 7 nm Au absorption
    mu_Au = log(Au_abs_5.^-1)/t1; % attenuation coefficient
    delta_L_Au = t1./cos(alpha_1)-t1; % extra path length (c.f. normal incidence)
    Au_abs_5 = Au_abs_5.*exp(-mu_Au.*delta_L_Au); % decreases transmission from non-normal incidence

    Si_abs_140 = Si_abs_140(:,2); % 100 nm Si absorption
    mu_Si = log(Si_abs_140.^-1)/t2; % attenuation coefficient of silicon
    delta_L_Si = t2./cos(alpha_2)-t2;
    Si_abs_140 = Si_abs_140.*exp(-mu_Si.*delta_L_Si);

    %
    % need to be adjusted for differeny AOI
    R1 = Au_R; % Front surface capping layer
    R2 = Si_Au_R.*(1-Au_R).^2.*Au_abs_5.^2; % Reflection from capping Au-Si interface
    R3 = Si_Au_R.*(1-Au_R).^2.*Au_abs_5.^2.*(1-Si_Au_R).^2.*Si_abs_140.^2; % Reflection from Si-rear buried Au interface
    R4 = Si_Au_R.*(1-Au_R).^2.*Au_abs_5.^2.*(1-Si_Au_R).^2.*Si_abs_140.^2.*(1-Au_R).^2.*Au_abs_5.^2; % reflection from buried Au Si substrate interface 

    % Cosine arguments
    Arg_1 = (omeg*2/c).*Au_n*t1.*cos(alpha_1);
    Arg_2 = (omeg*2/c).*Si_n*t2.*cos(alpha_2);
    Arg_3 = (omeg*2/c).*Au_n*t1.*cos(alpha_3);

    %
    Term_1 = 2*sqrt(R1.*R2).*cos(Arg_1); % i = 2, j = 1 
    Term_2 = 2*sqrt(R2.*R3).*cos(Arg_2)+2*sqrt(R1.*R3).*cos(Arg_1+Arg_2); %i=3, j=1,2
    Term_3 = 2*sqrt(R4.*R3).*cos(Arg_3)+2*sqrt(R4.*R2).*cos(Arg_2+Arg_3)...
        +2*sqrt(R4.*R1).*cos(Arg_1+Arg_2+Arg_3);% i = 4, j=1,2,3
    %
    Spec_out = Spec_in'.*(R1+R2+R3+R4+Term_1+Term_2+Term_3);

%     figure(1)
%     plot(E,Spec_in*max(Spec_out),':')
%     hold on
%     plot(E,Spec_out)
%     hold off
%     xlabel('Photon Energy (eV)')
%     ylabel('Intensity (arb.)')
%     set(gca,'linewidth',4)
%     legend('Incident spectrum','Reflected spectrum','location','best')
%     title('Spectrum In - Spectrum Out')
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
    
    if g.rsquare<0.99
        clear f
        f = @(x) 9.44.*exp(-x./7.85);
        Recon_pad = abs(fft(Spec_out_pad'./feval(f,omeg_pad/conv_E_omeg)'));
    else
        Recon_pad = abs(fft(Spec_out_pad'./feval(f,omeg_pad/conv_E_omeg)));
    end
    %Recon_pad = abs(fft(Spec_out_pad));
    %
    Delay_pad = [1:L2]*2*pi/(L2*d_omeg); % padded delay axis
    Depth = mean(Si_n)*1E9*Delay_pad*3E8/2;% Fudge factor 2 in here. Also not the correct way to correct for dispersion.
    [~,b1] = min(abs(Depth*1E-9-(t1+t2))); % find where buried layer begins
    [~,b2] = min(abs(Depth*1E-9-(t1+t2+t3))); % find where buried layer ends
    [~,b3] = min(abs(Depth*1E-9-(t1))); % find where capping layer ends
    
    R(i,:) = Recon_pad/max(Recon_pad);

end
%
degrees = 90*angles/pi;
imagesc(Depth(1:200),degrees,R(:,1:200))
view(2)
axis tight
xlabel('Depth (nm)')
ylabel('AOI (degrees)')
line([t1+t2 t1+t2]*1E9,[min(degrees) max(degrees)],[100 100],'color','k','linestyle',':','linewidth',2)
line([t1+t2+t3 t1+t2+t3]*1E9,[min(degrees) max(degrees)],[100 100],'color','k','linestyle',':','linewidth',2)
title('OCT Reconstruction as a function of AOI')
% figure(2)
% plot(Depth, Recon_pad) 
% hold on
% plot(ones(1,10)*Depth(b1),linspace(0, max(Recon_pad),10),'k:','linewidth',3)
% plot(ones(1,10)*Depth(b2),linspace(0, max(Recon_pad),10),'k:','linewidth',3)
% plot(ones(1,10)*0,linspace(0, max(Recon_pad),10),'k:','linewidth',3)
% plot(ones(1,10)*Depth(b3),linspace(0, max(Recon_pad),10),'-.','color',0.7*ones(1,3),'linewidth',3)
% hold off
% xlim([-1 200])
% xlabel('Depth (nm)')
% ylabel('|F[S(\omega)]|')
% set(gca,'linewidth',4)
% title('Reconstructed Depth Profile')
% legend('Reconstruction','Buried layer location','Capping layer')