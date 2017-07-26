% Two target calculation using analytic HHG model
% sense of z matters for CEP phase, whereas dipole phase is symmetric wrt z
% "visibility" best when moving target through focal plane, encountering low
% values of dI/dz;
% Normalised longitudinal coordinates (z/zr) used throughout.
% see equation (5) in Catoire et al. PRA 94,063401 (2016) for info on
% model.
% alpha = 4E-14 corresponds to short trajectory.
% First part shows two target fringe patterns for three different fixed
% target posistions.
% Second part shows modualation "k-vector" linear dependence on q and how
% this differs from Gouy only model (as per Laban PRL 2012 paper).
% Third part shows the plot derived for a single order as a function of all
% possible combinations of two cell locations.
% Fourth part just shows a typical modulation pattern.

%%
zr = 3; % Rayleigh range (mm)
I_0 = 2E14; % peak intensity (W/cm^2)
q = 23; % harmonic order
q_eff = 2.5; % effective nonlinearity
alpha = 4E-14; % alpha coefficient (cm^2/W)

z = 2.1*linspace(-zr,zr,1000); % span of z locations
Iz = 1./(1+(z./zr).^2); % Fundamental intensity variation with z
CEP = -atan(z/zr); % Gouy phase - sign important

% three choices of static target location
[~,pix_targ(1)] =min(abs(z-0)); % at the focal plane
[~,pix_targ(2)] =min(abs(z-zr/3)); % z_r/3 toward laser
[~,pix_targ(3)] =min(abs(z+zr/3)); % z_r/3 away from laser

% some pretty colors
c{1} = [0.1 0.3 0.4]*2;
c{2} = [0.4 0.3 0.2]*2;
c{3} = [0.1 0.3 0.2]*2;

% Harmonic field according to Catoire et al.
A_q = Iz.^q_eff.*exp(-1i*(q*CEP+alpha*I_0*Iz));
% calculate interference pattern for three different static target lcoations
figure(1)
for i = 1:3
    A_s = A_q(pix_targ(i));

    I_total = abs(A_s+A_q).^2;
    
    plot(z,I_total,'color',c{i})
    hold on
end
% plot locations of static target with dotted lines
plot(ones(1,20)*z(pix_targ(1)),linspace(0,4,20),':','color',c{1})
plot(ones(1,20)*z(pix_targ(2)),linspace(0,4,20),':','color',c{2})
plot(ones(1,20)*z(pix_targ(3)),linspace(0,4,20),':','color',c{3})
hold off
legend('z_1 = 0','z_1 = z_r/3','z_1 = -z_r/3','location','northwest')
axis tight
set(gca,'xtick',[-2*zr -zr 0 zr 2*zr],'xticklabel',{'-2z_r', '-z_r', '0', 'z_r', '2z_r'},'linewidth',4)
title(['\alpha\timesI_0 = ' num2str(alpha*I_0) '; q = ' num2str(q) '; --> toward spectrometer.'])
ylabel('Harmonic Intensity (arb.)')
xlabel('Moving Gas Cell Location')
%% Same model as above, just varying q.
% Fix target at focal plane, find first maximia location as a fn. of q
clear all
zr = 3; % Rayleigh range (mm)
I_0 = 2E14; % peak intensity (W/cm^2)
q_eff = 2.5; % effective nonlinearity
alpha = 4E-14; % alpha coefficient (cm^2/W)

z = linspace(0,zr,1000); % span of z locations
Iz = 1./(1+(z./zr).^2); % Fundamental intensity variation with z
CEP = -atan(z/zr); % Gouy phase - sign important

% three choices of static target location
%[~,pix_targ] =min(abs(z-0)); % at the focal plane
%[~,pix_targ] =min(abs(z-zr/3)); % z_r/3 toward laser
[~,pix_targ] =min(abs(z+zr/3)); % z_r/3 away from laser

% Harmonic field according to Catoire et al.
q_all = 17:2:29;
for i=1:length(q_all)
    q = q_all(i);
    A_q = Iz.^q_eff.*exp(-1i*(q*CEP+alpha*I_0*Iz));
    A_s = A_q(pix_targ);
    I_total = abs(A_s+A_q).^2;
    [a,b] = findpeaks(I_total,'Npeaks',1);
    mod_x(i) = z(b);
end
% plot "k-vector" as defined in Kevin's document (15/03/2017)
figure(2)
plot(q_all,2*pi./mod_x,'o')
line_fit = polyfit(q_all,2*pi./mod_x,1);
hold on
plot(q_all,polyval(line_fit,q_all))
plot(q_all,q_all/zr)
hold off
legend('Calculated values','Linear fit','q/z_r','location','northwest')
xlabel('q')
ylabel('"k-vector" (mm^{-1})')
set(gca,'xtick',q_all,'linewidth',4)
%%
%set(gcf, 'PaperPositionMode', 'auto');
%print -depsc2 Two-source-fringes.eps
%close
%% plot as function of z_1 and z_2
clear all
zr = 3; % Rayleigh range (mm)
I_0 = 2E14; % peak intensity (W/cm^2)
q = 23; % harmonic order
q_eff = 2.5; % effective nonlinearity
alpha = 4E-14; % alpha coefficient (cm^2/W)

z = 2.1*linspace(-zr,zr,1000); % span of z locations
Iz = 1./(1+(z./zr).^2); % Fundamental intensity variation with z
CEP = -atan(z/zr); % Gouy phase - sign important
A_q = Iz.^q_eff.*exp(-1i*(q*CEP+alpha*I_0*Iz));
for i = 1:length(z)
    for j = 1:length(z)
          A_q = Iz(i).^q_eff.*exp(-1i*(q*CEP(i)+alpha*I_0*Iz(i)));
          A_s = Iz(j).^q_eff.*exp(-1i*(q*CEP(j)+alpha*I_0*Iz(j)));
          I_total(i,j) = abs(A_s+A_q).^2;
    end
end
figure(3)
mesh(z/zr,z/zr,I_total)
view(2)
axis tight
daspect([1 1 1])
xlabel('z_1/z_r')
ylabel('z_2/z_r')
title(['q=' num2str(q) '; All combinations of {z_1,z_2}'])
%%
figure(4)
plot(z/zr,I_total(:,500))
axis tight
xlabel('z_2/z_r')
title('Typical Interference Pattern') 