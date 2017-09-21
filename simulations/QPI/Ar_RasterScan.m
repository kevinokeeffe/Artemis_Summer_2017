%	Scientific / Atomic units' constants
[SI, AU] = Units;
C = SI.C*1e3*1e-15;			%	Speed of light [mm/fs]

%	Radial axes
%	----------------------------------------------------------------------------
W0 = .03;					%	Beam waist (e^-2 intensity half width) [mm]
Nr = 2^8;					%	Radial samples
rmax = 5*W0;				%	Radial boundary [mm]

%	Generate Hankel transform matrix & associated items
if ~exist('H', 'var') || ~isequal(H.Nr, Nr)
	H = hankel_matrix(0, rmax, Nr);
end

%	Create radial intensity distribution (unity peak intensity)
Er = exp(-2*(H.r./W0).^2);

%	Radial data subset used for HHG SFA simulations - harmonics only generated
%	near centre of beam
indrq = (H.r<=3*W0);
Nrq = sum(indrq);

% plot(H.r, Er);



%	Spectral axes
%	----------------------------------------------------------------------------
l0 = 780;									%	Central wavelength [nm]
zR = pi*W0.^2 ./ (l0*1e-6);					%	Rayleigh length [mm]
w0 = chng_rng(l0);							%	Central frequency [rad/fs]
T = 2*pi/w0;								%	Period [fs]
T_au = T * (1e-15/AU.time);					%	Period [au]

Nwq = 2^13;
qmax = 75;									%	Max harmonic order
wmax = qmax * w0;							%	Max frequency [rad/fs]
dw = wmax/(Nwq-1);							%	Frequency resolution [rad/fs]
dt = 2*pi/(Nwq*dw);							%	Temporal resolution [fs]

wq = (0:Nwq-1)'*dw;							%	Frequency axis [rad/fs]
l = chng_rng(wq);							%	Wavelength axis [nm]
t = ((0:Nwq-1).'-floor(Nwq/2))*dt;			%	Temporal axis [fs]
t_au = t*1e-15/AU.time;						%	Temporal axis [au]
wq_eV = wq * 1e15*AU.h_bar/AU.e;				%	Frequency axis [eV]
q = wq/w0;									%	Harmonic order axis
k = wq/C;									%	Wave vector [rad/mm]
wq2 = wq.^2;

%	Frequency subset for fundamental
%	----------------------------------------------------------------------------

DT = 30;
DW = 4*log(2)/DT;
Et = .5*gauss1D(t, 0, sqrt(2)*DT, 1).*exp(-1i*w0*t);
Ew = Time2Freq(Et);

indw1 = (abs(wq-w0)<4*DW);
w1 = wq(indw1);
k1 = k(indw1);
Nw1 = length(w1);

%	Longitudinal wavevector for harmonics/fundamental
kz = sqrt(k.'.^2 - H.kr.^2) - k.'; kz(k.'<H.kr) = 0;
kz1 = sqrt(k1.'.^2 - H.kr.^2) - k1.'; kz1(k1.'<H.kr) = 0;

%% Propagate fundamental
%  =============================================================================
Nz = 100;
z = linspace(-1, 1, Nz).'*zR*.75;

%	Calculate field
E1 = (Er./H.JR).*(Ew(indw1).');

%	Propagate field in k-space
E1 = (H.T*E1).*exp(1i*kz1.*reshape(z, 1, 1, Nz));

%	Convert back to real-space
E1 = reshape(H.T*E1(:, :), Nr, Nw1, Nz).*H.JR;

%% HHG SFA simulations
%  =============================================================================

%	Ionisation potential
Ip_eV =	15.759;								%	Argon [eV]
Ip_au = Ip_eV * AU.e / AU.En;				%	[au]

%	Peak intensity / *real* field strength 
I0p = 1e14;								%	Peak intensity [W/cm^2]
e0 = sqrt(2e4*I0p / SI.eps_0 / SI.C);		%	Peak real field strength [V/m]
e0_au = e0 / (AU.force / AU.e);				%	Peak real field strength [au]

%	Harmonic spectrum at peak intensity
Eq0 = Time2Freq(QSFA1(t_au, e0_au*real(Et), Ip_au, T_au)) .* wq2;
Iq0 = abs(Eq0).^2;

Eq1 = Time2Freq(QSFA(t_au, e0_au*real(Et), Ip_au, T_au)) .* wq2;
Iq1 = abs(Eq1).^2;
% Eq0 = Freq2Time(Time2Freq(QSFA1(t_au, e0_au*Et, Ip_au, T_au)) ...
% 	.* [wq2; flipud(wq2)]);
% Iq0 = abs(Time2Freq(Eq0)).^2;
% Iq0 = Iq0(indwq);
semilogy(q, [Iq0 Iq1]);
return

%	Simple check to see if simulation results already exist and of correct size
if exist('Eq.mat', 'file')
	load('Eq.mat');
	if ~isequal(size(Eq), [Nrq, Nw, Nz])
		clear Eq;
	end
end

%	Run simuulation if either simulation doesn't exist, or of incorrect size
if ~exist('Eq', 'var')

	fprintf('\tSimulating HHG ...\n');
	E = zeros(Nw, 1, 'like', 1i);
	Eq = zeros(Nrq, Nw, Nz, 'like', 1i);
	nz = 1;
	multiWaitbar('Longitudinal Position', 'Reset');
	for nz=1:Nz
		multiWaitbar('Radial Position', 'Reset');
		for nr=1:Nrq
			E(indw1) = Freq2Time(e0_au*E1(nr, :, nz));
			Eq(nr, :, nz) = Row(Time2Freq(QSFA1(t_au, E, Ip_au, T_au)) .*w2);
			multiWaitbar('Radial Position', 'Increment', 1/Nrq);
		end
		multiWaitbar('Longitudinal Position', 'Increment', 1/Nz);
	end
	multiWaitbar('CloseAll');
	
	save('Eq.mat', 'Eq');
end

return
%% Propagate HHG
%  =============================================================================

%	Fill in gaps - set HHG intensity outside central region to zero
if size(Eq, 1)~=H.Nr
	Eq = [Eq; zeros(H.Nr-Nrq, Nwq, Nz, 'like', 1i)];
end

%	Convert to reciprocal space
Eq1 = reshape(H.T*(Eq(:, :)./H.JR), H.Nr, Nwq, Nz) .* H.JV;

%	Back-propagate to same longitudinal plane (z=0: location of focus)
Eq1 = Eq1 .* exp(-1i*kzq.*reshape(z, 1, 1, Nz));

%% Calculate harmonic locations
%  =============================================================================

%	Number of harmonics
Nq = 17;

%	Find frequency index of harmonic peak
I = mean(mean(abs(Eq1).^2, 3), 1);
nwq = sort(FindPeaks(I, 10, [], Nq));
semilogy(qq, I, qq(nwq), I(nwq), '+');

%% Calculate interferograms (raster scan)
%  =============================================================================
I = @(nth, nw) abs(Column(Eq1(nth, nw, :)) + Row(Eq1(nth, nw, :))).^2;
% I = squeeze(Eq1(1, nwq, :)).';
% I = abs(reshape(I, Nz, 1, Nq) + reshape(I, 1, Nz, Nq)).^2;

%	Plot on-axis interferograms (two different harmonics)
clf
subplot(1,2,1)
imagesc(z, z, I(1, nwq(8)));

subplot(1,2,2)
imagesc(z, z, I(1, nwq(12)));
return
%%


% clf
% subplot(2,1,1);
% h = pcolor(qq, H.r, log10(abs(Eq).^2));
% h.LineStyle = 'None';
% 
% subplot(2,1,2);
% h = pcolor(qq, H.kr, log10(abs((H.T*(Eq./H.JR)).*H.JV).^2));
% h.LineStyle = 'None';
% 
% % semilogy(q, Iq);
% 

