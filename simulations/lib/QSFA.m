%% QUICK STRONG FIELD APPROXIMATION (QSFA)
% Simulates high harmonic generation (HHG) using a quick SFA model:
%	The action is calculated using the stationary momentum, but only 
%
%	Usage:
%		xt = QSFA(t, Et, Ip, tau)
%
%		xt = Time-dependent atomic dipole moment
%		t  = Temporal axis [au]
%		Et = Input real electric field [au]
%		Ip = Ionization potential [au]
%		tau_max = Maximum excursion time in continuum [au]
%%
function xt = QSFA(t, Et, Ip, tau_max, dipole, fft_int)

if ~exist('dipole', 'var')
	dipole = @Hydrogenlike;
elseif ischar(dipole)
	dipole = eval(dipole);
end

if ~exist('fft_int', 'var')
	fft_int = false;
end

Nt = length(t);							%	Number of time samples
dt = diff(t(1:2));						%	Temporal resolution

DN = ceil(tau_max/dt);					%	Number of delay samples to integrate

t1 = (DN:-1:1)'*dt;						%	Temporal padding for delay int
Et1 = Et(1)*exp(-(3*t1/tau_max).^2);	%	Append exponential to initial field
Et = [Et1; Et];

t = [t(1)-t1; t];						%	Padded temporal axis
Nt = length(t);							%	Number of time samples

At = -cumtrapz(t, Et);					%	Vector potential
Bt = cumtrapz(t, At);					%	Inegral of VP
B2t = .5*cumtrapz(t, At.^2);			%	Integral of .5*VP^2

indr = (DN+1:Nt)';						%	Recombination indices
indi = (1:DN)' + (0:Nt-DN-1);			%	Ionization indices

tau = t(DN+1)-t(1:DN);					%	Delays

Bt = Bt(indr).' - Bt(indi);

B2t = B2t(indr).' - B2t(indi);


%	Stationary momentum
ps = Bt ./ tau;

%	Action
S = ((.5*ps.^2 + Ip) .* tau) + B2t - ps.*Bt;

%	Dipole moment at ionization time t_f
dpf = dipole(ps - At(indr).', Ip);

%	Dipole moment at ionization time t_i
dpi = dipole(ps - At(indi), Ip);

A = (pi ./ (.5i*tau)).^1.5;				%	Amplitude factor

A = A .* dpf .* dpi .* Et(indi) .* exp(-1i*S);

%	Dipole moment
xt = -2*dt*trapz(imag(A)).';
end

function dp = Hydrogenlike(p, Ip)
alpha = 2*Ip;
A = 2^3.5 * alpha^1.25 / pi;
dp = A * p .* (p.^2 + alpha).^-3;
end

function dp = Gaussian(p, Ip) %#ok<*DEFNU>
alpha = 2*Ip;
A = (pi*alpha).^-.75 ./ alpha;
dp = A * p .* exp(-p.^2/(2*alpha));
end
