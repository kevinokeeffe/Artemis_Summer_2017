%	Get folder names
dname = GetDataFolder('Artemis');
dnames = GetFilenames(dname, '*');
dnames = ExtractStrings(dnames, 'Scan_F10pt37_long_exp', true);
FNAME = @(dname1) cell2mat(ExtractStrings(GetFilenames(fullfile(dname, dname1)), ...
	'(Scan|Image)Data', true));
load(fullfile(dname, dnames{1}, FNAME(dnames{1})), 'I')

I = mean(I, 3);

%	Filter image and convert to binary
M = (medfilt2(I, [1 1]*4)>21);
M(end, :) = [];
I(end, :) = [];

%%

%	Get size of image
[Nth, Nw] = size(M);

%	Generate axis indices
nw = (1:Nw).';
nth = (1:Nth).';


%	Inegrate vertically & find edges of MCP
m = mean(M);
nw_edge = FindPeaks(1-m, 23, .4, 2);
nw0 = mean(nw_edge);

%	Remove outside MCP
M = M .* (nw>=nw_edge(1) & nw<=nw_edge(2)).';

%	Find width of MCP as function of vertical position
width = sum(M, 2);
[width, nth0] = max(width);

% %	Generate circle
% [nx, ny, th] = Geometry2D.Circle(.5*diff(nw_edge), 'Samples', 128);

MASK = @(nw0, nth0, r) (((nth-nth0).^2 + (nw.'-nw0).^2) <= r.^2);
ERR = @(m) sum(sum((~m).*M + m.*(~M)));

v = fminsearch(@(v) ERR(MASK(v(1), v(2), v(3))), [nw0, nth0, width/2]);

[nwc, nthc] = Geometry2D.Circle(v(3), 'Samples', 128);
nwc = nwc + v(1);
nthc = nthc + v(2);

fig = gcf;
clf(fig);
fig.Position(3:4) = [1200 600];
imagesc(nw, nth, M + MASK(v(1), v(2), v(3)) + scale(log10(I), 0))
grid on
axis equal

h = arrow2([v(1) v(2)], [v(1) v(2)] + v(3)*[cosd(10) -sind(10)], ...
	v(3)/100*[1 1 4 0], [0 .8 .8 1], [0 0 0], false);
line(nwc, nthc, 'Color', 'r', 'LineStyle', ':');
line(v(1) + [0 -1; 0 1]*v(3)/20, v(2) + [-1 0; 1 0]*v(3)/20, 'Color', 'r');
text(v(1)+v(3)/2, v(2)-10, sprintf('r = %.1fpxl = 20mm', v(3)), ...
	'FontSize', 20, 'Rotation', 10, 'HorizontalAlignment', 'Center');


ylim([-100 400]);
title('MCP Divergence Calibration')
xlabel('Frequency [pxl]');
ylabel('Divergence [pxl]');


L = 430 + 260;
Res_Space = 20/v(3) * 1e3;
Res_Div = Res_Space / L * 1e3;
text(100, -50, sprintf( ...
	['Spatial resolution: %.2f\\mum/pxl\n' ...
	'Divergence resolution (L\\sim%dmm): {\\bf\\color{red}%.2f\\murad/pxl}'], ...
	Res_Space, L, Res_Div), 'FontSize', 20);

print(fig, 'DivergenceCalibration', '-dpng', '-r300');
