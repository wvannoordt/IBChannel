clear
clc
close all

%grid
glob.nguard = 2;
glob.nx = 100;
glob.ny = 100;
glob.L = 1.0;
glob.dx = glob.L/glob.nx;
glob.dy = glob.L/glob.ny;

%Variable names
glob.names{1} = 'P';
glob.names{2} = 'U';
glob.names{3} = 'V';
glob.names{4} = 'T';

%plot colors
colors{1} = [0.8 0.3 0.1];
colors{2} = [0.8 0.3 0.8];
colors{3} = [0.5 0.1 0.1];
colors{4} = [0.3 0.7 0.5];

%Initialize flow field
glob.nvars = length(glob.names);
glob.prims = zeros(glob.nvars, glob.nx+2*glob.nguard, glob.ny+2*glob.nguard);
glob.cons  = zeros(glob.nvars, glob.nx+2*glob.nguard, glob.ny+2*glob.nguard);

%properties
glob.mu = 1e-4;
glob.gamma = 1.4;
glob.cp = 1005;
glob.P0 = 1e5;
glob.umax = 40.0;
glob.Pr = 0.72;
glob.Twall = 300;
glob.R = glob.cp*(glob.gamma - 1) / glob.gamma;

%analytical solution, 1D
glob.ana1D = zeros(glob.nvars, glob.ny + 2*glob.nguard);
glob.y1D   = linspace(0-glob.nguard*glob.dy + 0.5*glob.dy, glob.L+glob.nguard*glob.dy-0.5*glob.dy, glob.ny + 2*glob.nguard);
H = glob.umax^2*glob.Pr/(2*glob.L^4*glob.cp);
M_mw = glob.umax/sqrt(glob.gamma*glob.R*glob.Twall);
glob.ana1D(1, :) = glob.P0;
glob.ana1D(2, :) = (glob.L^2 - glob.y1D.^2)*glob.umax;
glob.ana1D(3, :) = 0*glob.ana1D(3, :);
glob.ana1D(4, :) = glob.Twall - H*(glob.L^2-glob.y1D.^2).^2;


glob.psi = 0.5;
glob.conforms = (abs(glob.psi-0.5) < 1e-9);

%initialize flow
for ix = 1:glob.nx+2*glob.nguard
    glob.prims(1, ix, :) = glob.ana1D(1, :);
    glob.prims(2, ix, :) = glob.ana1D(2, :);
    glob.prims(3, ix, :) = glob.ana1D(2, :)*0.005;
    glob.prims(4, ix, :) = glob.ana1D(4, :);
end


Ntsteps = 100;
cfl = 0.5;
dt = cfl*min([glob.dx, glob.dy])/(glob.umax + sqrt(glob.gamma*glob.R*glob.Twall));
for nt = 1:Ntsteps
    [ghostswall, ghostscenter] = GetGhostValues(glob);
    glob.prims(:, :, (glob.nguard+glob.ny+1):(glob.nguard+glob.ny+glob.nguard)) = ghostswall(:,:,:);
    glob.prims(:, :, 1:glob.nguard) = ghostscenter(:,:,:);
    glob.prims(:, ((glob.nguard+glob.ny+1):(glob.nguard+glob.ny+glob.nguard)), :) = glob.prims(:, ((glob.nguard+1):(glob.nguard+glob.nguard)), :);
    glob.prims(:, (1:glob.nguard), :) = glob.prims(:, ((glob.ny+1):(glob.ny+glob.nguard)), :);
    glob.cons = Prims2Cons(glob);
    rhsConv = GetConvRhs(glob);
    rhsVisc = GetViscRhs(glob);
    rhsDiss = GetDissRhs(glob);
    glob.cons = glob.cons + dt*(rhsConv+rhsVisc+rhsDiss);
    glob.prims= Cons2Prims(glob);
end

%plots
%==========================================================================
%plot analytical solution
nxSampl = round(glob.nx/2);
prot = [1 glob.ny+2*glob.nguard];
lwAna = 3;
lwNum = 2;
figure
hold on
plot(glob.ana1D(1, :)/glob.P0, glob.y1D, '--', 'color', colors{1}, 'linewidth', lwAna)
plot(glob.ana1D(2, :)/glob.umax, glob.y1D, '--', 'color', colors{2}, 'linewidth', lwAna)
plot(glob.ana1D(3, :)/glob.umax, glob.y1D, '--', 'color', colors{3}, 'linewidth', lwAna)
plot((glob.ana1D(4, :)-glob.Twall)/glob.Twall, glob.y1D, '--', 'color', colors{4}, 'linewidth', lwAna)
plot(reshape(glob.prims(1, nxSampl, :), prot)/glob.P0, glob.y1D, 'color', colors{1}, 'HandleVisibility', 'off', 'linewidth', lwNum)
plot(reshape(glob.prims(2, nxSampl, :), prot)/glob.umax, glob.y1D, 'color', colors{2}, 'HandleVisibility', 'off', 'linewidth', lwNum)
plot(reshape(glob.prims(3, nxSampl, :), prot)/glob.umax, glob.y1D, 'color', colors{3}, 'HandleVisibility', 'off', 'linewidth', lwNum)
plot((reshape(glob.prims(4, nxSampl, :), prot)-glob.Twall)/glob.Twall, glob.y1D, 'color', colors{4}, 'HandleVisibility', 'off', 'linewidth', lwNum)
plot([-1 1], [glob.L glob.L], 'color', [0 0 0], 'HandleVisibility', 'off');

legend('P/P_0', 'U/U_{max}', 'V/U_{max}', '\Delta T/T_{wall}');
