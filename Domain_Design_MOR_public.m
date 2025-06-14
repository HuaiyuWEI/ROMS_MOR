% Huaiyu Wei
% This script is to set parameters and generate input files for an idealized
% mid-ocean ridge simulation using UCLA-ROMS to investigate eddy-driven topostrophy.


close all
clear all

% Add necessary paths for ROMS tools and postprocessing scripts
addpath(genpath('D:\OneDrive - University of California\Projects\Neptune\roms_tools_others'))
addpath('D:\OneDrive - University of California\MATLAB Codes\Neptune')
addpath('D:\OneDrive - University of California\MATLAB Codes\Neptune\Postprocess_V2')

% Unit conversions and basic settings
m1km = 1000;         % 1 km in meters
s1day = 86400;       % 1 day in seconds
fontsize = 15;       % Default font size for plots
LW = 1;              % Default line width



%% Set parameters of the simulation

% Define experiment name and create output directory
expname ='MOR_Res4_hr2000_tanhSm_Fm075Fs80Ft10'
output_dir = fullfile('E:\Data_ROMS\Neptune',expname);
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end
cd(output_dir)


% Forcing options: select one only
use_randForc = 1;    % 1 = use random wind forcing
use_wind = 0;        % 1 = use along-ridge wind forcing

if (use_randForc + use_wind > 1)
    error('Choose only one type of forcing: random or along-ridge wind.')
end


% Restore stratification at eastern/western boundaries (0 = no restoring)
restoring = 0;

% Barotropic run (0 = baroclinic simulation; 1 = barotropic simulation)
barotropic = 0;

% Mid-ocean ridge topography shape: choose only one
MOR_tanh        = 1;
MOR_linear      = 0;
MOR_halftanh    = 0;
MOR_doubletanh  = 0;
MOR_cos         = 0;
MOR_cos_v2      = 0;
if (MOR_cos + MOR_tanh + MOR_halftanh + MOR_doubletanh + MOR_cos_v2 + MOR_linear > 1)
    error('Select only one ridge topography shape.')
end

% Whether to apply smoothing to ensure smoothness at the ridge crest
topo_smooth = 1;

% Topographic parameters
hmax = 5000;         % Maximum domain depth (m)
h_ridge = 2000;      % Ridge height (m)

% Set to 1 if running a spin-down (unforced) simulation
spindown = 0;

% Use linearly stratified initial condition (0 = default: linear+exponential profile)
LinearT = 0;

% Wind forcing magnitude (only used when use_wind = 1)
% Positive: prograde; Negative: retrograde
tau0 = -0.05;

% Random forcing parameters (only used when use_randForc = 1)
RF_F0_rot = 0.75;              % Forcing amplitude
lambdaK   = 80 * m1km;         % Peak forcing length scale (m)
RF_tau    = 10 * s1day;        % Forcing autocorrelation timescale (s)

% Reference density
rho0 = 1000;

% Thermal expansion coefficient
alpha0 = 0.2/rho0;

%%% other important parameters
zob = 1e-3; %1mm bottom roughness
fconst = 8e-5; 
beta = 0; 

%%% Resolution and domain size
dx = 4 * m1km;
dy = dx;
if dx == 4 * m1km
    Nx = 502;
    Ny = 252;
end
Nlayer = 80;


%% generate the grid
Lx = (Nx-2)* dx ;
Ly = (Ny-2)* dx ;

% this function sets the ROMS grid Cartesian coordinates and metrics
G= uniform_grid (dx,dy,Nx-1,Ny-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Copyright (c) 2002-2023 The ROMS/TOMS Group
% Notice that ROMS has a staggered Arakawa's C-grid:
%
%    L          M                Number of PSI-points
%    Lr = L+1   Mr = M+1         Number of RHO-points
%    Lu = L     Mu = M+1         Number of U-points
%    Lv = L+1   Mv = M           Number of V-points
%
% On Input:
%
%    dx           Grid spacing in the X-direction (meters)
%    dy           Grid spacing in the Y-direction (meters)
%    L=Nx-1       Number of PSI grid-points in the X-direction
%    M=Ny-1       Number of PSI grid-points in the Y-direction
%
% On Output:
%
%    G            Uniform grid coordinates and metrics (struct array)
%
%                   G.x_psi    X-coordinates at PSI-points (meters)
%                   G.y_psi    Y-coordinates at PSI-points (meters)
%                   G.x_rho    X-coordinates at RHO-points (meters)
%                   G.y_rho    Y-coordinates at RHO-points (meters)
%                   G.x_u      X-coordinates at U-points (meters)
%                   G.y_u      Y-coordinates at U-points (meters)
%                   G.x_v      X-coordinates at V-points (meters)
%                   G.y_v      Y-coordinates at V-points (meters)
%                   G.pm       X-coordinates metric "m" (1/dx, 1/meters)
%                   G.pn       Y-coordinates metric "n" (1/dy, 1/meters)
%                   G.dndx     Inverse metric, d(1/pn)/d(x)
%                   G.dmde     Inverse metric, d(1/pm)/d(y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Lx ~= G.x_psi(end)
    error('Domain size in x is not equal to G.x_psi(end); check grid generation');
end

if Ly ~= G.y_psi(end)
    error('Domain size in y is not equal to G.y_psi(end); check grid generation');
end


X = G.x_rho;
Y = G.y_rho;

X_u = G.x_u;
Y_u = G.y_u;

X_v = G.x_v;
Y_v = G.y_v;

X_psi = G.x_psi;
Y_psi = G.y_psi;

%% Topography setup

% Define topography (etab) based on chosen ridge shape

if(MOR_tanh)
   
    etab = -hmax*ones(Nx,Ny);
    Ws = 1200*m1km*(1); % width of the mid-ocean ridge
    etab= -hmax + h_ridge/2 .*(1-tanh((abs(X-Lx/2)-Ws/4)/Ws*8));

elseif(MOR_cos)

    etab = -hmax + h_ridge/2 - h_ridge/2*cos(2*pi*X/Lx);

elseif(MOR_cos_v2)

    L_s = 1200*m1km;
    X_s = (Lx-L_s)/2;
    etab = -hmax + h_ridge/2 - h_ridge/2*cos(Lx/L_s*2*pi*(X-X_s)/Lx);
    etab(X<X_s)=-hmax;
    etab(X>Lx-X_s)=-hmax;

elseif(MOR_halftanh)

    etab = -hmax*ones(Nx,Ny);
    Ws = 300*m1km*(1); % slope half-width
    Rs_center = 998*m1km*(1); % slope mid-position
    X_away = abs(X - Rs_center);
    etab = -hmax + h_ridge - h_ridge*tanh((X_away)/Ws);

elseif(MOR_doubletanh)

    etab = -hmax*ones(Nx,Ny);
    Ws = 75*m1km; % slope half-width
    Rs_left1 = 850*m1km; % slope mid-position
    Rs_left2 = 500*m1km; % slope mid-position

    Rs_right1 = Lx-Rs_left1; % slope mid-position
    Rs_right2 = Lx-Rs_left2; % slope mid-position

    idx_left = X <= Lx/2;
    idx_right = X > Lx/2;

    etab(idx_left) = -hmax + h_ridge/2 + h_ridge/4*tanh((X(idx_left)-Rs_left1)/Ws) + h_ridge/4*tanh((X(idx_left)-Rs_left2)/Ws);
    etab(idx_right) = -hmax + h_ridge/2 - h_ridge/4*tanh((X(idx_right)-Rs_right1)/Ws)- h_ridge/4*tanh((X(idx_right)-Rs_right2)/Ws);;

elseif(MOR_linear)

    Ws = 600*m1km*(1);
    etab = -hmax + h_ridge - h_ridge*   abs(X - Lx/2)/(Ws);
    etab(etab<-hmax) = -hmax;

end


%%% Optional: Smooth the topography cross-section

if(topo_smooth)
    etab_raw = etab(:,1);
    etab_smooth = smoothn(etab_raw,50,'robust');
    figure
    plot(etab_raw)
    hold on
    plot(etab_smooth)
    figure
    plot(diff(etab_raw),'.')
    hold on
    plot(diff(etab_smooth),'.')
    title('diff h')
    legend('original','smoothed')
    etab = repmat(etab_smooth,[1,Ny]);
end


%%% Plot 3D view of the bathymetry
fignum = 100;
fignum = fignum+1;
figure(fignum);
surf(X,Y,etab);
shading interp;
colorbar;
% plot(yy_h,etab(1,:))
pbaspect([2,1,1]);



%%% Plot cross-section of the ridge (x-z view)
figure('Position',[150 50 800 400]);
plot(X(:,1)/ 1000,etab(:,1),'--','LineWidth',LW+1);
xlabel('Offshore distance (km)','Interpreter','Latex');
ylabel('-Depth (m)','Interpreter','Latex');
title('Bathymetry;  cross-section view','Interpreter','Latex');
xlim([X(2) X(end-1)]/ 1000)
set(gca,'fontsize',fontsize-3)
set(gca,'layer','top','linewidth',LW)
grid on
set(gca,'GridLineStyle','--')
exportgraphics(gcf,'bathy_crossview.png','resolution',200)


%%% Plot full-domain view of the bathymetry
h_rho = etab;
mask_rho = ones(Nx,Ny);
% mask_rho(:,[1:2 Ny-1:Ny]) = 0;
if(use_wind)
    mask_rho([1:2 Nx-1:Nx],:) = 0;
end
figure('Position',[50 50 800 800]);
pcolor(X / 1000, Y / 1000, h_rho.* mask_rho);
hold on
[C,hh]=contour(X / 1000, Y / 1000, abs(h_rho).* mask_rho, [0:250:6000],'-k','linewidth',LW-0.5);
clabel(C,hh,'color','k','fontsize',fontsize-8,'labelspacing',3000)
colormap('haxby')
axis equal
set(gca,'fontsize',fontsize)
set(gca,'layer','top','linewidth',LW)
shading interp;
colorbar;
xlabel('X (km)','Interpreter','Latex');
ylabel('Y (km)','Interpreter','Latex');
clim([min(h_rho,[],'all') 0])
title('Full domain view','Interpreter','Latex');
exportgraphics(gcf,'bathy_fullview.png','resolution',200)





%% generate and plot vertical grid


% Vertical stretching parameters (ROMS terrain-following coordinates)
theta_s = 5;
theta_b = 7;
hc = 300;

% Generate vertical levels (r- and w-points) using ROMS stretching
[z3D_r,Cs_r] = zlevs4(-h_rho,0*h_rho,theta_s,theta_b,hc,Nlayer,'r','new2012');
[z3D_w,Cs_w] = zlevs4(-h_rho,0*h_rho,theta_s,theta_b,hc,Nlayer,'w','new2012');


%Compute vertical grid cell thickness at each level (in meters)
Dz = (z3D_w(2:end,:,:) - z3D_w(1:end-1,:,:)); % dimension of z3D is Nlayer, Nx, Ny



%%% Plot Vertical Grid Cross-section
index=1; % the y index at which the section will be plotted
plt=2;
Zzoom=120;
x = G.x_rho;
z = squeeze(z3D_r(:,:,index));
[Lp,Mp]=size(h_rho);
figure;
set(gcf,'Units','Normalized',...
    'Position',[0.2 0.1 0.6 0.8],...
    'PaperOrientation', 'landscape', ...
    'PaperUnits','Normalized',...
    'PaperPosition',[0.2 0.1 0.6 0.8]);

xi=x(:,index)';
xi2=[xi(1) xi xi(Lp)];
hplot=h_rho(:,index)';
zmin=min(hplot);
hplot=[zmin hplot zmin];
if (plt == 2)
    h1=subplot(2,1,1);
    p1=get(h1,'pos');
    p1(2)=p1(2)+0.1;
    p1(4)=p1(4)-0.1;
    set (h1,'pos',p1);
end
hold off;
fill(xi2,hplot,[0.6 0.7 0.6]);
hold on;
han1=plot(xi,z);
set(han1,'color', [0.5 0.5 0.5]);
if (plt == 2)
    set(gca,'xlim',[-Inf Inf],'ylim',[-abs(Zzoom) 0]);
    ylabel('-depth  (m)');
else
    set(gca,'xlim',[-Inf Inf],'ylim',[zmin 0]);
end

title(['Grid (\rho-points) Section at  y = ',num2str(index*dy/m1km),'km']);

if (plt == 2)
    h2=subplot(2,1,2);
    p2=get(h2,'pos');
    p2(4)=p2(4)+0.2;
    set (h2,'pos',p2);

    hold off;
    fill(xi2,hplot,[0.6 0.7 0.6]);
    hold on;
    han2=plot(xi,z);
    set(han2,'color', [0.5 0.5 0.5]);
    set(gca,'xlim',[-Inf Inf],'ylim',[zmin 0]);
end

xlabel({['   \theta_s = ' num2str(theta_s), ...
    '   \theta_b = ' num2str(theta_b), ...
    '   hc  = ' num2str(hc), ...
    '   N = ' num2str(Nlayer)],['x-axis']});
ylabel('-depth  (m)');

exportgraphics(gcf,'VerticalGrid.png','resolution',200)
saveas(gcf,'VerticalGrid.fig')


%%% Plot vertical cell thickness profile at slope center
figure('Position',[50 50 600 800]);
indx = round(Nx/4);
plot(Dz(:,indx,index),[1:Nlayer],'-','LineWidth',1.5)
hold on
plot(Dz(:,indx,index),[1:Nlayer],'.','MarkerSize',10)
set(gca,'fontsize',fontsize)
set(gca,'layer','top','linewidth',LW)
grid on
xlabel('Cell thickness (m)','Interpreter','Latex');
ylabel('Sigma level; 0 stands for bottom','Interpreter','Latex');
title('Cell thickness at slope center','Interpreter','Latex');
exportgraphics(gcf,'VerticalGrid_dz.png','resolution',200)


%%% Estimate and Plot Bottom Drag Coefficient (Cd)
Dz3D_r = z3D_w(2:end,:,:) - z3D_w(1:end-1,:,:);
figure('Position',[50 50 800 800]);
Dz_bottom = squeeze(Dz3D_r(1,:,:));
Cd = (0.41./log(1+0.5*Dz_bottom./zob)).^2;
pcolor(X / 1000, Y / 1000, Cd );
shading flat
colorbar
% clim([2 5]*1e-2)
cmocean('shad',100)
axis equal
set(gca,'fontsize',fontsize)
set(gca,'layer','top','linewidth',LW)
xlabel('X (km)','Interpreter','Latex');
ylabel('Y (km)','Interpreter','Latex');
title('Cd','Interpreter','Latex');
xlim([0 Lx]/1000)
ylim([0 Ly]/1000)
exportgraphics(gcf,'Cd.png','resolution',200)



%%% Coriolis parameter
f_rho = fconst + Y * beta;
figure('Position',[50 50 800 800]);
pcolor(X / 1000, Y / 1000, f_rho)
hold on
[C,hh]=contour(X / 1000, Y / 1000, f_rho, 5,'-k','linewidth',LW);
clabel(C,hh,'color','k','fontsize',fontsize-5,'labelspacing',3000)
shading interp
colorbar
clim([0.6 1]*1e-4)
cmocean('shad',100)
axis equal
set(gca,'fontsize',fontsize)
set(gca,'layer','top','linewidth',LW)
xlabel('X (km)','Interpreter','Latex');
ylabel('Y (km)','Interpreter','Latex');
title('Coriolis frequency','Interpreter','Latex');
exportgraphics(gcf,'Coriolis.png','resolution',200)



%% Check Grid Stiffness Ratios in the terrain-following coordinate (rx0 and rx1)

% Generate U, V, and PSI point masks from the RHO-point mask
[mask_u, mask_v, mask_psi] = rt.uvp_masks(mask_rho);

% Initialize stiffness ratio arrays
my_rx0_x = zeros(size(mask_u));  % rx0 in x-direction
my_rx1_x = zeros(size(mask_u));  % rx1 in x-direction
my_rx0_y = zeros(size(mask_v));  % rx0 in y-direction
my_rx1_y = zeros(size(mask_v));  % rx1 in y-direction

% Rearrange dimensions for z_w: [Nx, Ny, Nlayer+1]
z_w = permute(z3D_w, [2, 3, 1]);

% -------- Compute rx0 and rx1 in x-direction --------
for j = 1:Ny-1
    for i = 2:Nx-1
        if mask_u(i, j) > 0.5
            % rx0: surface slope normalized by total depth
            numerator = z_w(i, j, 1) - z_w(i-1, j, 1);
            denominator = z_w(i, j, 1) + z_w(i-1, j, 1);
            my_rx0_x(i, j) = abs(numerator / denominator);

            % rx1: maximum vertical slope across layers
            temp = 0;
            for k = 1:Nlayer
                num = (z_w(i, j, k+1) - z_w(i-1, j, k+1)) + ...
                      (z_w(i, j, k)   - z_w(i-1, j, k));
                denom = (z_w(i, j, k+1) + z_w(i-1, j, k+1)) - ...
                        (z_w(i, j, k)   + z_w(i-1, j, k));
                temp = max(temp, abs(num / denom));
            end
            my_rx1_x(i, j) = temp;
        end
    end
end

% -------- Compute rx0 and rx1 in y-direction --------
for j = 3:Ny-1
    for i = 1:Nx-1
        if mask_v(i, j) > 0.5
            numerator = z_w(i, j, 1) - z_w(i, j-1, 1);
            denominator = z_w(i, j, 1) + z_w(i, j-1, 1);
            my_rx0_y(i, j) = abs(numerator / denominator);

            temp = 0;
            for k = 1:Nlayer
                num = (z_w(i, j, k+1) - z_w(i, j-1, k+1)) + ...
                      (z_w(i, j, k)   - z_w(i, j-1, k));
                denom = (z_w(i, j, k+1) + z_w(i, j-1, k+1)) - ...
                        (z_w(i, j, k)   + z_w(i, j-1, k));
                temp = max(temp, abs(num / denom));
            end
            my_rx1_y(i, j) = temp;
        end
    end
end


%%% Plot stiffness ratios (rx0, rx1)

figure('Position', [50 50 1600 800]);

% --- rx0 in y-direction ---
subplot(2,2,1)
pcolor(X_v / 1000, Y_v / 1000, my_rx0_y); shading flat; hold on
[C, hh] = contour(X / 1000, Y / 1000, abs(h_rho) .* mask_rho, ...
                  [110 300 1500 3000 3500 4000], '-k', 'LineWidth', LW - 0.5);
clabel(C, hh, 'Color', 'k', 'FontSize', fontsize - 8, 'LabelSpacing', 3000)
cmocean('shad', 20)
axis equal
clim([0, max(my_rx0_y, [], 'all') + 0.001])
set(gca, 'FontSize', fontsize, 'Layer', 'top', 'LineWidth', LW)
colorbar
xlabel('X (km)', 'Interpreter', 'Latex')
ylabel('Y (km)', 'Interpreter', 'Latex')
title(['rx0 in y-direction; max = ', num2str(max(my_rx0_y, [], 'all'), 4)], 'Interpreter', 'Latex')

% --- rx1 in y-direction ---
subplot(2,2,2)
pcolor(X_v / 1000, Y_v / 1000, my_rx1_y); shading flat; hold on
[C, hh] = contour(X / 1000, Y / 1000, abs(h_rho) .* mask_rho, ...
                  [110 300 1500 3000 3500 4000], '-k', 'LineWidth', LW - 0.5);
clabel(C, hh, 'Color', 'k', 'FontSize', fontsize - 8, 'LabelSpacing', 3000)
cmocean('shad', 20)
axis equal
clim([0, max(my_rx1_y, [], 'all') + 0.001])
set(gca, 'FontSize', fontsize, 'Layer', 'top', 'LineWidth', LW)
colorbar
xlabel('X (km)', 'Interpreter', 'Latex')
ylabel('Y (km)', 'Interpreter', 'Latex')
title(['rx1 in y-direction; max = ', num2str(max(my_rx1_y, [], 'all'), 4)], 'Interpreter', 'Latex')

% --- rx0 in x-direction ---
subplot(2,2,3)
pcolor(X_u / 1000, Y_u / 1000, my_rx0_x); shading flat; hold on
[C, hh] = contour(X / 1000, Y / 1000, abs(h_rho) .* mask_rho, ...
                  [110 300 1500 3000 3500 4000], '-k', 'LineWidth', LW - 0.5);
clabel(C, hh, 'Color', 'k', 'FontSize', fontsize - 8, 'LabelSpacing', 3000)
cmocean('shad', 20)
axis equal
clim([0, max(my_rx0_x, [], 'all') + 0.001])
set(gca, 'FontSize', fontsize, 'Layer', 'top', 'LineWidth', LW)
colorbar
xlabel('X (km)', 'Interpreter', 'Latex')
ylabel('Y (km)', 'Interpreter', 'Latex')
title(['rx0 in x-direction; max = ', num2str(max(my_rx0_x, [], 'all'), 4)], 'Interpreter', 'Latex')

% --- rx1 in x-direction ---
subplot(2,2,4)
pcolor(X_u / 1000, Y_u / 1000, my_rx1_x); shading flat; hold on
[C, hh] = contour(X / 1000, Y / 1000, abs(h_rho) .* mask_rho, ...
                  [110 300 1500 3000 3500 4000], '-k', 'LineWidth', LW - 0.5);
clabel(C, hh, 'Color', 'k', 'FontSize', fontsize - 8, 'LabelSpacing', 3000)
cmocean('shad', 20)
axis equal
clim([0, max(my_rx1_x, [], 'all') + 0.001])
set(gca, 'FontSize', fontsize, 'Layer', 'top', 'LineWidth', LW)
colorbar
xlabel('X (km)', 'Interpreter', 'Latex')
ylabel('Y (km)', 'Interpreter', 'Latex')
title(['rx1 in x-direction; max = ', num2str(max(my_rx1_x, [], 'all'), 4)], 'Interpreter', 'Latex')

% Save the figure
exportgraphics(gcf, 'Stiffness.png', 'resolution', 200)



%% generate random wind forcing function
%%% please refer to Stewart et al. (2024) for the fomulation of the random forcing
%%% https://doi.org/10.1175/JPO-D-23-0196.1

if(use_randForc)
    
    % Total number of forcing snapshots (e.g., daily forcing for ~5 years)
    Nfrc = 365*5 + 2; % +2 to account for leap years

    %%% set random seed
    rng(42)

    % Define zonal wavenumbers
    if Nx/2 == round(Nx/2) 
        k = [0:1:Nx/2-1,-Nx/2:1:-1]; %%% Zonal wavenumber
    else
        k = [0:1:round(Nx/2)-1,-round(Nx/2)+1:1:-1]; %%% Zonal wavenumber
    end
    K_xk = 2*pi.*(k)./Lx;

    % Define meridional wavenumbers
    if Ny/2 == round(Ny/2)
        l = [0:1:Ny/2-1,-Ny/2:1:-1]; %%% Meridional wavenumber
    else
        l = [0:1:round(Ny/2)-1,-round(Ny/2)+1:1:-1]; %%% Meridional wavenumber
    end
    K_yl = 2*pi.*(l)./Ly;

    [K_ykl,K_xkl] = meshgrid(K_yl, K_xk);   
    K = sqrt(K_xkl.^2 + K_ykl.^2); %%% Absolute wavenumber

    K_0 = 2*pi/lambdaK; %%% Most strongly forced wavenumber
    W = K_0/8; %%% Exponential width of energy band in wavenumber space
    RF_mask_fft = 2*exp(-((K-K_0)/W).^2);
    RF_mask_norm = sqrt(0.5*sum(RF_mask_fft(:).^2));
    RF_kk_recip = 1.0 ./ K;
    RF_kk_recip(1,1) = 0;
    RFrot_fft = complex(zeros(Nx,Ny));
    RFamp_rot_fft = RF_F0_rot .* RF_kk_recip .* RF_mask_fft ./ RF_mask_norm;

    %%% Spin-up of random forcing function (to reach statistical equilibrium)
    %%% Evolve over many forcing decorrelation time scales to get a representative forcing field
    %%% The random forcing function during this "spin up" phase will not be used
    deltaT = 864;
    Tratio = deltaT/RF_tau;
    for n=1:20000
        RF_phase = 2*pi*rand(Nx,Ny);
        RFrot_fft = RFrot_fft - Tratio.*RFrot_fft + sqrt(2*Tratio).*exp(1i*RF_phase).*RFamp_rot_fft;
        temp2 = Nx*Ny*real(ifft2(RFrot_fft));
        temp3 = ((temp2(1:Nx,1:Ny,:)-temp2([Nx 1:Nx-1],1:Ny,:))/dx).^2 + ...
            ((temp2(1:Nx,1:Ny,:)-temp2(1:Nx,[Ny 1:Ny-1],:))/dy).^2;
        temp(n) = sqrt(mean(temp3,[1,2]));
    end


    % Plot spin-up series of the forcing magnitude
    figure
    plot(squeeze(temp))
    xlabel('Iteration')
    ylabel('Domain averaged RMS forcing magnitude')
    exportgraphics(gcf,'RF_spinup.png','resolution',200)
    clear temp


    %%% Generate Random Forcing that will be used in the simulation
    deltaT = 86400;
    Tratio = deltaT/RF_tau;
    for n=1:Nfrc
        RF_phase = 2*pi*rand(Nx,Ny);
        RFrot_fft = RFrot_fft - Tratio.*RFrot_fft + sqrt(2*Tratio).*exp(1i*RF_phase).*RFamp_rot_fft;
        RFrot_alltime(:,:,n) = Nx*Ny*real(ifft2(RFrot_fft)) ;
    end


    %%% Check if the RMS forcing magnitude matches specified forcing amplitude
    dRFrot_dx = (RFrot_alltime(1:Nx,1:Ny,:)-RFrot_alltime([Nx 1:Nx-1],1:Ny,:))/dx;
    dRFrot_dy = (RFrot_alltime(1:Nx,1:Ny,:)-RFrot_alltime(1:Nx,[Ny 1:Ny-1],:))/dy;
    gradsq_dRFrot = dRFrot_dx.^2+dRFrot_dy.^2;
    sqrt(mean(gradsq_dRFrot(:)))
    figure
    plot(squeeze(sqrt(mean(gradsq_dRFrot,[1,2]))))
    xlabel('Iteration')
    ylabel('Domain averaged RMS forcing magnitude')
    xlim([0 Nfrc])
    exportgraphics(gcf,'RF_used.png','resolution',200)

    clear tauy_dx taux_dy dRFrot_dx dRFrot_dy gradsq_dRFrot
   

    %% Convert random forcing function to Wind Stress Components (tau_u, tau_v)

    % Assume RFrot_alltime is defined on PSI points with ghost points
    taux_u = - (RFrot_alltime(1:Nx,1:Ny,:)-RFrot_alltime(1:Nx,[Ny 1:Ny-1],:))/dy;
    gradsq_dRFrot = taux_u.^2;
    taux_u = taux_u(2:end,:,:);
    
    tauy_v = (RFrot_alltime(1:Nx,1:Ny,:)-RFrot_alltime([Nx 1:Nx-1],1:Ny,:))/dx;
    gradsq_dRFrot = gradsq_dRFrot + tauy_v.^2;
    tauy_v = tauy_v(:,2:end,:);

    % Verify RMS wind stress magnitude
    sqrt(mean(gradsq_dRFrot(:)))


    %%% Apply Forcing Mask Near Restoring Boundaries (if restoring is active)
    if(restoring)
        L_nudg = 50*m1km;
        Lramp_RF = 30*m1km;

        RF_mask_u = ones(size(mask_u));
        idx_0 = (X_u <= L_nudg);
        idx_cos = (X_u >= L_nudg  & X_u < L_nudg+Lramp_RF);
        RF_mask_u(idx_0) = 0;
        RF_mask_u(idx_cos) = 0.5.*( 1-cos(pi*(X_u(idx_cos)-L_nudg)/Lramp_RF));

        idx_0 = (X_u >= Lx-L_nudg);
        idx_cos = (X_u >= Lx-L_nudg-Lramp_RF  & X_u < Lx-L_nudg );
        RF_mask_u(idx_0) = 0;
        RF_mask_u(idx_cos) = 0.5.*( 1-cos(pi*(X_u(idx_cos)-(Lx-L_nudg))/Lramp_RF));

        RF_mask_v = ones(size(mask_v));
        idx_0 = (X_v <= L_nudg);
        idx_cos = (X_v >= L_nudg  & X_v < L_nudg+Lramp_RF);
        RF_mask_v(idx_0) = 0;
        RF_mask_v(idx_cos) = 0.5.*( 1-cos(pi*(X_v(idx_cos)-L_nudg)/Lramp_RF));

        idx_0 = (X_v >= Lx-L_nudg);
        idx_cos = (X_v >= Lx-L_nudg-Lramp_RF  & X_v < Lx-L_nudg );
        RF_mask_v(idx_0) = 0;
        RF_mask_v(idx_cos) = 0.5.*( 1-cos(pi*(X_v(idx_cos)-(Lx-L_nudg))/Lramp_RF));

        taux_u = taux_u.* repmat(RF_mask_u,[1,1,Nfrc]);
        tauy_v = tauy_v.* repmat(RF_mask_v,[1,1,Nfrc]);

    end


    % Check if there is net mean momentum input
    figure
    plot(squeeze(mean(taux_u,[1,2])))
    hold on
    plot(squeeze(mean(tauy_v,[1,2])))
    title('domain averaged wind stress')
    legend('tau_x','tau_y')
    exportgraphics(gcf,'RF_sum.png','resolution',200)
    % time-mean domain-averaged tau_x and tau_y
    mean_tau_x = mean(taux_u, [1, 2, 3]);
    mean_tau_y = mean(tauy_v, [1, 2, 3]);
    disp(['Time-mean domain-averaged tau_x: ', num2str(mean_tau_x)]);
    disp(['Time-mean domain-averaged tau_y: ', num2str(mean_tau_y)]);


    %% Make an animation of the random forcing function and forcing components


    crange3 = 10 * dx * [-RF_F0_rot RF_F0_rot];  % for RFrot_alltime
    fontsize = 18;
    LW = 1;

    figure('Position', [50 50 1600 800]);

    for n = 1 %:360
        clf

        % Plot RFrot_alltime as background color
        pcolor(X_psi / 1000, Y_psi / 1000, RFrot_alltime(1:end-1,1:end-1,n));
        hold on
        clim(crange3)
        cmocean('red', 100, 'pivot', 0)
        shading interp
        axis equal

        grid off;box on
        set(gca,'layer','top','linewidth',LW,'FontSize',fontsize,'tickLabelinterpreter', 'latex');

        xlabel('Cross-ridge distance (km)', 'Interpreter', 'Latex')
        ylabel('Along-ridge distance (km)', 'Interpreter', 'Latex')
        xlim([0 Lx] / 1000)
        ylim([0 Ly] / 1000)
        % Overlay wind stress vectors
        step = 5; % Adjust for vector density
        quiver(X_u(1:step:end,1:step:end) / 1000, Y_u(1:step:end,1:step:end) / 1000, ...
            taux_u(1:step:end,1:step:end,n), tauy_v(1:step:end,1:step:end,n), ...
            1.5, 'color',[.4 .4 .4], 'LineWidth', LW-0.2,'MaxHeadSize',1);
        colorbar('Position',[0.915 0.13 0.01 0.775],'tickLabelinterpreter', 'latex')
        title(['Random force function; Day ' num2str(n)  ], 'FontSize', fontsize+2, 'Interpreter', 'latex');
        % Save figure
        folderName = './RF';
        if ~exist(folderName, 'dir')
            mkdir(folderName)
        end
        exportgraphics(gcf, [folderName, '/', num2str(n), '.png'], 'Resolution', 200)
    end

    
    clear tau_curl_psi  RFrot_alltime


    % %%% save as mp4
    % n = 1;
    % vidobj = VideoWriter([folderName '/RandomForcing'],'MPEG-4');
    % vidobj.FrameRate = 15;
    % vidobj.Quality = 100;
    % open(vidobj);
    % for i = 1:1:360
    %    imgname = strcat([folderName,'/',num2str(n),'.png'])
    %   if (exist(imgname,'file') ~= 2)
    %     'Ran out of images!'
    %     break;
    %   else
    %     imgdata = imread(imgname);
    %     writeVideo(vidobj,imgdata);
    %     n = n+1;
    %   end
    % end
    % close(vidobj);
    % 'Save avi successfully'


%%  if use along-ridge wind forcing instead of random wind

elseif(use_wind)

    Nfrc = 2;
    tauy_v = 0*Y_v;
    Lwind = Ws*2;
    idx_left = and(X_v <= Lx/2,X_v>=Lx/2-Lwind);
    idx_right = and(X_v > Lx/2,X_v <=Lx/2+Lwind);
    tauy_v(idx_left) = tau0*sin(pi*(Lx/2-X_v(idx_left))/Lwind).^2;
    tauy_v(idx_right) = -tau0*sin(pi*(X_v(idx_right)-Lx/2)/Lwind).^2;

    taux_u = X_u*0;
    figure
    plot(tauy_v(:,1),'.')
    figure
    plot(diff(tauy_v(:,1)),'.')

    tau_x_rho = rt.u2rho(taux_u);
    tau_y_rho = rt.v2rho(tauy_v);
    figure

    figure('Position',[50 50 800 800]);
    pcolor(X / 1000, Y / 1000, h_rho .* mask_rho);
    hold on
    colormap('haxby')
    axis equal
    set(gca,'fontsize',fontsize)
    set(gca,'layer','top','linewidth',LW)
    shading interp;
    colorbar;
    xlabel('X (km)','Interpreter','Latex');
    ylabel('Y (km)','Interpreter','Latex');
    title('Wind stress; full domain view','Interpreter','Latex');
    hold on
    magnitude = sqrt(tau_x_rho.^2 + tau_y_rho.^2);
    quiver(X(1:10:end,1:10:end)/ 1000,Y(1:10:end,1:10:end)/ 1000,tau_x_rho(1:10:end,1:10:end),tau_y_rho(1:10:end,1:10:end),2, ...
        'Color', [0 0 0],...
        'LineWidth', 1, 'MaxHeadSize', 0.8, 'AutoScale', 'on')
    exportgraphics(gcf,'wind.png','resolution',200)

    taux_u = repmat(taux_u,[1,1,Nfrc])
    tauy_v = repmat(tauy_v,[1,1,Nfrc])

end



%% Initial Stratification Setup

g0 = 9.81;               %%% Gravitational acceleration (m/s^2)
hs = 400;                %%% Exponential decay scale for temperature profile
hl = 40*m1km;            %%% Linear decay scale for temperature profile
TBot = 0;                %%% Sea floor temperature
TSurf = 25;              %%% Sea surface temperature

% For barotropic run, set constant temperature
if(barotropic)
    TSurf = TBot;
end

% Total temperature difference (negative)
TRange = TBot-TSurf;

%%% Temperature profile on geopotential coordinate
hmax_rand = max(abs(z3D_w),[],'all');
if(LinearT)
    T_init = TSurf - TRange*z3D_r/hmax_rand;
    T_init_w = TSurf - TRange*z3D_w/hmax_rand;
else
    T_init = TBot - TRange*(exp(z3D_r/hs)+z3D_r/hl-exp(-hmax_rand/hs)+hmax_rand/hl)/(1-exp(-hmax_rand/hs)+hmax_rand/hl);
    T_init_w =  TBot - TRange*(exp(z3D_w/hs)+z3D_w/hl-exp(-hmax_rand/hs)+hmax_rand/hl)/(1-exp(-hmax_rand/hs)+hmax_rand/hl);
end


%%% Plot Initial Temperature Profile at Deepest Point

indzplot = find(abs(z3D_w(1,:,:)) == hmax_rand,1);

figure('Position',[150 50 1200 600]);
subplot(1,3,1)
plot(T_init(:,indzplot),z3D_r(:,indzplot),'LineWidth',LW+1);
title('Initial temperature','Interpreter','Latex');
xlabel('$\theta$ ($^\circ$C)','Interpreter','Latex');
ylabel('Depth (m)','Interpreter','Latex');
y_ticks = [-hmax:500:0];
yticks(y_ticks);
yticklabels(arrayfun(@(x) num2str(-x), y_ticks, 'UniformOutput', false));
ylim([-hmax_rand 0])
set(gca,'fontsize',fontsize)
set(gca,'layer','top','linewidth',LW)
grid on
set(gca,'GridLineStyle','--')
set(gca,'Position',[0.09 0.14 0.18 0.8])
xticks([TBot:5:TSurf])
text(0.8,0.05,'(a)','Units','normalized','FontSize',fontsize+5,'Interpreter','latex');


%%% calculate buoyancy frequency and deformation radius

N2_init =alpha0*g0 *( T_init_w(2:end,:,:)-T_init_w(1:end-1,:,:) ) ./ (z3D_w(2:end,:,:)-z3D_w(1:end-1,:,:));
N_init = sqrt(N2_init);

%%% Plot Buoyancy Frequency Profile
subplot(1,3,2)
plot(N_init(:,indzplot),z3D_r(:,indzplot),'LineWidth',LW+1);
xlabel('$N$ (s$^{-1}$)','Interpreter','Latex');
title('Buoyancy frequency','Interpreter','Latex');
yticks(y_ticks);
ylim([-hmax_rand 0])
yticklabels('')
set(gca,'fontsize',fontsize)
set(gca,'layer','top','linewidth',LW)
set(gca,'xscale','log')
grid on
set(gca,'GridLineStyle','--')
xlim([1e-3 2e-2])
ylim([-hmax 0])
set(gca,'Position',[0.3 0.14 0.18 0.8])
text(0.8,0.05,'(b)','Units','normalized','FontSize',fontsize+5,'Interpreter','latex');

%%% Compute First Baroclinic Deformation Radius Ld (km)
Ld = squeeze(sum(N_init .* (z3D_w(2:end,:,:)-z3D_w(1:end-1,:,:)),1)/fconst/pi/1000);
Ld(Ld==0)=nan ;%for plotting purpose

%%% Plot Zonal Cross-Section of Deformation Radius
subplot(1,3,3)
plot(X(:,1) / 1000, mean(Ld(:,:),2),'LineWidth',LW+1);
hold on
xlabel('X (km)','Interpreter','Latex');
ylabel('$L_\mathrm{d}$ (km)','Interpreter','Latex');
title('Deformation radius;  cross-section view','Interpreter','Latex');
xlim([0 Lx/1000])
set(gca,'fontsize',fontsize)
set(gca,'layer','top','linewidth',LW)
grid on
set(gca,'GridLineStyle','--')
set(gca,'Position',[0.55 0.14 0.4 0.8])
text(0.02,0.05,'(c)','Units','normalized','FontSize',fontsize+5,'Interpreter','latex');

exportgraphics(gcf,'Init_stratification.png','resolution',200)




% (Nlay,Nx,Ny) => (Nx,Ny,Nlay)
T_init = permute(T_init,[2,3,1]);


%% resotring
if(restoring)
    nudg_Time = 1 * s1day;

    t_targ = T_init;
    nudg_mask_2d = zeros(Nx,Ny);

    %%% restore temperature at zonal boundaries
    idx_r = (X(1:Nx,1) >= Lx-L_nudg);
    nudg_mask_2d(idx_r,:) =( 1- (Lx - X(idx_r,:))./L_nudg ) .* mask_rho(idx_r,:);
    idx_l = (X(1:Nx,1) <= L_nudg);
    nudg_mask_2d(idx_l,:) =( (L_nudg-(X(idx_l,:)))./L_nudg ) .* mask_rho(idx_l,:);

    nudg_mask_3d = repmat(nudg_mask_2d,[1 1 Nlayer]);
    nudgc = nudg_mask_3d .* 1/nudg_Time;
    figure
    pcolor(squeeze(nudg_mask_3d(:,:,1))')
    colorbar
    title('nudging mask')

    figure
    pcolor(squeeze(nudgc(:,:,1))')
    colorbar
    title('nudging coefficient s-1')
end


%% Generate an empty netcdf file of grid
if ~exist('./Neptune_input', 'dir')
    mkdir('./Neptune_input')
end
filename = './Neptune_input/neptune_grid.nc';

if exist(filename, 'file')
    disp(['Deleting existing grid file: ', filename]);
    delete(filename);
end

% Create NetCDF file
ncid = netcdf.create(filename, 'CLOBBER');
%the 'CLOBBER' mode allows the new file to overwrite any existing file without warning.

% Define dimensions
dimid_xi_psi = netcdf.defDim(ncid, 'xi_psi', Nx-1);
dimid_xi_rho = netcdf.defDim(ncid, 'xi_rho', Nx);
dimid_xi_u = netcdf.defDim(ncid, 'xi_u', Nx-1);
dimid_xi_v = netcdf.defDim(ncid, 'xi_v', Nx);
dimid_eta_psi = netcdf.defDim(ncid, 'eta_psi', Ny-1);
dimid_eta_rho = netcdf.defDim(ncid, 'eta_rho', Ny);
dimid_eta_u = netcdf.defDim(ncid, 'eta_u', Ny);
dimid_eta_v = netcdf.defDim(ncid, 'eta_v', Ny-1);

if(restoring)
    dimid_s_rho = netcdf.defDim(ncid, 's_rho', Nlayer);
    dimid_s_w = netcdf.defDim(ncid, 's_w', Nlayer+1);
end

% Define variables
% Variable 'spherical'
varid_spherical = netcdf.defVar(ncid, 'spherical', 'NC_CHAR', []);
netcdf.putAtt(ncid, varid_spherical, 'long_name', 'Grid type logical switch');
netcdf.putAtt(ncid, varid_spherical, 'option(T)', 'spherical');
netcdf.putAtt(ncid, varid_spherical, 'option(F)', 'Cartesian');

% Variable 'xl'
varid_xl = netcdf.defVar(ncid, 'xl', 'NC_DOUBLE', []);
netcdf.putAtt(ncid, varid_xl, 'long_name', 'domain length in the XI-direction');
netcdf.putAtt(ncid, varid_xl, 'units', 'meter');

% Variable 'el'
varid_el = netcdf.defVar(ncid, 'el', 'NC_DOUBLE', []);
netcdf.putAtt(ncid, varid_el, 'long_name', 'domain length in the ETA-direction');
netcdf.putAtt(ncid, varid_el, 'units', 'meter');

% Variable 'h'
varid_h = netcdf.defVar(ncid, 'h', 'NC_DOUBLE', [dimid_xi_rho, dimid_eta_rho]);
netcdf.putAtt(ncid, varid_h, 'long_name', 'Final bathymetry at RHO-points');
netcdf.putAtt(ncid, varid_h, 'units', 'meter');
netcdf.putAtt(ncid, varid_h, 'coordinates', 'lon_rho lat_rho');
netcdf.putAtt(ncid, varid_h, 'field', 'bath, scalar');

% Variable 'f'
varid_f = netcdf.defVar(ncid, 'f', 'NC_DOUBLE', [dimid_xi_rho, dimid_eta_rho]);
netcdf.putAtt(ncid, varid_f, 'long_name', 'Coriolis parameter at RHO-points');
netcdf.putAtt(ncid, varid_f, 'units', 'second-1');
netcdf.putAtt(ncid, varid_f, 'coordinates', 'lon_rho lat_rho');
netcdf.putAtt(ncid, varid_f, 'field', 'Coriolis, scalar');

% Variable 'pm'
varid_pm = netcdf.defVar(ncid, 'pm', 'NC_DOUBLE', [dimid_xi_rho, dimid_eta_rho]);
netcdf.putAtt(ncid, varid_pm, 'long_name', 'curvilinear coordinate metric in XI');
netcdf.putAtt(ncid, varid_pm, 'units', 'meter-1');
netcdf.putAtt(ncid, varid_pm, 'field', 'pm, scalar');

% Variable 'pn'
varid_pn = netcdf.defVar(ncid, 'pn', 'NC_DOUBLE', [dimid_xi_rho, dimid_eta_rho]);
netcdf.putAtt(ncid, varid_pn, 'long_name', 'curvilinear coordinate metric in ETA');
netcdf.putAtt(ncid, varid_pn, 'units', 'meter-1');
netcdf.putAtt(ncid, varid_pn, 'field', 'pn, scalar');

% Variable 'dndx'
varid_dndx = netcdf.defVar(ncid, 'dndx', 'NC_DOUBLE', [dimid_xi_rho, dimid_eta_rho]);
netcdf.putAtt(ncid, varid_dndx, 'long_name', 'xi derivative of inverse metric factor pn');
netcdf.putAtt(ncid, varid_dndx, 'units', 'meter');
netcdf.putAtt(ncid, varid_dndx, 'field', 'dndx, scalar');

% Variable 'dmde'
varid_dmde = netcdf.defVar(ncid, 'dmde', 'NC_DOUBLE', [dimid_xi_rho, dimid_eta_rho]);
netcdf.putAtt(ncid, varid_dmde, 'long_name', 'eta derivative of inverse metric factor pm');
netcdf.putAtt(ncid, varid_dmde, 'units', 'meter');
netcdf.putAtt(ncid, varid_dmde, 'field', 'dmde, scalar');

% Variable 'x_u'
varid_x_u = netcdf.defVar(ncid, 'x_u', 'NC_DOUBLE', [dimid_xi_u, dimid_eta_u]);
netcdf.putAtt(ncid, varid_x_u, 'long_name', 'x location of U-points');
netcdf.putAtt(ncid, varid_x_u, 'units', 'meter');

% Variable 'x_rho'
varid_x_rho = netcdf.defVar(ncid, 'x_rho', 'NC_DOUBLE', [dimid_xi_rho, dimid_eta_rho]);
netcdf.putAtt(ncid, varid_x_rho, 'long_name', 'x location of RHO-points');
netcdf.putAtt(ncid, varid_x_rho, 'units', 'meter');

% Variable 'y_rho'
varid_y_rho = netcdf.defVar(ncid, 'y_rho', 'NC_DOUBLE', [dimid_xi_rho, dimid_eta_rho]);
netcdf.putAtt(ncid, varid_y_rho, 'long_name', 'y location of RHO-points');
netcdf.putAtt(ncid, varid_y_rho, 'units', 'meter');


% Variable 'x_psi'
varid_x_psi = netcdf.defVar(ncid, 'x_psi', 'NC_DOUBLE', [dimid_xi_psi, dimid_eta_psi]);
netcdf.putAtt(ncid, varid_x_psi, 'long_name', 'x location of PSI-points');
netcdf.putAtt(ncid, varid_x_psi, 'units', 'meter');

% Variable 'y_psi'
varid_y_psi = netcdf.defVar(ncid, 'y_psi', 'NC_DOUBLE', [dimid_xi_psi, dimid_eta_psi]);
netcdf.putAtt(ncid, varid_y_psi, 'long_name', 'y location of PSI-points');
netcdf.putAtt(ncid, varid_y_psi, 'units', 'meter');


% Variable 'y_u'
varid_y_u = netcdf.defVar(ncid, 'y_u', 'NC_DOUBLE', [dimid_xi_u, dimid_eta_u]);
netcdf.putAtt(ncid, varid_y_u, 'long_name', 'y location of U-points');
netcdf.putAtt(ncid, varid_y_u, 'units', 'meter');

% Variable 'x_v'
varid_x_v = netcdf.defVar(ncid, 'x_v', 'NC_DOUBLE', [dimid_xi_v, dimid_eta_v]);
netcdf.putAtt(ncid, varid_x_v, 'long_name', 'x location of V-points');
netcdf.putAtt(ncid, varid_x_v, 'units', 'meter');

% Variable 'y_v'
varid_y_v = netcdf.defVar(ncid, 'y_v', 'NC_DOUBLE', [dimid_xi_v, dimid_eta_v]);
netcdf.putAtt(ncid, varid_y_v, 'long_name', 'y location of V-points');
netcdf.putAtt(ncid, varid_y_v, 'units', 'meter');



% Variable 'mask_rho'
varid_mask_rho = netcdf.defVar(ncid, 'mask_rho', 'NC_DOUBLE', [dimid_xi_rho, dimid_eta_rho]);
netcdf.putAtt(ncid, varid_mask_rho, 'long_name', 'mask on RHO-points');
netcdf.putAtt(ncid, varid_mask_rho, 'option(0)', 'land');
netcdf.putAtt(ncid, varid_mask_rho, 'option(1)', 'water');

if(restoring)
    varid_theta_s = netcdf.defVar(ncid, 'theta_s', 'NC_DOUBLE', []);
    netcdf.putAtt(ncid, varid_theta_s, 'long_name', 'S-coordinate surface control parameter');
    netcdf.putAtt(ncid, varid_theta_s, 'units', 'nondimensional');

    varid_theta_b = netcdf.defVar(ncid, 'theta_b', 'NC_DOUBLE', []);
    netcdf.putAtt(ncid, varid_theta_b, 'long_name', 'S-coordinate bottom control parameter');
    netcdf.putAtt(ncid, varid_theta_b, 'units', 'nondimensional');

    varid_Tcline = netcdf.defVar(ncid, 'Tcline', 'NC_DOUBLE', []);
    netcdf.putAtt(ncid, varid_Tcline, 'long_name', 'S-coordinate surface/bottom layer width');
    netcdf.putAtt(ncid, varid_Tcline, 'units', 'meter');

    varid_hc = netcdf.defVar(ncid, 'hc', 'NC_DOUBLE', []);
    netcdf.putAtt(ncid, varid_hc, 'long_name', 'S-coordinate parameter, critical depth');
    netcdf.putAtt(ncid, varid_hc, 'units', 'meter');

    varid_Cs_r = netcdf.defVar(ncid, 'Cs_r', 'NC_DOUBLE', [dimid_s_rho]);
    netcdf.putAtt(ncid, varid_Cs_r, 'long_name', 'S-coordinate stretching curves at RHO-points');
    netcdf.putAtt(ncid, varid_Cs_r, 'units', 'nondimensional');
    netcdf.putAtt(ncid, varid_Cs_r, 'valid_min', -1);
    netcdf.putAtt(ncid, varid_Cs_r, 'valid_max', 0);
    netcdf.putAtt(ncid, varid_Cs_r, 'field', 'Cs_r, scalar');

    varid_Cs_w = netcdf.defVar(ncid, 'Cs_w', 'NC_DOUBLE', [dimid_s_w]);
    netcdf.putAtt(ncid, varid_Cs_w, 'long_name', 'S-coordinate stretching curves at W-points');
    netcdf.putAtt(ncid, varid_Cs_w, 'units', 'nondimensional');
    netcdf.putAtt(ncid, varid_Cs_w, 'valid_min', -1);
    netcdf.putAtt(ncid, varid_Cs_w, 'valid_max', 0);
    netcdf.putAtt(ncid, varid_Cs_w, 'field', 'Cs_w, scalar');

    varid_t_targ = netcdf.defVar(ncid, 't_targ', 'NC_DOUBLE', [dimid_xi_rho, dimid_eta_rho, dimid_s_rho]);
    netcdf.putAtt(ncid, varid_t_targ, 'long_name', 'Targeted temperature field for nudging');
    netcdf.putAtt(ncid, varid_t_targ, 'units', 'Celsius');

    varid_nudgc = netcdf.defVar(ncid, 'nudgc', 'NC_DOUBLE', [dimid_xi_rho, dimid_eta_rho, dimid_s_rho]);
    netcdf.putAtt(ncid, varid_nudgc, 'long_name', 'Nudging coefficient');
    netcdf.putAtt(ncid, varid_nudgc, 'units', 'second-1');
end


% Define global attributes
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'type', 'Gridpak file');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'VertCoordType', 'SM09');



% End NetCDF file definition
netcdf.endDef(ncid);


% Close NetCDF file
netcdf.close(ncid);

fprintf('NetCDF file created: %s\n', filename);



%% write data into the netcdf grid file

% open the empty netcdf file
finfo = ncinfo(filename);
varname = {finfo.Variables.Name};
dim = {finfo.Variables.Size};


% write data into the empty netcdf file
ivar=find(strcmp(varname,'spherical'))
varname{ivar}
clear tmp;tmp='F';
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error('Check ', [varname{ivar},' data'])
end

ivar=find(strcmp(varname,'xl'))
varname{ivar}
clear tmp;tmp=Lx*ones(dim{ivar});
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error('Check ', [varname{ivar},' data'])
end

ivar=find(strcmp(varname,'el'))
varname{ivar}
clear tmp;tmp=Ly*ones(dim{ivar});
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error('Check ', [varname{ivar},' data'])
end

ivar=find(strcmp(varname,'f'))
varname{ivar}
clear tmp;tmp=f_rho;
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error('Check ', [varname{ivar},' data'])
end

ivar=find(strcmp(varname,'h'))
varname{ivar}
clear tmp;tmp=-h_rho; % h should be positive in roms
tmp(tmp==0)=100; % h=0 causes nan in roms ??
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error('Check ', [varname{ivar},' data'])
end


list=fieldnames(G);
for i=1:length(list)
    ivar=find(strcmp(varname,list{i}))
    varname{ivar}
    clear tmp*;eval(['tmp=G.',list{i},';']);
    ncwrite(filename,varname{ivar},tmp);
    tmp2=ncread(filename,varname{ivar})-tmp;
    if length(find(tmp2~=0))>=1
        error('Check ', [varname{ivar},' data'])
    end
end

listmask={'mask_rho'};
for i=1:length(listmask)
    ivar=find(strcmp(varname,listmask{i}))
    varname{ivar}
    clear tmp*;eval(['tmp=',listmask{i},';']);
    ncwrite(filename,varname{ivar},tmp);
    tmp2=ncread(filename,varname{ivar})-tmp;
    if length(find(tmp2~=0))>=1
        error('Check ', [varname{ivar},' data'])
    end
end

if(restoring)

    ivar=find(strcmp(varname,'Tcline'))
    varname{ivar}
    clear tmp;tmp=hc;
    ncwrite(filename,varname{ivar},tmp);
    tmp2=ncread(filename,varname{ivar})-tmp;
    if length(find(tmp2~=0))>=1
        error('Check ', [varname{ivar},' data'])
    end

    %
    ivar=find(strcmp(varname,'theta_b'))
    varname{ivar}
    clear tmp;tmp=theta_b;
    ncwrite(filename,varname{ivar},tmp);
    tmp2=ncread(filename,varname{ivar})-tmp;
    if length(find(tmp2~=0))>=1
        error('Check ', [varname{ivar},' data'])
    end

    %
    ivar=find(strcmp(varname,'theta_s'))
    varname{ivar}
    clear tmp;tmp=theta_s;
    ncwrite(filename,varname{ivar},tmp);
    tmp2=ncread(filename,varname{ivar})-tmp;
    if length(find(tmp2~=0))>=1
        error('Check ', [varname{ivar},' data'])
    end

    ivar=find(strcmp(varname,'hc'))
    varname{ivar}
    clear tmp;tmp=hc;
    ncwrite(filename,varname{ivar},tmp);
    tmp2=ncread(filename,varname{ivar})-tmp;
    if length(find(tmp2~=0))>=1
        error('Check ', [varname{ivar},' data'])
    end


    ivar=find(strcmp(varname,'Cs_w'))
    varname{ivar}
    clear tmp;tmp=Cs_w';
    ncwrite(filename,varname{ivar},tmp);
    tmp2=ncread(filename,varname{ivar})-tmp;
    if length(find(tmp2~=0))>=1
        error('Check ', [varname{ivar},' data'])
    end

    ivar=find(strcmp(varname,'Cs_r'))
    varname{ivar}
    clear tmp;tmp=Cs_r';
    ncwrite(filename,varname{ivar},tmp);
    tmp2=ncread(filename,varname{ivar})-tmp;
    if length(find(tmp2~=0))>=1
        error('Check ', [varname{ivar},' data'])
    end


    ivar=find(strcmp(varname,'t_targ'))
    varname{ivar}
    clear tmp;tmp=t_targ;
    ncwrite(filename,varname{ivar},tmp);
    tmp2=ncread(filename,varname{ivar})-tmp;
    if length(find(tmp2~=0))>=1
        error('Check ', [varname{ivar},' data'])
    end

    ivar=find(strcmp(varname,'nudgc'))
    varname{ivar}
    clear tmp;tmp=nudgc;
    ncwrite(filename,varname{ivar},tmp);
    tmp2=ncread(filename,varname{ivar})-tmp;
    if length(find(tmp2~=0))>=1
        error('Check ', [varname{ivar},' data'])
    end


end




%%  Generate an empty netcdf file of forcing

% Define the filename
filename = './Neptune_input/neptune_frc.nc';

% Create NetCDF file with CLOBBER mode to overwrite any existing file
% ncid = netcdf.create(filename, 'CLOBBER');
ncid = netcdf.create(filename, 'NETCDF4');

% Define dimensions
% dimid_frc_time = netcdf.defDim(ncid, 'frc_time', netcdf.getConstant('NC_UNLIMITED'));
dimid_frc_time = netcdf.defDim(ncid, 'frc_time', Nfrc);
dimid_xi_u = netcdf.defDim(ncid, 'xi_u', Nx-1);
dimid_eta_rho = netcdf.defDim(ncid, 'eta_rho', Ny);
dimid_xi_rho = netcdf.defDim(ncid, 'xi_rho', Nx);
dimid_eta_v = netcdf.defDim(ncid, 'eta_v', Ny-1);

% Define variables with their dimensions and attributes
varid_frc_time = netcdf.defVar(ncid, 'frc_time', 'NC_FLOAT', dimid_frc_time);
netcdf.putAtt(ncid, varid_frc_time, 'long_name', 'surface forcing time');
netcdf.putAtt(ncid, varid_frc_time, 'units', 'day');
netcdf.putAtt(ncid, varid_frc_time, 'cycle_length', Nfrc);

varid_sustr = netcdf.defVar(ncid, 'sustr', 'NC_FLOAT', [dimid_xi_u, dimid_eta_rho, dimid_frc_time]);
netcdf.putAtt(ncid, varid_sustr, 'long_name', 'surface u-momentum stress');
netcdf.putAtt(ncid, varid_sustr, 'units', 'Newton meter-2');

varid_svstr = netcdf.defVar(ncid, 'svstr', 'NC_FLOAT', [dimid_xi_rho, dimid_eta_v, dimid_frc_time]);
netcdf.putAtt(ncid, varid_svstr, 'long_name', 'surface v-momentum stress');
netcdf.putAtt(ncid, varid_svstr, 'units', 'Newton meter-2');

varid_shflux = netcdf.defVar(ncid, 'shflux', 'NC_FLOAT', [dimid_xi_rho, dimid_eta_rho, dimid_frc_time]);
netcdf.putAtt(ncid, varid_shflux, 'long_name', 'surface net heat flux');
netcdf.putAtt(ncid, varid_shflux, 'units', 'Watts meter-2');

varid_swflux = netcdf.defVar(ncid, 'swflux', 'NC_FLOAT', [dimid_xi_rho, dimid_eta_rho, dimid_frc_time]);
netcdf.putAtt(ncid, varid_swflux, 'long_name', 'surface freshwater flux (E-P)');
netcdf.putAtt(ncid, varid_swflux, 'units', 'centimeter day-1');
netcdf.putAtt(ncid, varid_swflux, 'positive', 'net evaporation');
netcdf.putAtt(ncid, varid_swflux, 'negative', 'net precipitation');


varid_SST = netcdf.defVar(ncid, 'SST', 'NC_FLOAT', [dimid_xi_rho, dimid_eta_rho, dimid_frc_time]);
netcdf.putAtt(ncid, varid_SST, 'long_name', 'sea surface temperature');
netcdf.putAtt(ncid, varid_SST, 'units', 'Celsius');

varid_SSS = netcdf.defVar(ncid, 'SSS', 'NC_FLOAT', [dimid_xi_rho, dimid_eta_rho, dimid_frc_time]);
netcdf.putAtt(ncid, varid_SSS, 'long_name', 'sea surface salinity');
netcdf.putAtt(ncid, varid_SSS, 'units', 'PSU');

varid_dQdSST = netcdf.defVar(ncid, 'dQdSST', 'NC_FLOAT', [dimid_xi_rho, dimid_eta_rho, dimid_frc_time]);
netcdf.putAtt(ncid, varid_dQdSST, 'long_name', 'sea heat flux sensitivity to SST');
netcdf.putAtt(ncid, varid_dQdSST, 'units', 'Watts meter-2 Celsius-1');

varid_swrad = netcdf.defVar(ncid, 'swrad', 'NC_FLOAT', [dimid_xi_rho, dimid_eta_rho, dimid_frc_time]);
netcdf.putAtt(ncid, varid_swrad, 'long_name', 'solar shortwave radiation');
netcdf.putAtt(ncid, varid_swrad, 'units', 'Watts meter-2');
netcdf.putAtt(ncid, varid_swrad, 'positive', 'downward flux, heating');
netcdf.putAtt(ncid, varid_swrad, 'negative', 'upward flux, cooling');

% End NetCDF file definition
netcdf.endDef(ncid);

% Close NetCDF file
netcdf.close(ncid);

disp('NetCDF file created successfully.');



%% write random forcings into the netcdf forcing file

finfo = ncinfo(filename);
varname = {finfo.Variables.Name};
dim = {finfo.Variables.Size};
ncdisp(filename)


ivar=find(strcmp(varname,'frc_time'))
varname{ivar}
clear tmp;tmp=[0:1:Nfrc-1]';
tmp = single(tmp);
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error('Check ', [varname{ivar},' data'])
end


ivar=find(strcmp(varname,'sustr'))
varname{ivar}
clear tmp;
if(spindown)
    tmp=zeros(dim{ivar});
else
    tmp=taux_u;
end
tmp = single(tmp);
ncwrite(filename,varname{ivar},tmp);
% tmp2=ncread(filename,varname{ivar})-tmp;
% if length(find(tmp2~=0))>=1
%     error('Check ', [varname{ivar},' data'])
% end


ivar=find(strcmp(varname,'svstr'))
varname{ivar}
clear tmp;
if(spindown)
    tmp=zeros(dim{ivar});
else
    tmp=tauy_v;
end
tmp = single(tmp);
ncwrite(filename,varname{ivar},tmp);
% tmp2=ncread(filename,varname{ivar})-tmp;
% if length(find(tmp2~=0))>=1
%     error('Check ', [varname{ivar},' data'])
% end


listfrc={'shflux','swflux','SST','SSS','dQdSST','swrad'};

for i=1:length(listfrc)
    ivar=find(strcmp(varname,listfrc{i}))
    varname{ivar}
    clear tmp;tmp=zeros(dim{ivar});
    tmp = single(tmp);
    ncwrite(filename,varname{ivar},tmp);
    % tmp2=ncread(filename,varname{ivar})-tmp;
    % if length(find(tmp2~=0))>=1
    %     error('Check ', [varname{ivar},' data'])
    % end
end

% ncdispmore(filename)
%% %% Generate an empty netcdf file of inital conditions

% Define the filename
filename = './Neptune_input/neptune_init.nc';

% Create a new NetCDF file, overwriting any existing file
ncid = netcdf.create(filename, 'CLOBBER');

% Define dimensions
dimid_xi_rho = netcdf.defDim(ncid, 'xi_rho', Nx);
dimid_xi_u = netcdf.defDim(ncid, 'xi_u', Nx-1);
dimid_eta_rho = netcdf.defDim(ncid, 'eta_rho', Ny);
dimid_eta_v = netcdf.defDim(ncid, 'eta_v', Ny-1);
dimid_s_rho = netcdf.defDim(ncid, 's_rho', Nlayer);
dimid_s_w = netcdf.defDim(ncid, 's_w', Nlayer+1);

% Define variables and their attributes

varid_ocean_time = netcdf.defVar(ncid, 'ocean_time', 'NC_DOUBLE', []);
netcdf.putAtt(ncid, varid_ocean_time, 'long_name', 'Time since initialization');
netcdf.putAtt(ncid, varid_ocean_time, 'units', 'second');

varid_theta_s = netcdf.defVar(ncid, 'theta_s', 'NC_DOUBLE', []);
netcdf.putAtt(ncid, varid_theta_s, 'long_name', 'S-coordinate surface control parameter');
netcdf.putAtt(ncid, varid_theta_s, 'units', 'nondimensional');

varid_theta_b = netcdf.defVar(ncid, 'theta_b', 'NC_DOUBLE', []);
netcdf.putAtt(ncid, varid_theta_b, 'long_name', 'S-coordinate bottom control parameter');
netcdf.putAtt(ncid, varid_theta_b, 'units', 'nondimensional');

varid_Tcline = netcdf.defVar(ncid, 'Tcline', 'NC_DOUBLE', []);
netcdf.putAtt(ncid, varid_Tcline, 'long_name', 'S-coordinate surface/bottom layer width');
netcdf.putAtt(ncid, varid_Tcline, 'units', 'meter');

varid_hc = netcdf.defVar(ncid, 'hc', 'NC_DOUBLE', []);
netcdf.putAtt(ncid, varid_hc, 'long_name', 'S-coordinate parameter, critical depth');
netcdf.putAtt(ncid, varid_hc, 'units', 'meter');

varid_Cs_r = netcdf.defVar(ncid, 'Cs_r', 'NC_DOUBLE', [dimid_s_rho]);
netcdf.putAtt(ncid, varid_Cs_r, 'long_name', 'S-coordinate stretching curves at RHO-points');
netcdf.putAtt(ncid, varid_Cs_r, 'units', 'nondimensional');
netcdf.putAtt(ncid, varid_Cs_r, 'valid_min', -1);
netcdf.putAtt(ncid, varid_Cs_r, 'valid_max', 0);
netcdf.putAtt(ncid, varid_Cs_r, 'field', 'Cs_r, scalar');

varid_Cs_w = netcdf.defVar(ncid, 'Cs_w', 'NC_DOUBLE', [dimid_s_w]);
netcdf.putAtt(ncid, varid_Cs_w, 'long_name', 'S-coordinate stretching curves at W-points');
netcdf.putAtt(ncid, varid_Cs_w, 'units', 'nondimensional');
netcdf.putAtt(ncid, varid_Cs_w, 'valid_min', -1);
netcdf.putAtt(ncid, varid_Cs_w, 'valid_max', 0);
netcdf.putAtt(ncid, varid_Cs_w, 'field', 'Cs_w, scalar');

varid_zeta = netcdf.defVar(ncid, 'zeta', 'NC_FLOAT', [dimid_xi_rho,dimid_eta_rho]);
netcdf.putAtt(ncid, varid_zeta, 'long_name', 'free-surface elevation');
netcdf.putAtt(ncid, varid_zeta, 'units', 'meter');
netcdf.putAtt(ncid, varid_zeta, '_FillValue', single(1.e+33));

% Variable ubar
varid_ubar = netcdf.defVar(ncid, 'ubar', 'NC_FLOAT', [dimid_xi_u, dimid_eta_rho]);
netcdf.putAtt(ncid, varid_ubar, 'long_name', 'barotropic XI-velocity');
netcdf.putAtt(ncid, varid_ubar, 'units', 'meter/second');
netcdf.putAtt(ncid, varid_ubar, '_FillValue', single(1.e+33));

% Variable vbar
varid_vbar = netcdf.defVar(ncid, 'vbar', 'NC_FLOAT', [dimid_xi_rho, dimid_eta_v]);
netcdf.putAtt(ncid, varid_vbar, 'long_name', 'barotropic ETA-velocity');
netcdf.putAtt(ncid, varid_vbar, 'units', 'meter/second');
netcdf.putAtt(ncid, varid_vbar, '_FillValue', single(1.e+33));

% Variable u
varid_u = netcdf.defVar(ncid, 'u', 'NC_FLOAT', [dimid_xi_u, dimid_eta_rho, dimid_s_rho]);
netcdf.putAtt(ncid, varid_u, 'long_name', 'XI-velocity component');
netcdf.putAtt(ncid, varid_u, 'units', 'meter/second');
netcdf.putAtt(ncid, varid_u, '_FillValue', single(1.e+33));

% Variable v
varid_v = netcdf.defVar(ncid, 'v', 'NC_FLOAT', [dimid_xi_rho, dimid_eta_v, dimid_s_rho]);
netcdf.putAtt(ncid, varid_v, 'long_name', 'ETA-velocity component');
netcdf.putAtt(ncid, varid_v, 'units', 'meter/second');
netcdf.putAtt(ncid, varid_v, '_FillValue', single(1.e+33));

% Variable temp
varid_temp = netcdf.defVar(ncid, 'temp', 'NC_FLOAT', [dimid_xi_rho, dimid_eta_rho, dimid_s_rho]);
netcdf.putAtt(ncid, varid_temp, 'long_name', 'potential temperature');
netcdf.putAtt(ncid, varid_temp, 'units', 'Celsius');
netcdf.putAtt(ncid, varid_temp, '_FillValue', single(1.e+33));

% Variable salt
varid_salt = netcdf.defVar(ncid, 'salt', 'NC_FLOAT', [dimid_xi_rho, dimid_eta_rho, dimid_s_rho]);
netcdf.putAtt(ncid, varid_salt, 'long_name', 'salinity');
netcdf.putAtt(ncid, varid_salt, 'units', 'PSU');
netcdf.putAtt(ncid, varid_salt, '_FillValue', single(1.e+33));

% End NetCDF file definition
netcdf.endDef(ncid);

% Closing the NetCDF file after writing is crucial to ensure data integrity
netcdf.close(ncid);

fprintf('NetCDF file created: %s\n', filename);


%% write the initial conditions into the empty nc file

if(use_wind)
    %%% Random noise amplitude
    tNoise = 0.1;
    %%% Add some random noise
    T_init = T_init + tNoise*(2*rand(Nx,Ny,Nlayer)-1);
end


finfo = ncinfo(filename);
varname = {finfo.Variables.Name};
dim = {finfo.Variables.Size};

ivar=find(strcmp(varname,'temp'))
varname{ivar}
clear tmp;tmp=single(T_init);
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end

ivar=find(strcmp(varname,'salt'))
varname{ivar}
clear tmp;tmp=0*T_init;
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end

ivar=find(strcmp(varname,'u'))
varname{ivar}
clear tmp;tmp=zeros(dim{ivar});
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end

ivar=find(strcmp(varname,'v'))
varname{ivar}
clear tmp;tmp=zeros(dim{ivar});
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end


ivar=find(strcmp(varname,'hc'))
varname{ivar}
clear tmp;tmp=hc;
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end
%
ivar=find(strcmp(varname,'Tcline'))
varname{ivar}
clear tmp;tmp=hc;
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end
%
ivar=find(strcmp(varname,'theta_b'))
varname{ivar}
clear tmp;tmp=theta_b;
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end
%
ivar=find(strcmp(varname,'theta_s'))
varname{ivar}
clear tmp;tmp=theta_s;
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end


ivar=find(strcmp(varname,'ubar'))
varname{ivar}
clear tmp;tmp=0*ones(dim{ivar});
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end

ivar=find(strcmp(varname,'vbar'))
varname{ivar}
clear tmp;tmp=0*ones(dim{ivar});
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end


ivar=find(strcmp(varname,'zeta'))
varname{ivar}
clear tmp;tmp=0*ones(dim{ivar});
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end

ivar=find(strcmp(varname,'ocean_time'))
varname{ivar}
clear tmp;tmp=0;
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end

ivar=find(strcmp(varname,'Cs_w'))
varname{ivar}
clear tmp;tmp=Cs_w';
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end

ivar=find(strcmp(varname,'Cs_r'))
varname{ivar}
clear tmp;tmp=Cs_r';
ncwrite(filename,varname{ivar},tmp);
tmp2=ncread(filename,varname{ivar})-tmp;
if length(find(tmp2~=0))>=1
    error(['Check ', varname{ivar},' data'])
end

% ncdispmore(filename)


%%  Check input files

checkNaN('.\Neptune_input\neptune_init.nc')
checkNaN('.\Neptune_input\neptune_grid.nc')
checkNaN('.\Neptune_input\neptune_frc.nc')

function checkNaN(ncfile)
% Open the NetCDF file
ncid = netcdf.open(ncfile, 'NC_NOWRITE');

% Get information about the file
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);

ind =0;
% Loop through each variable to check for NaNs
for varid = 0:numvars-1
    % Get variable name and other properties
    [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid, varid);

    % Read data from the current variable
    data = netcdf.getVar(ncid, varid);

    % Check if there are any NaN values
    if any(isnan(data(:)))
        fprintf('!!!!!!!!!!!!!!! Variable "%s" contains NaN values.\n', varname);
        ind = ind + 1
    else
        fprintf('Variable "%s" does not contain NaN values.\n', varname);
    end
end

if ind~=0
    fprintf('Some varibles in "%s" contain NaN values.\n',ncfile);
else
    fprintf('All varibles in "%s" do not contain NaN values.\n',ncfile);
end
% Close the NetCDF file
netcdf.close(ncid);
end


