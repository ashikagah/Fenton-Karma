%% Fenton-Karma Model - parameter set 06 - mech 4

tic

clear all
close all
sympref('HeavisideAtOrigin',1); % Set heaviside(zero) = 1

%% Model geometry
ncols = 500; 
nrows = 500; 
h = 0.025; % grid spacing = 0.025 cm; 12.5 x 12.5 cm lattice (Figure 17)
h2 = h^2;

%% State variables

V = zeros(nrows,ncols); 
u = ones(nrows,ncols);
w = ones(nrows,ncols);

%% Constants

% Parameter set 6b in Fenton 2002 paper
Cm = 1; % Capacitance (uF_per_cm2)
% V_0 = -85; % Baseline potetial (mV)
% V_fi = 15; % Peak potential (mV)
V_c = 0.13; % (dimensionless)
V_v = 0.04; % (dimensionless)
tau_d = 0.388; % fast_inward_current (ms)
tau_v1_minus = 9; % fast_inward_current_v_gate (ms)
tau_v2_minus = 8; % fast_inward_current_v_gate (ms)
tau_v_plus = 3.33; % fast_inward_current_v_gate (ms)
tau_0 = 9; % slow_outward_current (ms)
tau_r = 33.33; % slow_outward_current (ms)
tau_si = 29; % slow_inward_current (ms)
V_csi = 0.50;  % slow_inward_current (dimensionless)
k = 15; % slow_inward_current (dimensionless)
tau_w_minus = 60; % slow_inward_current_w_gate (ms)
tau_w_plus = 250; % slow_inward_current_w_gate (ms)
IstimAmplitude = -0.2;
Dx = 0.001; % diffusivity (cm^2/ms)
Dy = 0.001; % diffusivity (cm^2/ms)

%% Integration Parameters
t_end = 100000;                 % Total time = 100,000 ms = 100 sec
dt = 0.1;                       % Duration of each time step = 0.1 ms 
T = 0:dt:t_end; T(end) = [];    % time vector (ms)
si = 10/dt;                     % Final sampling interval; 1/dt/10 = 1 kHz

% Input Signal
Iex = IstimAmplitude;

% Set initial stim current and pattern
iex = zeros(nrows,ncols);
iex(end-10:end,:) = IstimAmplitude;
iex(:,end-10:end) = IstimAmplitude;

% Final matrices
Vts = zeros(ncols,nrows,floor(numel(T)/si));
% uts = zeros(ncols,nrows,floor(numel(T)/si));
% wts = zeros(ncols,nrows,floor(numel(T)/si));
s=1;

% % Setup image
% V0 = ones(nrows,ncols);
% ih = imagesc(V0); colorbar; caxis([0 1]);
% colormap(lce); axis image off; th=title('');
% set(gcf,'position',[500 600 512 512],'color',[1 1 1]);

Ts = 5620;

n1e=2/dt;       % Step at which to end 1st stimulus
n2b=Ts;         % Step at which to begin 2nd stimulus
n2e=Ts+2/dt;    % Step at which to end 2nd stimulus

for t = 1:numel(T)
     %fprintf('Working on frame %d sec...\n',T(t)/1000);
%%%%%%%%%%%%%%%%% Cross-field stimulation %%%%%%%%%%%%%%%%%
    if t == n1e; iex=zeros(nrows,ncols); end% End 1st stimulus
    if t == n2b; iex(nrows/2:nrows/2+10,1:end/2)=Iex; end % Begin 2nd stimulus
    if t == n2e; iex=zeros(nrows,ncols); end % End 2nd stimulus
    
%%%%%%%%%%%%%%%%% FK Main %%%%%%%%%%%%%%%%%
    % Create padded v matrix to incorporate Neuman boundary conditions 
    VV=[[0 V(2,:) 0];[V(:,2) V V(:,end-1)];[0 V(end-1,:) 0]];
 
    %% Algebraic variables
    p = heaviside(V - V_c); % (dimensionless)
    q = heaviside(V - V_v); % q (dimensionless)
    J_fi = - u.*p.*(1 - V).*(V - V_c)/tau_d; % fast_inward_current (per_ms)
    J_so = V.*(1 - p)/tau_0 + p/tau_r; % slow_outward_current (per_ms)
    tau_v_minus =  q * tau_v1_minus + (1 - q) * tau_v2_minus; % fast_inward_current_v_gate (ms)
    %tau_v_minus = 18.2;
    J_si  = - w.*(1 + tanh(k*(V - V_csi)))/(2*tau_si); % slow_inward_current (per_ms)

    %% Rate variables
 
    % du/dt (dimensionless)
    Vxx = (VV(2:end-1,1:end-2) + VV(2:end-1,3:end) -2*V)/h2; 
    Vyy = (VV(1:end-2,2:end-1) + VV(3:end,2:end-1) -2*V)/h2;    
    Vdot =  - (J_fi + J_so + J_si + iex) + Dx * Vxx + Dy * Vyy;
    V_new = V + Vdot * dt;

    % dv/dt v fast_inward_current_v_gate (dimensionless)
    udot = ((1 - p).*(1 - u)./tau_v_minus) - (p.*u)/tau_v_plus;
    u_new = u + udot * dt;
    
    % dw/dt w slow_inward_current_w_gate (dimensionless)
    wdot = ((1 - p).*(1 - w)./tau_w_minus) - (p.*w)/tau_w_plus;
    w_new = w + wdot * dt;
    
    % Update rate variables
    V = V_new; clear V_new;
    u = u_new; clear u_new;
    w = w_new; clear w_new;
    
%     % Update image and text 
%     % Vm = V_0 + u*(V_fi - V_0); % membrane potential (mV)
%     set(ih,'cdata',V);
%     set(th,'string',sprintf('%0.3f',T(t)/1000)); % in sec
%     drawnow
        
    % Create a matrix for downsampled time series
    if rem(t,si) == 0
        l = floor(t/si);
        Vts(:,:,l) = V;
        % uts(:,:,l) = u;
        % wts(:,:,l) = w;
        fprintf('%.02f percent completed ...\n',100*l/floor(numel(T)/si));
    end
end

Vts(:,:,1:floor(Ts/si)) = [];

% Save the time series data
save mech04_00_90.mat Vts -v7.3; % large size data (7.5GB)

% Make a movie
k = 100; % grid size - change this parameter for k=30, k=15, k=8 and k=4
[X,Y] = meshgrid(linspace(1,k,size(Vts,1)),linspace(1,k,size(Vts,2)));
[x,y] = meshgrid(1:k,1:k);
skip = 2; % compression ratio in time dimension
figure;
i0 = zeros(size(Vts(:,:,1)));
ih = imagesc(i0(:,:,1)); caxis([0 1]);
colormap(lce); axis image off; colorbar
set(gcf,'position',[500 600 512 512],'color',[1 1 1]);
k = 1;
for frame=1:500 % only first 5 seconds (size<10MB)
    V = double(Vts(:,:,frame));
    vm = interp2(X,Y,V,x,y,'cubic'); % downsample to a smaller grid size (kxk)
    if rem(frame,skip)==1
        set(ih,'cdata',vm);
        drawnow
        mov(k) = getframe;
        k = k + 1;
    end
end
writerObj = VideoWriter(['mech04_00_05.avi'],'Motion JPEG AVI');
writerObj.FrameRate = 20;
open(writerObj);
writeVideo(writerObj,mov);
close(writerObj);
close all
clear mov

%% long, high-res movie to make sure spirals persist
figure;
i0 = zeros(size(Vts(:,:,1)));
ih = imagesc(i0(:,:,1)); caxis([0 1]);
colormap(lce); axis image off; colorbar
set(gcf,'position',[500 600 512 512],'color',[1 1 1]);
for frame=1:size(Vts,3) % only first 5 seconds (size<10MB)
    set(ih,'cdata',Vts(:,:,frame));
    drawnow
    mov(frame) = getframe;
end
writerObj = VideoWriter(['mech04_00_90.avi'],'Motion JPEG AVI');
writerObj.FrameRate = 20;
open(writerObj);
writeVideo(writerObj,mov);
close(writerObj);

toc