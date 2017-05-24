%       ***************************************************
%       *  Copyright (C) 2017, Hiroshi Ashikaga, MD, PhD  *
%       *  hashika1@jhmi.edu                              *
%       *  Cardiac Arrhythmia Service                    *
%       *  Johns Hopkins University School of Medicine    *
%       *  Baltimore, Maryland, USA                       *
%       *  5/24/2017                                      *
%       ***************************************************

%% Generate spiral wave breakup using Fenton-Karma model of cardiac action potentials
% Fenton FH and Karma A AD. Vortex dynamics in three-dimensional continuous myocardium with 
% fiber rotation: Filament instability and fibrillationa. Chaos 8: 20-47, 1998

function ts = fk_breakup(time_units)
% INPUT:    
%   time_units   ... Total time units of time series (1ms/unit); 
%                    e.g. 10,000units = 10,000ms = 10sec
% OUTPUT:
%   ts           ... 2-D time series of excitation variable V [N x M x time]

sympref('HeavisideAtOrigin',1);             % Set heaviside(zero) = 1

% Model geometry
ncols = 500;                                % Number of columns; e.g. 500unit x 0.025mm/unit = 12.5cm
nrows = 500;                                % Number of rows; e.g. 500unit x 0.025mm/unit = 12.5cm
h = 0.025;                                  % Grid spacing [mm/unit]; e.g. 12.5 x 12.5 cm lattice
h2 = h^2;

% State variables
V = zeros(nrows,ncols);                     % V: transmembrane voltage
u = ones(nrows,ncols);                      % u: gate variable for inactivation of a fast inward current (Ifi) after depolarization and its reactivation after repolarization
w = ones(nrows,ncols);                      % w: gate variable for inactivation and reactivation of a slow inward current (Isi)

% Model parameters
Cm = 1;                                     % Capacitance (uF_per_cm2)
V_c = 0.13;                                 % (dimensionless)
V_v = 0.04;                                 % (dimensionless)
tau_d = 0.395;                              % Fast_inward_current (ms)
tau_v1_minus = 9;                           % Fast_inward_current_v_gate (ms)
tau_v2_minus = 8;                           % Fast_inward_current_v_gate (ms)
tau_v_plus = 3.33;                          % Fast_inward_current_v_gate (ms)
tau_0 = 9;                                  % Slow_outward_current (ms)
tau_r = 33.33;                              % Slow_outward_current (ms)
tau_si = 29;                                % Slow_inward_current (ms)
V_csi = 0.50;                               % Slow_inward_current (dimensionless)
k = 15;                                     % Slow_inward_current (dimensionless)
tau_w_minus = 60;                           % Slow_inward_current_w_gate (ms)
tau_w_plus = 250;                           % Slow_inward_current_w_gate (ms)
Dx = 0.001;                                 % X diffusivity (cm^2/ms)
Dy = 0.001;                                 % Y diffusivity (cm^2/ms)

% Integration Parameters
dt = 0.1;                                   % Duration of each time step = 0.1 ms 
T = 0:dt:time_units; T(end) = [];           % Time vector (ms)
si = 10/dt;                                 % Final sampling interval; 10/0.1 = 100time steps = 100 x 0.1ms = 10ms/frame
                                            % Final sampling rate = 1,000ms/10ms ~ 100Hz

% External current parameters
Iex = -0.2;                                 % Stimulation amplitude
iex = zeros(nrows,ncols);
iex(end-10:end,:) = Iex;
iex(:,end-10:end) = Iex;

% Output matrix
ts = zeros(ncols,nrows,floor(numel(T)/si));

% % Setup image
% V0 = ones(nrows,ncols);
% ih = imagesc(V0); colorbar; caxis([0 1]);
% colormap(lce); axis image off; th=title('');
% set(gcf,'position',[500 600 512 512],'color',[1 1 1]);

% s1-s2 stimulation to induce spiral waves
Ts = 3400;
n1e = 2/dt;                                 % Step at which to end s1
n2b = Ts;                                   % Step at which to begin s2
n2e = Ts+2/dt;                              % Step at which to end s2

% Fenton-Karma main
s = 1;

for t = 1:numel(T)

%%%%%%%%%%%%%%%%%    s1-s2 stimulation    %%%%%%%%%%%%%%%%%

    if t == n1e; iex=zeros(nrows,ncols); end % s1 ends
    if t == n2b; iex(nrows/2:nrows/2+10,1:end/2)=Iex; end % s2 begins
    if t == n2e; iex=zeros(nrows,ncols); end % s2 ends
    
%%%%%%%%%%%%%%%%%         FK Main         %%%%%%%%%%%%%%%%%

    % Padded v matrix for Neuman boundary conditions 
    VV=[[0 V(2,:) 0];[V(:,2) V V(:,end-1)];[0 V(end-1,:) 0]];
 
    % Algebraic variables
    p = heaviside(V - V_c);                 % (dimensionless)
    q = heaviside(V - V_v);                 % q (dimensionless)
    J_fi = - u.*p.*(1 - V).*(V - V_c)/tau_d; % fast_inward_current (per_ms)
    J_so = V.*(1 - p)/tau_0 + p/tau_r;      % slow_outward_current (per_ms)
    tau_v_minus =  q * tau_v1_minus + (1 - q) * tau_v2_minus; % fast_inward_current_v_gate (ms)
    J_si  = - w.*(1 + tanh(k*(V - V_csi)))/(2*tau_si); % slow_inward_current (per_ms)
 
    % dV/dt (dimensionless)
    Vxx = (VV(2:end-1,1:end-2) + VV(2:end-1,3:end) -2*V)/h2; 
    Vyy = (VV(1:end-2,2:end-1) + VV(3:end,2:end-1) -2*V)/h2;    
    Vdot =  - (J_fi + J_so + J_si + iex) + Dx * Vxx + Dy * Vyy;
    V_new = V + Vdot * dt;                  % Explicit Euler

    % du/dt u fast_inward_current_u_gate (dimensionless)
    udot = ((1 - p).*(1 - u)./tau_v_minus) - (p.*u)/tau_v_plus;
    u_new = u + udot * dt;                  % Explicit Euler
    
    % dw/dt w slow_inward_current_w_gate (dimensionless)
    wdot = ((1 - p).*(1 - w)./tau_w_minus) - (p.*w)/tau_w_plus;
    w_new = w + wdot * dt;                  % Explicit Euler
    
    % Update rate variables
    V = V_new; clear V_new;
    u = u_new; clear u_new;
    w = w_new; clear w_new;
   
%     % Update image and text 
%     set(ih,'cdata',V);
%     set(th,'string',sprintf('%0.3f',T(t)/1000)); % in sec
%     drawnow
    
    % Downsample to create output matrix
    if rem(t,si) == 0
        l = floor(t/si);
        ts(:,:,l) = V;
        fprintf('%1.0f percent completed ...\n',100*l/floor(numel(T)/si));
    end
end

ts(:,:,1:floor(Ts/si)) = [];