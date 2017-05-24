%       ***************************************************
%       *  Copyright (C) 2017, Hiroshi Ashikaga, MD, PhD  *
%       *  hashika1@jhmi.edu                              *
%       *  Cardiac Arrhythmia Service                     *
%       *  Johns Hopkins University School of Medicine    *
%       *  Baltimore, Maryland, USA                       *
%       *  5/21/2017                                      *
%       ***************************************************

%% Generate phase map of spiral waves

function p = phase_map(ts)
% INPUT:    
%   ts          ... 2-D time series ts [N x M x time]
%
% OUTPUT:
%   p           ... 2-D time series of phase [N x M x time]

offset = 0.4;               % Subtract 0.4 (half the max amplitude) to intentionally 
                            % create negative values to make the Hilbert transform phase [-pi pi]
H = zeros(size(ts));        % Final matrix size = 120x120x4500 (~500MB)
p = zeros(size(ts));

% Hilbert transform of excitation variable   
for m=1:size(ts,1)
    for n=1:size(ts,2)
        H(m,n,:) = hilbert(squeeze(ts(m,n,:))-offset);
    end
    fprintf('%1.0f percent completed ...\n',100*m/size(ts,1));
end

% Calculate phase at each time frame
for frame=1:size(H,3)
    p(:,:,frame) = angle(H(:,:,frame));
end
