function [image,imSize,fs,epsilon] = loadRIRs(room,u)
% LOADRIRS Load RIR reference images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IN    room:       string room choice ('munin','freja','balder')
%       u:          under-sampling factor (3, 5)
% -----------------------------------------------------------------------------------
% OUT   image:      RIR image
%       imSize:     2-array element with image size
%       fs:         temporal sampling frequency (Hz)
%       epsilon:    noise level
% -----------------------------------------------------------------------------------
% AUTHOR: E. Zea (zea@kth.se)
% DATE: 2019-03-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch u % hard-wired noise levels (estimated with ITA Toolbox...)
    case 3 
        allEps = [8.7e-6; 7.1e-4; 9.4e-6]; % munin; freja; balder
    case 5
        allEps = [8.5e-6; 7.0e-4; 9.1e-6]; % munin; freja; balder
    
end
if strcmpi(room,'Munin')
    epsilon = allEps(1);
elseif strcmpi(room,'Freja')
    epsilon = allEps(2);
elseif strcmpi(room,'Balder')
    epsilon = allEps(3);
else 
    warning('Are you sure this is a room?');
end
% load RIR image
load(['measurementData/' room 'RIR.mat'],'out');
% generate outputs
image  = out.image;
imSize = [out.T,out.M];
fs = out.fs;
end