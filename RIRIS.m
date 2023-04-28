function RIRIS(room,u,tau,nu_max,saveFlag)
%RIRIS Main script to interpolate RIRs with shearlet dictionaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       This script exectutes the interpolation of RIR images measured in the rooms 
%  investigated in [1], using a dictionary of cone-adapted shearlets [2], 
%   via the iterative soft-thresholding algorithm (ISTA) [3]. 
% -----------------------------------------------------------------------------------
% IN    room:       string room choice ('munin','freja','balder')
%       u:          under-sampling factor (3, 5)
%       tau:        no. decomposition scales (2, 3, 4, 5)
%       nu_max:     max. no. ISTA iterations
%       saveFlag:   binary flag to store data (default 0)
% -----------------------------------------------------------------------------------
% OUT   void
% -----------------------------------------------------------------------------------
% DEPENDENCIES(!)
%       Folders:    'measurementData'                   measured RIRs
%                   'basisFunctions'                    shearlet dictionary
%                   'utilityFunctions'                  utility functions
%                   'regularizationData' (optional)     Pareto thresholds
%
%   The shearlet dictionary bases are obtained with the Fast Finite
%   Shearlet Transform (FFST) Matlab toolbox [2]. 
% -----------------------------------------------------------------------------------
% EXAMPLE OF USAGE: 
%   Interpolate RIRs in room 'Munin', provided under-sampling ratio of 3, a shearlet 
%   dictionary with 4 scales, and running 20 thresholding (ISTA [3]) iterations:
%
%       >> RIRIS('Munin',3,4,20);
%
%   Running this line of code in the Command Window will output the interpolated 
%   results in Figure 1, and the modal assurance criterion (MAC) versus frequency in 
%   Figure 2. The normalized mean-squared error (NMSE) and the estimated computation 
%   time (CT) are output in the Command Window. See also Fig. 16 in [1] for further 
%   details on computation times versus no. thresholding iterations.
% -----------------------------------------------------------------------------------
% Copyright 2019 Elias Zea
% RIRIS is covered by a GPL v3 license (see COPYING for license terms)
% CURRENT RELEASE: VERSION 1.0
% DATE: 2023-04-27
% -----------------------------------------------------------------------------------
% CONTACT: Elias Zea (zea@kth.se)
%          Marcus Wallenberg Laboratory for Sound and Vibration Research
%          KTH Royal Institute of Technology
%          Teknikringen 8
%          10044 Stockholm, Sweden
% -----------------------------------------------------------------------------------
% REFERENCES: 
% [1] E. Zea, 'Compressed sensing of impulse responses in rooms of unknown properties 
%     and contents', J. Sound Vib. 459, 114871 (2019)
% [2] S. Häuser, G. Steidl, 'Fast finite shearlet transform: a tutorial', ArXiv 
%     1202.1773, 1-41 (2012)
% [3] I. Daubechies, M. Defrise, C. De Mol, 'An iterative thresholding algorithm for
%     linear inverse problems with a sparsity constraint', Comm. Pure Appl. Math. 
%     57, 1413-1457 (2004)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%               BEGIN CODE...               %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add dependencies...
addpath(genpath('dependencies'));
% set default parameters (running script from Editor)
if nargin == 0 
    room = 'Balder'; % room choice: 'Munin','Freja','Balder'
    u = 3;           % under-sampling ratio (integer)
    tau = 3;         % no. decomposition scales (2,...,5)
    nu_max = 15;     % thresholding iterations (positive integer)
    saveFlag = 0;    % store data (true/false)
elseif nargin == 4   
    saveFlag = 0;
end
fprintf(['\n---------------- INTERPOLATION SETUP ----------------' ...
         '\n     >> Room: ' room...
         '\n     >> u = %1d'...
         '\n     >> tau = %1d decomposition scales'...
         '\n     >> nu_max = %3d thresholding iterations\n'],u,tau,nu_max);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%           INTERPOLATION SETUP             %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load image
[image,imSize,~,epsilon] = loadRIRs(room,u);
% mask image
[lambda,ranPos,~] = selectionMtrx(u,imSize);
image_under = lambda.*image; 
% extrapolate image
[image_under,imSize,~,LAMBDA,origAprtr] = ...
    extRIRimage(image_under,lambda,ranPos);
% matrix-vector masking of observations
y_hat = LAMBDA*image_under(:);
% expansion matrix (remove vertical-like shearlets)
E = expansionMtrx(imSize,tau);
% load basis functions... (assuming FFST toolbox is not locally installed)
load(['basisFunctions/' room '_tau' num2str(tau) '.mat'],'Psi'); Sk = Psi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%             CHOOSING THRESHOLD            %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paretoFlag = 0; % change to 1 to actually calculate beta_star
if ~paretoFlag
    load(['regularizationData/' room '_u' num2str(u) '_tau' ...
          num2str(tau) '.mat'],'reguThresh');
    beta_star = reguThresh.beta_star;
else
    [beta_star,Jcurve,beta_set] = computePareto(y_hat,LAMBDA,E,imSize,Sk);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%               ISTA RECOVERY               %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; alpha = ista(y_hat,LAMBDA,E,nu_max,imSize,Sk,beta_star,epsilon); 
CT = toc; % estimate computation time (CT)
fprintf('\n    >> Computation time (CT): %3.1f minutes\n',CT/60);
% recover inpainted image from sparse coefficients (Eq. 19)
image_recov = reshape(iffst(E*alpha,imSize,Sk),imSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%            PERFORMANCE METRICS            %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NMSE,MAC,frqMAC] = perforMetrics(image,image_recov(:,origAprtr),...
                                  image_under(:,origAprtr));
fprintf('-----------------------------------------------------\n');
fprintf('-----------------------------------------------------\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%               STORE RESULTS               %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFlag
    if paretoFlag % only if computePareto runs...
        % regularization data...
        reguThresh.beta_star = beta_star;
        reguThresh.Jcurve    = Jcurve;
        reguThresh.beta_set  = beta_set;
        save(['regularizationData/' room '_u' num2str(u) '_tau' ...
              num2str(tau) '.mat'],'reguThresh');
    end
	% interpolation results...
    if ~exist(['results/' room],'dir')
        mkdir(['results/' room]);
    end
    results.NMSE = NMSE;
    results.MAC = MAC;
    results.frqMAC = frqMAC;
    results.image = image;
    results.image_recov = image_recov(:,origAprtr);
    results.image_under = image_under(:,origAprtr);
    save(['results/' room '/u' num2str(u) 'tau' num2str(tau) '.mat'],'results');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                END CODE...                %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove dependencies...
rmpath(genpath('dependencies'));
end
