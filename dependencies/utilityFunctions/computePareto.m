function [beta_star,Jcurve,beta_set] = computePareto(y_hat,LAMBDA,E,imSize,Sk)
% COMPUTEPARETO Compute L-curve corner with cubic spline interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IN    y_hat:      observation vector
%       LAMBDA:     masking matrix
%       E:          expansion matrix
%       imSize:     2-element array with image dimensions
%       Sk:         shearlet dictionary bases
% -----------------------------------------------------------------------------------
% OUT   beta_star:  optimal threshold
%       Jcurve:     L-curvature function
%       beta_set:   set of regularization parameters
% -----------------------------------------------------------------------------------
% AUTHOR: E. Zea (zea@kth.se)
% DATE: 2019-03-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nFinding optimal Pareto point...\n');
% define the set of possible reg. parameters
beta_set = sort(logspace(-1,-2.5,50),'descend');
% preallocate...
eta = zeros(numel(beta_set),1); % sparsity norm
rho = zeros(numel(beta_set),1); % residual norm
% iterate through beta_set...
tic;
alpha0 = E'*ffst(LAMBDA'*y_hat,imSize,Sk);  % initialize solution vector
for bb = 1:numel(beta_set)
    alpha1 = wthresh(alpha0 + E'*ffst( LAMBDA'*(y_hat - LAMBDA*iffst(E*alpha0,imSize,Sk)), ...
                     imSize,Sk ), 's', beta_set(bb)*norm(alpha0,Inf));
    eta(bb) = norm(alpha1,1);
    rho(bb) = norm(y_hat-LAMBDA*iffst(E*alpha1,imSize,Sk),2);
    disp(['Pareto iteration ' num2str(bb) '/' num2str(numel(beta_set)) '.']);
end
toc
% L-curve cubic spline interpolation
eta_sp   = spline(beta_set,log(eta));
rho_sp   = spline(beta_set,log(rho));
eta_der1 = fnder(eta_sp,1); eta_prime1 = ppval(eta_der1,beta_set); % eta'
eta_der2 = fnder(eta_sp,2); eta_prime2 = ppval(eta_der2,beta_set); % eta''
rho_der1 = fnder(rho_sp,1); rho_prime1 = ppval(rho_der1,beta_set); % rho'
rho_der2 = fnder(rho_sp,2); rho_prime2 = ppval(rho_der2,beta_set); % rho''
Jcurve   = ( rho_prime2.*eta_prime1 - rho_prime1.*eta_prime2 )./...
           ( rho_prime1.^2+eta_prime1.^2 ).^1.5;
% plot curvature function
figure(2);
subplot(1,2,2);
fig = semilogx(beta_set,Jcurve,'ko');
set(fig,'linewidth',2,'markersize',8);
axis tight; set(gca,'fontsize',20);
xlabel('$\beta$','interpreter','latex');
ylabel('$\mathcal{J}(\beta)$','interpreter','latex');
grid on; drawnow;
% find optimal beta
[~,idx_max_curv] = max(Jcurve);
beta_star = beta_set(idx_max_curv);
% plot L-curve
figure(2);
subplot(1,2,1);
fig = plot(log(rho),log(eta),'ko');
set(fig,'linewidth',2,'markersize',8);
axis tight; set(gca,'fontsize',20);
xlabel('$\rho(\beta) = \| \hat{\mathbf{y}} - \mathbf{\Phi}\alpha \|_2$',...
       'interpreter','latex');
ylabel('$\eta(\beta) = \| \alpha \|_1$','interpreter','latex');
vline(log(rho(idx_max_curv)));
hline(log(eta(idx_max_curv)));
grid on; drawnow;
end