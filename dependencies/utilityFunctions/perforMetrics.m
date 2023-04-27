function [NMSE,MAC,frqMAC] = perforMetrics(image,image_recov,image_under)
% PERFORMETRICS Assessment function with performance metrics and figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IN    image:       reference RIR image
%       image_recov: interpolated RIR image
%       image_under: under-sampled RIR image (LAMBDA^T*y_hat)
% -----------------------------------------------------------------------------------
% OUT   NMSE_dB:     normalized mean-squared error (dB)
%       MAC:         modal assurance criterion
%       frqMAC:      MAC's frequency axis vector (Hz)
% -----------------------------------------------------------------------------------
% AUTHOR: E. Zea (zea@kth.se)
% DATE: 2019-03-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T  = size(image,1);
M0 = size(image,2);
y_tru = image; 
y_rec = image_recov;
dx = 3e-2;  % mic array spacing (m)
fs = 11250; % sampling frequency (Hz)
% define space-time grid
t_idxs = 1:256;
[XX,TT] = meshgrid(dx*[0:M0-1],1000*(t_idxs-1)/fs);
% plot image results
aux = image(t_idxs,:);
zscale = [min(real(aux(:))) max(real(aux(:)))];
figure(1);
subplot(1,3,1);
fig = surf(XX,TT,image_under(t_idxs,:));
set(fig,'edgecolor','none'); 
set(gca,'fontsize',18,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','interpreter','latex'); 
ylabel('$t$ (ms)','interpreter','latex'); colormap gray;
title('Under-sampled','interpreter','latex');
axis tight; view(0,90); caxis(zscale); ylim([0 15]); drawnow; 
subplot(1,3,2);
fig = surf(XX,TT,image_recov(t_idxs,:));
set(fig,'edgecolor','none'); 
set(gca,'fontsize',18,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','interpreter','latex'); 
ylabel('$t$ (ms)','interpreter','latex'); colormap gray;
title('Interpolated','interpreter','latex');
axis tight; view(0,90); caxis(zscale); ylim([0 15]); drawnow; 
subplot(1,3,3);
fig = surf(XX,TT,image(t_idxs,:));
set(fig,'edgecolor','none'); 
set(gca,'fontsize',18,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','interpreter','latex'); 
ylabel('$t$ (ms)','interpreter','latex'); colormap gray;
title('Reference','interpreter','latex');
axis tight; view(0,90); caxis(zscale); ylim([0 15]); drawnow; 
% calculate normalized mean-squared error
NMSE = 0;
for mm = 1:M0
    NMSE = NMSE + norm( y_tru(:,mm)-y_rec(:,mm) )^2/norm( y_tru(:,mm) )^2;
end
NMSE = 10*log10( NMSE/M0 );
fprintf('\n    >> Normalized mean-squared error (NMSE): %2.1f dB\n\n',NMSE);
% calculate MAC
FFT_Ptru = fft(y_tru)/T; 
FFT_Prec = fft(y_rec)/T;
fs = 11250;
frqMAC = 0:fs/T:fs/2-1;
MAC    = zeros(ceil(T/2),1);
for ff = 1:ceil(T/2)
    MAC(ff) = abs(FFT_Ptru(ff,:)*(FFT_Prec(ff,:))')^2./...
            ( abs(FFT_Ptru(ff,:)*(FFT_Ptru(ff,:))').*...
              abs(FFT_Prec(ff,:)*(FFT_Prec(ff,:))') );
end
figure(2);
fig = plot(frqMAC,MAC,'Color',[0.38,0.38,0.38]);
set(fig,'linewidth',2); 
xlabel('$f$ (Hz)','interpreter','latex'); 
ylabel('MAC','interpreter','latex');
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
ylim([0 1]); xlim([100 4500]); 
set(gca,'YGrid','on','XGrid','on','XMinorTick','on');
xticks(500:500:4500);
xticklabels({'','1000','','2000','','3000','','4000'});
end