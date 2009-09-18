% Roe nofix
errv = [0.1786, 0.0971, 0.0543, 0.0279];
errp = [0.0614, 0.0404, 0.0264, 0.0167];
errr = [0.0644, 0.0424, 0.0302, 0.0202];
% HLL
errv = [0.1894, 0.1011, 0.0554, 0.0283];
errp = [0.0616, 0.0407, 0.0263, 0.0169];
errr = [0.06, 0.0426, 0.03, 0.0201];
% HLLC
errp = [0.0601, 0.04, 0.0265, 0.0168];
errv = [0.1729, 0.0955, 0.0558, 0.0276];
errr = [0.0574, 0.0414, 0.0299, 0.0201];
% ROE fix
errp = [ 0.0680, 0.0396, 0.0247, 0.0164];
errv = [ 0.1815, 0.0941, 0.0498, 0.0234];
errr = [ 0.0626, 0.04,   0.0271, 0.0191];

estp = -log(errp(2:end)./errp(1:end-1))/log(2.0);
estv = -log(errv(2:end)./errv(1:end-1))/log(2.0);
estr = -log(errr(2:end)./errr(1:end-1))/log(2.0);
loglog(h,errp,'-sq','Color','k','LineWidth',2.0,'MarkerFaceColor','w'); hold on;
loglog(h,errv,'-sq','Color','k','LineWidth',2.0,'MarkerFaceColor','w');
loglog(h,errr,'-sq','Color','k','LineWidth',2.0,'MarkerFaceColor','w');
text(max(h)+1/128,errv(1),['Velocita'''],'HorizontalAlignment','left','BackgroundColor','w');
text(max(h)+1/128,errp(1)/1.15,['Pressione'],'HorizontalAlignment','left','BackgroundColor','w');
text(max(h)+1/128,errr(1)*1.1,['Densita'''],'HorizontalAlignment','left','BackgroundColor','w');
axis([min(h)/1.5,max(h)+1/16, min(errp)/1.5, max(errv)*1.5]);
xlabel('h')
ylabel('Errore in norma {L_1}')
grid on