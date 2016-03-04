function [] = fwhm_lambda()

a=load('../output/fwhm_lambda/power_ugri.txt');
binsize = 0.05;
bin = -0.6:binsize:0;
nCamcol = size(a,2)-1;
figure(1);clf;
for camcol=1:nCamcol
    subplot(3,2,camcol);
    hist(a(:,camcol+1),bin);
    title(sprintf('camcol=%d, ave=%5.2f',camcol,mean(a(:,camcol+1))));
    xlim([bin(1)-binsize bin(end)+binsize]);
    y1=ylim;
    line([-0.3 -0.3],[y1(1) y1(2)],'color','r','linewidth',2);
    line([-0.2 -0.2],[y1(1) y1(2)],'color','r','linewidth',2);
    grid;
end

end
