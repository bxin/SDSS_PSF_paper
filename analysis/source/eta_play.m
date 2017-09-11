function [] = eta_play(writepdf)

L = importdata('../data/Stripe82RunList.dat',' ',12);
nrun = size(L.data,1);

y = zeros(5,6,nrun);
ye = zeros(5,6,nrun);

if writepdf
    outputfile=['../output/eta_change_' datestr(now,'yymmdd')];
    outputps = [outputfile '.ps'];
    outputpdf = [outputfile '.pdf'];
end

for irun=1:nrun
    dfile = sprintf('../SDSSdata/masterTXT/run%d.txt',L.data(irun,1));
    M = importdata(dfile, ' ', 1);
    
    %field 	 camCol  filter FWHMvK eta psf_width airmass mjd		 psf_nstar  neff_psf	 sky_frames
    clf;
    
    for iband = 0:4
        for icamcol = 1:6
            subplot(5,6,iband*6+icamcol);
            idxsub = (M.data(:,3)==iband & M.data(:,2)==icamcol); % & M.data(:,1)==0);
            myd = M.data(idxsub, 5);
            hist(myd);
            y(iband+1,icamcol,irun)=mean(myd);
            ye(iband+1,icamcol,irun)=std(myd);
            title(sprintf('%.2f +/- %.2f', y(iband+1,icamcol,irun),ye(iband+1,icamcol,irun)));
            % xlim([0 4]);
            if (iband==4 && icamcol == 6)
                pageTitle = sprintf('run%d', L.data(irun,1));
                text(.65,.85,pageTitle,'fontsize',010, 'Units','normalized');
            end
        end
    end
    if writepdf
        set(gcf,'PaperPositionMode', 'manual', 'PaperUnits','centimeters', 'Paperposition',[1 1 20 25])
        if irun==1
            print(outputps,'-dpsc');
        else
            print(outputps,'-append','-dpsc');
        end
    end
end

for iband = 0:4
    for icamcol = 1:6
        subplot(5,6,iband*6+icamcol);
        errorbar(1:nrun, y(iband+1,icamcol,:),ye(iband+1,icamcol,:));
        xlim([0 nrun+1]);
        ylim([0 3]);
    end
end
if writepdf
    print(outputps,'-append','-dpsc');
end

if writepdf
    syscmd = ['/usr/local/bin/ps2pdf ' outputps ' ' outputpdf];
    system(syscmd);
    syscmd = ['rm ' outputps];
    system(syscmd);
    % syscmd = ['open -a /Applications/Adobe\ Acrobat\ X\ Pro/Adobe\ Acrobat\ Pro.app ' outputpdf];
    syscmd = ['open ' outputpdf];
    system(syscmd);
end

end
