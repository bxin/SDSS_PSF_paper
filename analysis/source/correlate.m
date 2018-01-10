function [] = correlate(writeMaster, usebands, col, filter)

% 1/7/18, this is only about alpha and FWHM. It used to be about alpha, beta, tau, and FWHM.
% we no longer use beta for spatial correlation. 
% beta was the correlation coefficients between camcols.
% now we use spatial structure function.
% we still use tau for temporal correlation, but we only fit temporal PSD
% for runs longer than 100 fields. And we do not correlate alpha and tau anymore.

% correlate(1,'ugriz',-16,2) %camcol 1-6
% correlate(1,'ugri',-16,2)
% correlate(1,'ugriz',-25,2) %camcol 2-5
% correlate(1,'ugri',-25,2)

% This writes out files like SDSSdata/ugriz/Stripe82_u_col1.txt.
% the only reason the filter is a variable is that, we need to correlate alpha with FWHM_band.

figure(1);clf;

bandL = {'u','g','r','i','z'};
% camcol = 3; % use camcol #3, b/c it is least affected by ccd edge effects(?)
% band = 2; % use r band to examime the correlations
nrun = 108;
alpha = zeros(nrun, 6);%for 6 camcol, 
alphaErr =zeros(nrun, 6);
seeing = zeros(nrun, 6);
for band=0:4
    if (filter<0 || band == filter)
        for camcol = 1:6
            if (col<0 || col == camcol)
                fname = sprintf('../SDSSdata/%s/Stripe82_%s_col%d.txt',usebands,bandL{band+1},camcol);
                if writeMaster
                    % alpha (fwhm vs lambda)
                    figure(2);subplot(1,2,1);
                    a=load(sprintf('../output/fwhm_lambda/power_%s.txt',usebands));
                    runlist = a(:,1);
                    alpha(:,camcol) = a(:,camcol*2);
                    alphaErr(:,camcol) = a(:,camcol*2+1);
                    hist(alpha(:,camcol));
                    % title(sprintf('alpha in camcol#3: mean = %.2f, std=%.2f\n', mean(alpha), std(alpha)));
                
                    % average seeing for each run
                    figure(2);subplot(1,2,2);
                    nrun = length(runlist);
                    for irun = 1:nrun
                        run = runlist(irun);
                        a = importdata(sprintf('../SDSSdata/masterTXT/run%d.txt', run),' ', 1);
                        a = a.data;
                        idx  = (a(:,2)==camcol) & (a(:,3)==band);
                        seeing(irun,camcol) = mean(a(idx,4));
                    end
                    hist(seeing(:,camcol));
                    
                    % write out master file
                    fid = fopen(fname, 'w');
                    fprintf(fid,'run \t alpha \t alphaErr \t seeing\n');
                    for irun = 1:nrun
                        fprintf(fid, '%d \t %.2f \t %.4f \t %.3f\n', runlist(irun), ...
                            alpha(irun,camcol), alphaErr(irun,camcol), seeing(irun,camcol));
                    end
                    fclose(fid);
                else
                    a = importdata(fname,' ',1);
                    a = a.data;
                    runlist = a(:,1);
                    alpha(:,camcol) = a(:,2);
                    seeing(:,camcol) = a(:,3);
                end
                
                figure(1);
                if col<0  
                    subplot(5,6,band*6+camcol)
                end
                scatter(seeing(:,camcol), alpha(:,camcol),400,'.');grid;title(sprintf('%s, col%d',bandL{band+1},camcol));
                if band==4 || col>0
                    xlabel('seeing');
                end
                if camcol==1 || col>0
                    ylabel('\alpha (\lambda dependence)');
                end
            end
        end
        if col<0
            fname = sprintf('../SDSSdata/%s/Stripe82_%s_col%d.txt',usebands,bandL{band+1},col);
            i1 = 1;
            i2 = 6;
            if col==-25
               i1=2;
               i2=5;
            end
            if writeMaster
                % write out master file
                fid = fopen(fname, 'w');
                fprintf(fid,'run \t alpha \t alphaErr \t seeing\n');
                for irun = 1:nrun
                    fprintf(fid, '%d \t %.2f \t %.4f \t %.3f\n', runlist(irun), ...
                        sum(alpha(irun,i1:i2)./alphaErr(irun,i1:i2).^2)/sum(1./alphaErr(irun,i1:i2).^2), ...
                        1./sqrt(sum(1./alphaErr(irun,i1:i2).^2)), ...
                        mean(seeing(irun,i1:i2)));
                end
                fclose(fid);
            end
        end
    end
end
end
