function [] = correlate(writeMaster, col, filter)

usebands = 'ugriz';
% usebands = 'ugri';
figure(1);clf;
figure(2);clf;
figure(3);clf;
figure(4);clf;
figure(5);clf;
figure(6);clf;
bandL = {'u','g','r','i','z'};
% camcol = 3; % use camcol #3, b/c it is least affected by ccd edge effects(?)
% band = 2; % use r band to examime the correlations
for camcol = 1:6
    if (col<0 || col == camcol)
        for band=0:4
            if (filter<0 || band == filter)
                fname = sprintf('../SDSSdata/%s/Stripe82_%s_col%d.txt',usebands,bandL{band+1},camcol);
                if writeMaster
                    % alpha (fwhm vs lambda)
                    figure(7);subplot(2,2,1);
                    a=load(sprintf('../output/fwhm_lambda/power_%s.txt',usebands));
                    runlist = a(:,1);
                    alpha = a(:,camcol+1);
                    hist(alpha);
                    % title(sprintf('alpha in camcol#3: mean = %.2f, std=%.2f\n', mean(alpha), std(alpha)));
                    
                    % beta (spatial correlation slope)
                    figure(7);subplot(2,2,2);
                    a = load('../output/correlate_spatial/run-9_fwhm_fitp.txt');
                    beta = a(a(:,2)==band, 3);
                    hist(beta);
                    x=-0.12:0.002:0;
                    hist(beta,x);
                    grid;
                    
                    % tau (time constant that seeing decorrelates)
                    figure(7);subplot(2,2,3);
                    a = load('../output/correlate_temporal/run-9_fwhm_fitp.txt');
                    idx = (a(:,2)==band) & (a(:,3)==camcol);
                    tau = a(idx, 6);
                    x=0:20:1600;
                    hist(tau,x);
                    
                    % average seeing for each run
                    figure(7);subplot(2,2,4);
                    nrun = length(runlist);
                    seeing = zeros(nrun,1);
                    for irun = 1:nrun
                        run = runlist(irun);
                        a = importdata(sprintf('../SDSSdata/masterTXT/run%d.txt', run),' ', 1);
                        a = a.data;
                        idx  = (a(:,2)==camcol) & (a(:,3)==band);
                        seeing(irun) = mean(a(idx,4));
                    end
                    hist(seeing);
                    
                    % write out master file
                    fid = fopen(fname, 'w');
                    fprintf(fid,'run \t alpha \t \beta \t tau \t seeing\n');
                    for irun = 1:nrun
                        fprintf(fid, '%d \t %.2f \t %.4f \t %8.1f \t %.3f\n', runlist(irun), alpha(irun), beta(irun), tau(irun),seeing(irun));
                    end
                    fclose(fid);
                else
                    a = importdata(fname,' ',1);
                    a = a.data;
                    runlist = a(:,1);
                    alpha = a(:,2);
                    beta = a(:,3);
                    tau = a(:,4);
                    seeing = a(:,5);
                end
                
                figure(1);
                if col<0  
                    subplot(5,6,band*6+camcol)
                end
                scatter(alpha,beta,400,'.');grid;title(sprintf('%s, col%d',bandL{band+1},camcol));
                xlabel('alpha (wavelength dependence)'); ylabel('spatial correlation slope');
                figure(2);
                if col<0
                    subplot(5,6,band*6+camcol)
                end
                scatter(alpha,tau,400,'.');grid;title(sprintf('%s, col%d',bandL{band+1},camcol));
                xlabel('alpha (wavelength dependence)'); ylabel('decorrelation timescale');
                figure(3);
                if col<0
                    subplot(5,6,band*6+camcol)
                end
                scatter(beta,tau,400,'.');grid;title(sprintf('%s, col%d',bandL{band+1},camcol));
                xlabel('spatial correlation slope'); ylabel('decorrelation timescale');
                
                figure(4);
                if col<0
                    subplot(5,6,band*6+camcol)
                end
                scatter(seeing, beta,400,'.');grid;title(sprintf('%s, col%d',bandL{band+1},camcol));
                xlabel('seeing'); ylabel('spatial correlation slope');
                figure(5);
                if col<0
                    subplot(5,6,band*6+camcol)
                end
                scatter(alpha,tau,400,'.');grid;title(sprintf('%s, col%d',bandL{band+1},camcol));
                xlabel('seeing'); ylabel('decorrelation timescale');
                
                figure(6);
                if col<0  
                    subplot(5,6,band*6+camcol)
                end
                scatter(seeing, alpha,400,'.');grid;title(sprintf('%s, col%d',bandL{band+1},camcol));
                if band==4 || col>0
                    xlabel('seeing');
                end
                if camcol==1 || col>0
                    ylabel('\alpha (\lambda dependence)');
                end
            end
        end
    end
end
end
