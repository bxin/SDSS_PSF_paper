function [] = alpha_seeing(iband, instfwhm)

% iband: on the x-axis, we plot iband seeing
%
% to reproduce the plot in Exp Astron (2016)42:85, 
%       use the correct Leff and L0 below, and iband=2
%       and search for "Astron"
% to produce results for use by Python code for plotting
%      switch back Leff, and L0, and use iband = 3 (r-band)

% example (for r-band):
%      alpha_seeing(2)

% airmass is always 1. 
% Our input variables are: L0, r0(500), lambda.
% The output is the vK FWHM.
zen = 0;

if nargin<2
    instfwhm = 0;
end

% r0inm500 = 0.02:0.01:0.52;
r0inm500 = [0.02:0.002:0.12 0.14:0.02:0.5];
Leff = [355.1 468.6 616.6 748.0 893.2];
% Leff = [355.1 550 616.6 748.0 2200]; %%!!! 893.2];
% L0 = [10 20 25 30 40 50 80 100 1e9];
% L0 = [25 30 50 100 1e9];
L0 = [2 5 10 20 30 50 100 1e9];

nL0 = length(L0);
nr0 = length(r0inm500);
data = zeros(nL0, 3, nr0);

figure(1);clf;

for iL0=1:nL0 % for every L0
    yL0 = L0(iL0);
    for ir0=1:nr0 % for every r0(500)
        % we need to get a slope alpha (fwhm vs lambda)
        yr0inm500 = r0inm500(ir0);
        yr0=yr0inm500*cos(zen)^0.6.*(Leff/500).^1.2;

        fwhm0=0.98*Leff*1e-9./(yr0/3600/180*pi);
        f2 = 1-2.183*(yr0/yL0).^0.356;
        f2(f2<0)=0;
        fwhm = fwhm0.*sqrt(f2);
        fwhm = sqrt(fwhm.^2+instfwhm^2);
        if fwhm(iband)==0
            alpha = nan;
        else
            ratio=fwhm/fwhm(iband);
            
            alpha0=-0.2;
            x=reshape(Leff,1,[]);
            z=reshape(ratio,1,[]);
            lb=-0.8;
            ub=0.3;
            options = optimoptions('lsqcurvefit','Display','off');
            % fprintf('iL0=%d, ir0=%d\n', iL0, ir0);
            alpha=lsqcurvefit(@(alpha, lambda)lambdapower(alpha, lambda, Leff(iband)),...
                alpha0,x,z,lb,ub,options);
            %         lfit = 300:1000;
            %         rfit = lambdapower(alpha,lfit);
            %         plot(Leff, ratio,'ko', lfit, rfit, '-r');
            %         data(iL0, 1, ir0) = fwhm0(iband);
        end
        data(iL0, 1, ir0) = fwhm(iband);
        data(iL0, 2, ir0) = alpha;

        % below for comparing with Exp Astron (2016) 42:85-98, Fig 4.
        % data(iL0, 2, ir0) = (Leff(iband)/Leff(end)).^alpha;
        %fprintf('alpha = %f\n', alpha);
        % data(iL0, 3, ir0) = fwhm(iband)/fwhm(end);
    end
    x = squeeze(data(iL0,1,:));
    y = squeeze(data(iL0,2,:));
    y1 = squeeze(data(iL0,3,:));
    plot(x,y,'-r+');
    hold on;
    plot(x,y1,'b');
    if iL0==1
       legend({'fit','ratio'}); 
    end
end
hold off;
grid;
% below for comparing with Exp Astron (2016) 42:85-98, Fig 4.
% xlim([0.2, 1.6]);
% ylim([1.1, 2]);
% xlabel('Observed Seeing V [arcsec]');
% ylabel('Ratio of observed seeing V/K');

a= squeeze(data(:,1,:));
b= squeeze(data(:,2,:));
a= [L0' a];
b= [L0' b];
save('../data/alpha_seeing_x.txt','a','-ascii');
save('../data/alpha_seeing_y.txt','b','-ascii');

end

function F = lambdapower(alpha, lambda, lambda0)

F = (lambda/lambda0).^alpha;
% fprintf('alpha=%f, lambda0=%f \n',alpha, lambda0);
end