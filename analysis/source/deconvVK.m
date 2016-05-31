function [] = deconvVK(test)

% test=1, in test mode, convolve deconvolved-PSF back and see if we get g2 (doubleG)
% test=0, skip test to make it faster

x=1:501;
x=x*0.1;
[xm, ym] = meshgrid(-50:0.1:50);

% we take run94, field#1 as example, psf.scaleR = 1.4158208768975589
vk = load('~/wavefront/activeoptics/matlab/pssn/data/vonK1.0.txt');
scaleR = 1.4158208768975589;
scaleV = 1.0350848010744576;

%for crosschecking with python output
if test
    vkN=vk/sum(vk(:));
    fwhmeffCheck = sqrt(1/sum(vkN(:).^2))*0.664*0.1; %0.1arcsec/pixel
    fwhmeffvk = fwhmeffCheck*scaleR;
    fprintf('fwhmeffvk = %4.2f\n',fwhmeffvk);
end
xm1=xm*scaleR;
ym1=ym*scaleR;

vk=interp2(xm1,ym1,vk*scaleV,xm,ym);

% this file was generated specifically for this test, using fwhm_check_eff.py
% -m pdb, stop at line#83
g2 = load('../output/run94_psf_2G_2D.txt');
%for crosschecking with python output
if test
    g2N=g2/sum(g2(:));
    fwhmeffg2 = sqrt(1/sum(g2N(:).^2))*0.664*0.1; %0.1arcsec/pixel
    fprintf('fwhmeff of g2: %4.2f\n',fwhmeffg2);
end

figure(1);clf;

if test
    %proper normalization in order to reproduce the plot
    vk=vk/max(vk(:));
    g2=g2*scaleV;

    subplot(2,3,1);
    imagesc(vk);axis square;title('von Karman');
    subplot(2,3,2);
    imagesc(g2);axis square;title('Double Gaussian');
    subplot(2,3,3);
    % reproduce the python plot
    semilogy(x,vk(501,501:end),'-b',x,g2(501,501:end),'-r');
    legend('vonk','2G');
    grid on;
    ylim([1e-6 10]);
    xlim([0 30]);
    xlabel('arcsec');
end

% now we can deconv
g2=g2/sum(g2(:));
vk=vk/sum(vk(:));
m = 1001;
vkfft = fftshift(fft2(fftshift(vk),m,m)); 
g2fft = fftshift(fft2(fftshift(g2),m,m)); 
psfi = abs(fftshift(ifft2(fftshift(g2fft./vkfft),m,m))); % OTF to PSF
psfi = psfi/max(psfi(:));

% % %convolve back and check
if test
    subplot(2,3,3);
    % %using conv2 directly
    % g2rec = conv2(vk,extractArray(psfi,200),'same');
    g2rec = conv2(vk,psfi,'same');
    % % use fft to convolve back
    % prodfft=fftshift(fft2(fftshift(psfi),m,m)).*vkfft;
    % g2rec = abs(fftshift(ifft2(fftshift(prodfft),m,m))); % OTF to PSF
    
    g2rec=g2rec/max(g2rec(:));
    hold on;
    semilogy(x,g2rec(501,501:end),'-g','linewidth',5);
end

sigma=0.4;
A1=0; %1e-3; %3e-5;
sigma1=2;
G1=A1*exp(-x.^2/2/sigma1^2);
A2=2e-5;
sigma2=14;
G2=A2*exp(-x.^2/2/sigma2^2);
g=exp(-x.^2/2/sigma^2)+G1+G2;

subplot(2,3,4);
semilogy(x,psfi(501,501:end),'-ro',x,g,'-b');
ylim([1e-6 10]);
xlim([0 10]);
hold on;
semilogy(x,G1,'-go');
grid on;

xlabel('arcsec');
text(4,0.1,sprintf('Gaussian sigma=%3.1f arcsec',sigma));

subplot(2,3,5);

plot(x,psfi(501,501:end),'-ro',x,g,'-b');
grid on;
ylim([0 1]);
xlim([0 2]);
xlabel('arcsec');
text(0.4,0.8,'Gaussian sigma=0.2 arcsec','units','normalized');

% [xm, ym] = meshgrid(-5:0.1:5);
g2d=exp(-(xm.^2+ym.^2)/2/sigma^2)+A1*exp(-(xm.^2+ym.^2)/2/sigma1^2)+A2*exp(-(xm.^2+ym.^2)/2/sigma2^2);
% new = conv2(vk,g2d,'same');
prodfft=fftshift(fft2(fftshift(g2d),m,m)).*vkfft;
new = abs(fftshift(ifft2(fftshift(prodfft),m,m))); % OTF to PSF
 
new=new/max(new(:));
newslice = new(501,501:end);
subplot(2,3,6);
% reproduce the python plot
vk=vk/max(vk(:));
g2=g2/max(g2(:));
semilogy(x,vk(501,501:end),'-b',x,g2(501,501:end),'-r');
legend('vonk','2G');
grid on;
ylim([1e-6 10]);
xlim([0 30]);
xlabel('arcsec');
hold on;
semilogy(x,newslice(1:501),'-g');

end

