function [] = conv1d_vs_2d(test)

% test = 'gau'
%      = 'linear'

[x,y]= meshgrid(-50:0.1:50);
mid=501;
r2 = (x.^2+y.^2);
s1=2;
s2=5;
s12=sqrt(s1^2+s2^2);

g1=exp(-r2/2/s1^2);
if strcmp(test, 'gau')
    g2=exp(-r2/2/s2^2);
elseif strcmp(test, 'linear')
    g2 = -sqrt(r2)+s2*4;
    g2(g2<0)=0;
end
g12=conv2(g1,g2,'same');
g1=g1/sum(g1(:));
g2=g2/sum(g2(:));
g12=g12/sum(g12(:));

s1_cp=0.664*0.1*sqrt(1./sum(sum(g1.^2)))/2.3548;
s2_cp=0.664*0.1*sqrt(1./sum(sum(g2.^2)))/2.3548;
s12_cp=0.664*0.1*sqrt(1./sum(sum(g12.^2)))/2.3548;

fprintf('input s1=%.4f, measured s1=%.4f\n', s1, s1_cp);
fprintf('input s2=%.4f, measured s2=%.4f\n', s2, s2_cp);
fprintf('RSS s12 = %.4f, measured s12=%.4f\n', s12, s12_cp);

g11d=g1(mid,:);
g21d=g2(mid,:);
g121d=g12(mid,:);
g121d=g121d/sum(g121d);
g121d_cp = conv(g11d, g21d, 'same');
g121d_cp=g121d_cp/sum(g121d_cp);

xx=1:1001;
clf;plot(xx,g121d,'-r',xx,g121d_cp,'b.');

end
