M = importdata('../SDSSdata/masterTXT/run3427.txt', ' ', 1);

idxu=(M.data(:,2)==1 &M.data(:,3)==0);
idxg=(M.data(:,2)==1 &M.data(:,3)==1);
idxr=(M.data(:,2)==1 &M.data(:,3)==2);
idxi=(M.data(:,2)==1 &M.data(:,3)==3);
idxz=(M.data(:,2)==1 &M.data(:,3)==4);

subplot(2,1,1);
yu = M.data(idxu,4)/1.222./M.data(idxu,6).^0.6;
yg = M.data(idxg,4)/1.222./M.data(idxg,6).^0.6;
yr = M.data(idxr,4)/1.222./M.data(idxr,6).^0.6;
yi = M.data(idxi,4)/1.222./M.data(idxi,6).^0.6;
yz = M.data(idxz,4)/1.222./M.data(idxz,6).^0.6;
n=length(yu);
x=1:n;
plot(x,yu,'-b',x,yg,'-g',x,yr,'-r',x,yi,'-m',x,yz,'-k');
legend('u','g','r','i','z');

% subplot(2,2,3);
% plot(x,yu./yr);

subplot(2,2,3);
xm = 1:n-8;
rat = yu(5:end-4)./yr(9:end);
plot(xm, rat);