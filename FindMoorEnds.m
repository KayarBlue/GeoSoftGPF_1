%% helper script for manually finding scan start and stop for moored CTDs.
% importdata the -matlab.txt, then copy-paste below:

%YSI
depth = YSI(:,6);
sal = YSI(:,7);
temp = YSI(:,8);
scan = 1:length(depth);

%CTDs
temp = CTD(:,1);
sal = CTD(:,4);
depth = CTD(:,9);
scan = 1:length(depth);

% figure screensize = [left bottom width height]
figure('Position', [800 110 672 798]);
ax(1) = subplot(3,1,1);
plot(scan,depth,'.-k')
ax(2) = subplot(3,1,2);
plot(scan,sal,'.-b')
ax(3) = subplot(3,1,3);
plot(scan,temp,'.-r')
linkaxes(ax,'x');

% proceed by manipulating temp

%% thermistors
% Before starting: uiimport all data and name them temp1 thru temp5
figure('Position', [130 420 1115 500]);
plot(getfield(temp1,'data'),'k.-')
hold on
plot(getfield(temp2,'data'), 'g.-')
plot(getfield(temp3,'data'),'r.-')
plot(getfield(temp4,'data'),'y.-')
plot(getfield(temp5,'data'),'b.-')

maxscan = max([length(getfield(temp1,'data')) length(getfield(temp2,'data')) length(getfield(temp3,'data')) length(getfield(temp4,'data')) length(getfield(temp5,'data'))]);
axis([0 maxscan 10 26])

% proceed by changing axis
