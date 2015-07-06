%% helper script for manually finding maxdepth in CTD casts.

% importdata the -matlab.txt, check that the columns match the data, then copy-paste below:
scan = data(:,1);
depth = data(:,3);
sal = data(:, 4);
temp = data(:,5);
figure('Position', [1042 110 672 798]);
ax(1) = subplot(3,1,1);
plot(scan,depth,'.-k')
ax(2) = subplot(3,1,2);
plot(scan,sal,'.-b')
ax(3) = subplot(3,1,3);
plot(scan,temp,'.-r')
linkaxes(ax,'x');

% proceed by manipulating temp
