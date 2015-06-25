%% this script makes preliminary ADCP plots to determine time endpoints and bin endpoints.
% load data from WinADCP before starting.
% answers should be put into a file called endpoints.txt

Pitch = AnP100thDeg *.01;
Roll = AnR100thDeg *.01;

East_u = SerEmmpersec * 0.1; % convert mm per sec to cm per sec
North_v = SerNmmpersec * 0.1;
East_u(East_u == -3276.8) = NaN;
North_v(North_v == -3276.8) = NaN;

VelError = SerErmmpersec * 0.1;
VelError(VelError == -3276.8) = NaN;

CorrAvg = SerCAcnt;
EchoAvg = SerEAAcnt;
PG1 = SerPG1;
PG4 = SerPG4;

% flags
[i, j] = find(CorrAvg < 64);
badcorr = [i j];
clear i j
[i, j] = find(EchoAvg < 40);
badecho = [i j];
clear i j
[i, j] = find((PG1 + PG4) < 80);
badPG = [i j];
clear i j
[i, j] = find(VelError > 3.7); % .74 cm per sec st.dev * 5
badM = [i j];
clear i j

badmaxpitch = find(Pitch > 20 | Pitch < -20);
badmaxroll = find(Roll > 20 | Roll < -20);
badchgpitch = find(diff(Pitch>5));
badchgroll = find(diff(Roll>5));

% time
timestamps = datenum(SerYear+2000, SerMon, SerDay, SerHour, SerMin+2.5, 0);

excess = datenum(SerYear(1)+2000, 0, 1, 0, 0, 0); % Jan 1 is Day 0
for i = 1:length(SerYear)
    jday(i) = timestamps(i) - excess;
end
clear i

figure(100)
plot(Pitch, '.-b')
hold on
plot(Roll, '.-r')
% use "xlim([min max])" to zoom in on x axis.

%stickplots by bin
chg = floor(timestamps);
daychg=[1; find(diff(chg)>0); length(chg)];
rangeD = ceil(jday(end)) - floor(jday(1));

d = length(SerBins);
for A = 1:10
    % position: left bottom width height
    figure('Position', [400 80 1100 850]);
    orient landscape
    for s = 1:2
        if d > 0
            subplot(2,1,s)
            plot([jday(1)-rangeD*.05 jday(end)+rangeD*.05], [0 0], 'k')
            axis([jday(1)-rangeD*.05 jday(end)+rangeD*.05 -8 8])
            daspect([1 1 1])
            hold on
            
            stickcolor = 'r';
            for e = 1:length(daychg) - 1
                stickhandle(e) = stickplot(jday(daychg(e):daychg(e+1)), North_v(daychg(e):daychg(e+1),d)*.1, -East_u(daychg(e):daychg(e+1),d)*.1);
                set(stickhandle(e), 'color', stickcolor);
                if stickcolor == 'r'
                    stickcolor = 'g';
                elseif stickcolor == 'g'
                    stickcolor = 'b';
                elseif stickcolor == 'b'
                    stickcolor = 'r';
                end
            end
            
            for t = {badcorr, badecho, badPG, badM}
                bads = t{1};
                binbads = find(bads(:, 2) == d);
                badstick = stickplot(jday(bads(binbads, 1)), North_v(bads(binbads, 1), d)*.1, -East_u(bads(binbads, 1), d)*.1);
                set(badstick, 'Color', [.8 .8 .8])
                clear bads binbads badstick
            end
            
            for j = {badmaxpitch, badmaxroll, badchgpitch, badchgroll}
                bads = j{1};
                badstick = stickplot(jday(bads), North_v(bads, d)*.1, -East_u(bads, d)*.1);
                set(badstick, 'color', [.8 .8 .8])
                clear bads badstick
            end
            
            set(gca, 'ytick', [-8 -4 0 4 8])
            set(gca, 'yticklabel', [-80 -40 0 40 80])
            ylabel('(cm s^-^1)')

            daxis = axis;
            text(daxis(1)+1, 7, ['bin # ', num2str(SerBins(d))]);

            if s == 1
                title('up = W');
            end

            clear e stickcolor stickhandle t j
            d = d-1;
       end % end if bin exists
    end % for s subplots
xlabel(['Days since 00:00 01/01 ', num2str(SerYear(1)+2000), ' (', datestr(timestamps(1), 'mm/dd/yyyy'), ')'])
text(daxis(2) - (daxis(2)-daxis(1))*.2, -12, datestr(now));

clear daxis

end
clear d A s
