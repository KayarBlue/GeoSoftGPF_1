%% This script to process Mooring ADCP data. 
% Requires Matlab 2010b with Mapping Toolbox. Uses the function stickplot.m.
% input file is a *.MAT that includes everything under "Series Data Types"
% and Pitch, Roll, Heading, and Battery under "Anc Data Types". This is
% exported from the *.000 datafile with software WinADCP. 
% Also need to have depth data from CTDBOT loaded. This exists in the
% ctdburst cell array from the Mooring processing (MOORprocess_all.m).
% Also need to have endpoints.txt loaded.

%endpoints.txt: start and stop scan of timestamps, last bin.
timestart = data(1,1);
timestop = data(1,3);
lastbin = data(1,4);
lastbindex = find(SerBins == lastbin);

%dep = 2.07:RDIBinSize:SerBins(end)*RDIBinSize;

for i = 1:length(SerBins)
    dep(i) = 2.07 + .5*(i-1);
end

% ADCP variables
Pitch = AnP100thDeg *.01;
Roll = AnR100thDeg *.01;

East_u = SerEmmpersec * 0.1; % convert mm per sec to cm per sec
North_v = SerNmmpersec * 0.1;
Vert_w = SerVmmpersec * 0.1;
VelError = SerErmmpersec * 0.1;
Mag = SerMagmmpersec * 0.1;
DIR = SerDir10thDeg * 0.1;

East_u(East_u == -3276.8) = NaN;
North_v(North_v == -3276.8) = NaN;
Vert_w(Vert_w == -3276.8) = NaN;
VelError(VelError == -3276.8) = NaN;
Mag(Mag == -3276.8) = NaN;
DIR(DIR == -3276.8) = NaN;

interpnum = cell(1,6);
fixed{1} = East_u;
fixed{2} = North_v;
fixed{3} = Vert_w;
fixed{4} = VelError;
fixed{5} = Mag;
fixed{6} = DIR;
for s = 1:6
    blort = fixed{s};
    k = 1; %within blort-cell counter
    for i = 1:length(SerBins)
        badeast = find(isnan(blort(:,i)));
        for j = 1:length(badeast)
            m = badeast(j);
            while isnan(blort(m, i))
                if m > 1
                    m = m - 1;
                else
                    break
                end
            end
            n = badeast(j);
            while isnan(blort(n, i))
                if n < length(SerYear)
                    n = n + 1;
                else
                    break
                end
            end
            interpnum{s}(k, 1) = badeast(j);
            interpnum{s}(k, 2) = i;
            if n - m < 6 && ~isnan(blort(m, i)) && ~isnan(blort(n, i))
                interpnum{s}(k, 3) = interp1([m n], [blort(m, i) blort(n, i)], badeast(j));
            elseif n - m < 6 && isnan(blort(m, i)) && ~isnan(blort(n, i))
                interpnum{s}(k, 3) = blort(n, i);
            elseif n - m < 6 && ~isnan(blort(m, i)) && isnan(blort(n, i))
                interpnum{s}(k, 3) = blort(m, i);
            else
                interpnum{s}(k, 3) = NaN;
            end
            k = k + 1;
        end
        clear badeast j m n
    end
    clear i blort k
    for f = 1:length(interpnum{s})
        fixed{s}(interpnum{s}(f, 1), interpnum{s}(f, 2)) = interpnum{s}(f, 3);
    end
end
fixedEast = fixed{1};
fixedNorth = fixed{2};
fixedVert = fixed{3};
fixedErr = fixed{4};
fixedMag = fixed{5};
fixedDIR = fixed{6};
clear s f fixed

CorrAvg = SerCAcnt;
EchoAvg = SerEAAcnt;

PG1 = SerPG1;
PG4 = SerPG4;

% flags
[i, j] = find(CorrAvg(timestart:timestop, 1:lastbindex) < 64); %64+ is good
badcorr = [i j];
clear i j

[i, j] = find(EchoAvg(timestart:timestop, 1:lastbindex) < 40); %40+ is good
badecho = [i j];
clear i j

[i, j] = find((PG1(timestart:timestop, 1:lastbindex) + PG4(timestart:timestop, 1:lastbindex)) < 80); %80+ is good
badPG = [i j];
clear i j

[i, j] = find(VelError(timestart:timestop, 1:lastbindex) > 3.7); % .74 cm per sec st.dev * 5
badM = [i j];
clear i j

badmaxpitch = find(Pitch(timestart:timestop) > 20 | Pitch(timestart:timestop) < -20);
badmaxroll = find(Roll(timestart:timestop) > 20 | Roll(timestart:timestop) < -20);
badchgpitch = find(diff(Pitch(timestart:timestop)>5));
badchgroll = find(diff(Roll(timestart:timestop)>5));

% time
timestamps = datenum(SerYear+2000, SerMon, SerDay, SerHour, SerMin+2.5, SerSec);

excess = datenum(SerYear(1)+2000, 0, 1, 0, 0, 0); % Jan 1 is Day 0. No rollover on new year. 
for i = 1:length(SerYear)
    jday(i) = timestamps(i) - excess;
end
clear i excess

if ~isempty(ctdburst{1})
inbotstamps = ctdburst{1}(:,2);
inbotdepth = ctdburst{1}(:,3);
inbotjday = ctdburst{1}(:,11) - 1; % Jan 1 = 0

if inbotjday(1) < jday(timestart)
    firstbot = find(datenum(datestr(inbotjday)) > datenum(datestr(jday(timestart))) - datenum(0,0,0,0,5,0) & datenum(datestr(inbotjday)) < datenum(datestr(jday(timestart))) + datenum(0,0,0,0,5,0));
    firstbottrack = timestart; %refers to ADCP timetrack
else
    firstbot = 1; %refers to CTD timetrack
    firstbottrack = find(datenum(datestr(jday)) > datenum(datestr(inbotjday(1))) - datenum(0,0,0,0,5,0) & datenum(datestr(jday)) < datenum(datestr(inbotjday(1))) + datenum(0,0,0,0,5,0));
end
if inbotjday(end) > jday(timestop)
    lastbottrack = timestop;
    lastbot = find(datenum(datestr(inbotjday)) > datenum(datestr(jday(timestop))) - datenum(0,0,0,0,5,0) & datenum(datestr(inbotjday)) < datenum(datestr(jday(timestop))) + datenum(0,0,0,0,5,0));
else
    lastbottrack = find(datenum(datestr(jday)) > datenum(datestr(inbotjday(end))) - datenum(0,0,0,0,5,0) & datenum(datestr(jday)) < datenum(datestr(inbotjday(end))) + datenum(0,0,0,0,5,0));
    lastbot = length(inbotjday);
end
%for the color plots
botdepth = inbotdepth(firstbot:lastbot);
botjday = inbotjday(firstbot:lastbot);

% for the output dat file
bottrack(1:length(jday)) = NaN; 
bottrack(firstbottrack:lastbottrack) = botdepth;

else 
botdepth(1:length(SerYear)) = NaN;
botjday = jday;
bottrack(1:length(jday)) = NaN;
end

%% the big dat file
bigdat = fopen([datestr(timestamps(timestart), 'yyyymmdd'), '_', datestr(timestamps(timestop), 'yyyymmdd'), '.dat'],'w');
fprintf(bigdat, 'n,    bin,    y,    m,    d,    h,    min,    t_day,        tdep,    dep,    u(cm/s),   v,    spd,    dir,    ECMP    Interpolated\n');

bigcsv = fopen([datestr(timestamps(timestart), 'yyyymmdd'), '_', datestr(timestamps(timestop), 'yyyymmdd'), '.csv'], 'w');
fprintf(bigcsv, 'Scan, Bin, Timestamp, East_u, North_v, Vert_w, Speed, Direction\n');

nlength = length(num2str(timestop)) + 2;
mlength = length(num2str(lastbindex)) + 2;
timelength = 6;
jdaylength = length(num2str(floor(max(jday)))) + 4;
trklength = length(num2str(floor(max(bottrack)))) + 4;
deplength = length(num2str(floor(max(dep)))) + 4;
eastlength = length(num2str(floor(max(max(fixedEast))))) + 4;
northlength = length(num2str(floor(max(max(fixedNorth))))) + 4;
maglength = length(num2str(floor(max(max(fixedMag))))) + 4;
dirlength = length(num2str(floor(max(max(fixedDIR))))) + 4;

for n = timestart:timestop
        [yrj moj dayj hrj mnj ~] = datevec(timestamps(n));
for m = 1:lastbindex
    qpad = nlength - length(num2str(n));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%i', n);
    fprintf(bigcsv, '%i,', n);

    qpad = mlength - length(num2str(m));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%i', SerBins(m));
    fprintf(bigcsv, '%i,', SerBins(m));

    if ~isnan(timestamps(n))
        fprintf(bigcsv, '%s,', datestr(timestamps(n)));
    else
        fprintf(bigcsv, 'NaN,');
    end
    
    qpad = timelength - length(num2str(yrj));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%i', yrj);
    qpad = timelength - length(num2str(moj));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%i', moj);
    qpad = timelength - length(num2str(dayj));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%i', dayj);
    qpad = timelength - length(num2str(hrj));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%i', hrj);
    qpad = timelength - length(num2str(mnj));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%i.5', mnj);
   
    qpad = jdaylength - length(num2str(floor(jday(n))));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%10.6f', jday(n));

    qpad = trklength - length(num2str(floor(bottrack(n))));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%6.3f', bottrack(n));

    qpad = deplength - length(num2str(floor(dep(m))));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%4.2f', dep(m));

    qpad = eastlength - length(num2str(floor(fixedEast(n,m))));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%4.2f', fixedEast(n, m));
    fprintf(bigcsv, '%4.2f,', fixedEast(n, m));

    qpad = northlength - length(num2str(floor(fixedNorth(n,m))));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%4.2f', fixedNorth(n, m));
    fprintf(bigcsv, '%4.2f,', fixedNorth(n, m));

    fprintf(bigcsv, '%4.2f,', fixedVert(n, m));

    qpad = maglength - length(num2str(floor(fixedMag(n,m))));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%5.2f', fixedMag(n, m));
    fprintf(bigcsv, '%5.2f,', fixedMag(n, m));

    qpad = dirlength - length(num2str(floor(fixedDIR(n,m))));
    for pad = 1:qpad
        fprintf(bigdat, ' ');
    end
    clear *pad
    fprintf(bigdat, '%5.2f    ', fixedDIR(n, m));
    fprintf(bigcsv, '%5.2f\n', fixedDIR(n, m));

    if EchoAvg(n, m) < 40
        fprintf(bigdat, '*');
    else
        fprintf(bigdat, '.');
    end
    if CorrAvg(n, m) < 64
        fprintf(bigdat, '*');
    else
        fprintf(bigdat, '.');
    end
    if VelError(n, m) > 37
        fprintf(bigdat, '*');
    else
        fprintf(bigdat, '.');
    end
    if (PG1(n, m) + PG4(n, m)) < 80
        fprintf(bigdat, '*    ');
    else
        fprintf(bigdat, '.    ');
    end
    if isnan(East_u(n, m))
        fprintf(bigdat, '*');
    else
        fprintf(bigdat, '.');
    end
    if isnan(North_v(n, m))
        fprintf(bigdat, '*');
    else
        fprintf(bigdat, '.');
    end
    if isnan(Mag(n, m))
        fprintf(bigdat, '*');
    else
        fprintf(bigdat, '.');
    end
    if isnan(DIR(n, m))
        fprintf(bigdat, '*\n');
    else
        fprintf(bigdat, '.\n');
    end
end
    clear yrj moj dayj hrj mnj
end
fclose(bigdat);
fclose(bigcsv);
clear n m bigdat bigcsv *length *pad

%% figures by bin

chg = floor(timestamps);
chg2=[1; find(diff(chg)>0); length(chg)];
daychg = chg2(chg2 >= timestart & chg2 <= timestop);

rangeD = ceil(jday(end)) - floor(jday(1));

%for count of NaNs in scatterplot titles, zavg bar plot
totbins(1:length(jday), 1:length(SerBins)) = 0;
binbar(1:length(jday)) = 0;
timebar(1:length(SerBins)) = 0;

for c = 1:length(SerBins)
    timebar(c) = length(find(totbins(timestart:timestop,c) ==0));
end
clear a b c

a = lastbindex;
b = lastbindex;
c = lastbindex;
d = lastbindex;
numfigs = ceil((SerBins(lastbindex)-SerBins(1)+1)/2);
for A = 1:numfigs %start bins loop

% 0_scatter.pdf
figure(A)
orient landscape
for p = 1:2
if a > 0
subplot(1, 2, p)
plot(fixedEast(timestart:timestop, a), fixedNorth(timestart:timestop,a),'.r')
hold on
plot([-80 80], [0 0], 'k')
plot([0 0], [-80 80], 'k')
axis equal

for t = {badcorr, badecho, badPG, badM}
    bads = t{1};
    binbads = find(bads(:, 2) == a);
    baddot = plot(fixedEast(bads(binbads, 1), a), fixedNorth(bads(binbads, 1), a),'.');
    set(baddot, 'Color', [.8 .8 .8])
    clear bads binbads baddot
end

for j = {badmaxpitch, badmaxroll, badchgpitch, badchgroll}
    bads = j{1};
    baddot = plot(fixedEast(bads, a), fixedNorth(bads, a), '.');
    set(baddot, 'color', [.8 .8 .8])
    clear bads baddot
end

% title is supposed to be "bin: non-missing (all)"
title([num2str(a), ': ', num2str(timebar(a)), ' (', num2str(timestop-timestart+1), ')']);
xlabel('u (cm s^-^1)')
ylabel('v (cm s^-^1)')
axis([-80 80 -80 80])
grid

clear t j
a = a-1;
end % if a exists
end % for p of fig A
text(20, -120, datestr(now));
saveas(A, ['scatter_', num2str(A), '.pdf'], 'pdf');

% 4_uv.pdf
B = numfigs + A;
figure(B)
orient landscape
for q = 1:2
if b > 0
subplot(2,1,q)
plot(jday(timestart:timestop), fixedEast(timestart:timestop,b),'b')
hold on
plot(jday(timestart:timestop), fixedNorth(timestart:timestop,b), 'r')
plot(jday, 0, 'k')
axis([jday(timestart) jday(timestop) -80 80])
ylabel('u: blue (cm s^-^1)')
text(jday(daychg(2)), 70, ['bin # ', num2str(SerBins(b))]);

b = b-1;
end % end if b exists

end % end for q of fig B
xlabel(['Days since 00:00 01/01 ', num2str(SerYear(timestart)+2000), ' (', datestr(timestamps(timestart), 'mm/dd/yyyy'), ')'])
text(jday(end) - (0.2*(jday(end) - jday(1))), -120, datestr(now));

saveas(B, ['uv_', num2str(A), '.pdf'], 'pdf');

% 3_spd.pdf: mag, dir
C = numfigs*2 + A;
figure(C)
orient landscape
for r = 1:2
    if c > 0
        subplot(2,1,r)
        [Amagdir, Hmag, Hdir] = plotyy(jday(timestart:timestop), fixedMag(timestart:timestop,c), jday(timestart:timestop), fixedDIR(timestart:timestop,c));
        set(get(Amagdir(1), 'Ylabel'), 'String', 'Speed (cm s^-^1)')
        set(Amagdir(1), 'Ylim', [0 80])
        set(Amagdir(1), 'Ycolor', [0 0 1])
        set(Amagdir(1), 'ytick', [0 20 40 60 80])
        set(Hmag, 'color', 'b')
        set(get(Amagdir(2), 'Ylabel'), 'String', 'Direction')
        set(Amagdir(2), 'Ylim', [0 360])
        set(Amagdir(2), 'Ycolor', [1 0 0])
        set(Amagdir(2), 'ytick', [0 90 180 270 360])
        set(Hdir, 'color', 'r')
        text(jday(daychg(2)), 75, ['bin # ', num2str(SerBins(c))]);
        
        clear Amagdir Hmag Hdir
        c = c-1;
    end
end
xlabel(['Days since 00:00 01/01 ', num2str(SerYear(timestart)+2000), ' (', datestr(timestamps(timestart), 'mm/dd/yyyy'), ')'])

axisC = axis;
rangex = axisC(2) - axisC(1);
text(axisC(2) - 0.2*rangex, -20, datestr(now));
clear axisC rangex

saveas(C, ['spd_', num2str(A), '.pdf'], 'pdf');

% 2_stick.pdf
D = numfigs*3 + A;

figure(D)
set(D, 'Position', [800 110 1100 850])
orient landscape

for s = 1:2
if d > 0
subplot(2,1,s)
plot([jday(timestart)-rangeD*.05 jday(timestop)+rangeD*.05], [0 0], 'k')
axis([jday(timestart)-rangeD*.05 jday(timestop)+rangeD*.05 -8 8])
% stickplots MUST HAVE a data aspect ratio of x=1 to y=1 (to z=1) for the
% sticks to point in the right direction.
daspect([1 1 1])
hold on

stickcolor = 'r';
for e = 1:length(daychg)-1
    stickhandle(e) = stickplot(jday(daychg(e):daychg(e+1)), fixedNorth(daychg(e):daychg(e+1),d)*.1, -fixedEast(daychg(e):daychg(e+1),d)*.1);
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
    badstick = stickplot(jday(bads(binbads, 1)), fixedNorth(bads(binbads, 1), d)*.1, -fixedEast(bads(binbads, 1), d)*.1);
    set(badstick, 'Color', [.8 .8 .8])
    clear bads binbads badstick
end

for j = {badmaxpitch, badmaxroll, badchgpitch, badchgroll}
    bads = j{1};
    badstick = stickplot(jday(bads), fixedNorth(bads, d)*.1, -fixedEast(bads, d)*.1);
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
xlabel(['Days since 00:00 01/01 ', num2str(SerYear(timestart)+2000), ' (', datestr(timestamps(timestart), 'mm/dd/yyyy'), ')'])
text(daxis(2) - (daxis(2)-daxis(1))*.2, -12, datestr(now));

clear daxis
saveas(D, ['stick_', num2str(A), '.pdf'], 'pdf');

close(A)
close(B)
close(C)
close(D)
end % for A, bins loop

clear A B C D a b c d p q r s t e

%% contour figures 

colorfig = 101;
for t = {fixedMag(timestart:timestop,1:lastbindex), fixedDIR(timestart:timestop,1:lastbindex), fixedEast(timestart:timestop,1:lastbindex), fixedNorth(timestart:timestop,1:lastbindex)}

figure(colorfig)

for E = 1:2
subplot(2,1,E)
pcolor(jday(timestart:timestop), dep(1:lastbindex)', t{1}')
shading interp
hold on
ylabel('Distance from bottom (m)')
ylim([2 20])
vecbar = colorbar;

if colorfig == 101
    ylabel(vecbar,'Speed (cm sec ^-^1)')
    caxis([0 80])
end
if colorfig ==102
    ylabel(vecbar,'Direction (degrees)')
    caxis([0 360])
end
if colorfig ==103
    ylabel(vecbar,'East U (cm sec ^-^1)')
    caxis([-80 80])
end
if colorfig ==104
    ylabel(vecbar,'North V (cm sec ^-^1)')
    caxis([-80 80])
end
clear vecbar

if E == 1
    set(gca, 'xticklabel', [])
end
if E == 2
    xlabel(['Days since 00:00 01/01 ', num2str(SerYear(timestart)+2000), ' (', datestr(timestamps(timestart), 'mm/dd/yyyy'), ')'])
    for g = {badM, badecho, badcorr, badPG}
        baddot = plot(jday(g{1}(:,1)), dep(g{1}(:,2))','.');
        set(baddot, 'color', [.8 .8 .8])
        clear baddot
    end
    for h = {badmaxpitch, badmaxroll, badchgpitch, badchgroll}
        bads = h{1};
        if ~isempty(bads)
            for P = 1:length(bads)
                badline = plot(jday(bads(P)), [0 dep(lastbindex)]);
                set(badline, 'color', [.8 .8 .8]);
                set(badline, 'linewidth', 2);
                clear badline
            end
        end
        clear bads P
    end
    text(jday(timestop), -3, datestr(now))
    clear g h
end

plot(botjday, botdepth, 'r')

end % E - one full plot

orient landscape
saveas(colorfig, ['color_', num2str(colorfig), '.pdf'], 'pdf');
colorfig = colorfig + 1;
clear E 
end % t - color plots
clear colorfig t

%% 5_zavg.pdf

for r = 1:length(DIR)
    getnansE = isnan(fixedEast(r,1:lastbindex));
    avgEast(r) = mean(fixedEast(r, getnansE==0));
    getnansN = isnan(fixedNorth(r,1:lastbindex));
    avgNorth(r) = mean(fixedNorth(r, getnansN==0));
    getnansV = isnan(fixedVert(r, 1:lastbindex));
    avgVert(r) = mean(fixedVert(r, getnansV==0));
    [avgth, avgMag(r)] = cart2pol(-avgNorth(r), -avgEast(r));
    avgDIR(r) = radtodeg(avgth)+180;
    clear getnans* avgth
end
clear r

overall = figure;
set(overall, 'Position', [800 110 1100 850])
orient landscape
subplot(4,1,1)
plot([jday(timestart) jday(timestop)], [0 0],'k')
axis([jday(timestart)-rangeD*.05 jday(timestop)+rangeD*.05 -4 4])
%stickplots MUST HAVE a data aspect ratio of x=1 to y=1 (to z=1) for the
%sticks to point in the right direction. All other daspect() in zavg figure
%are then scaled to be the same size as the stickplot. (Divide y by
%stickplot y-range.)
daspect([1 1 1])
hold on

stickcolor = 'r';
for e = 1:length(daychg)-1
    stickhandle(e) = stickplot(jday(daychg(e):daychg(e+1)), avgNorth(daychg(e):daychg(e+1))*.05, -avgEast(daychg(e):daychg(e+1))*.05);
    set(stickhandle(e), 'color', stickcolor);
    if stickcolor == 'r'
        stickcolor = 'g';
    elseif stickcolor == 'g'
        stickcolor = 'b';
    elseif stickcolor == 'b'
        stickcolor = 'r';
    end
end
set(gca, 'ytick', [-4 -2 0 2 4])
set(gca, 'yticklabel', [-80 -40 0 40 80])
ylabel('(cm s^-^1)')
set(gca, 'xticklabel', [])
title('up = W')

clear e stickcolor stickhandle

subplot(4,1,2)

plot(jday(timestart:timestop), avgEast(timestart:timestop),'b')
hold on
plot(jday(timestart:timestop), avgNorth(timestart:timestop), 'r')
plot([jday(timestart) jday(timestop)], [0 0],'k')
axis([jday(timestart)-rangeD*.05 jday(timestop)+rangeD*.05 -80 80])
daspect([1 20 1])
set(gca, 'ytick', [-80 -40 0 40 80])
ylabel('u: blue (cm s^-^1)')
set(gca, 'xticklabel', [])

subplot(4,1,3)

plot(jday(timestart:timestop), avgMag(timestart:timestop),'k')
axis([jday(timestart)-rangeD*.05 jday(timestop)+rangeD*.05 0 80])
daspect([1 10 1])
set(gca, 'ytick', [0 20 40 60 80])
ylabel('Speed (cm s^-^1)')
set(gca, 'xticklabel', [])

subplot(4,1,4)

plot(jday(timestart:timestop), avgDIR(timestart:timestop),'k')
axis([jday(timestart)-rangeD*.05 jday(timestop)+rangeD*.05 0 360])
daspect([1 45 1])
set(gca, 'ytick', [0 90 180 270 360])
ylabel('Direction')
xlabel(['Days since 00:00 01/01 ', num2str(SerYear(timestart)+2000), ' (', datestr(timestamps(timestart), 'mm/dd/yyyy'), ')'])

axisG = axis;
rangex = axisG(2) - axisG(1);
text(axisG(2) - 0.2*rangex, -150, datestr(now));
clear axisG rangex

saveas(overall, [datestr(timestamps(timestart), 'yyyymmdd'), '_', datestr(timestamps(timestop), 'yyyymmdd'), '_zavg.pdf'],'pdf');

%% average.dat

avgcsv = fopen([datestr(timestamps(timestart), 'yyyymmdd'), '_', datestr(timestamps(timestop), 'yyyymmdd'), '_zavg.csv'], 'w');
fprintf(avgcsv, 'Scan, Timestamp, avgEast_u, avgNorth_v, avgVert_w, avgSpeed, avgDirection\n');

avgdat = fopen([datestr(timestamps(timestart), 'yyyymmdd'), '_', datestr(timestamps(timestop), 'yyyymmdd'), '_zavg.dat'], 'w');
fprintf(avgdat, 'n,    y,    m,    d,    h,    min,    t_day,        tdep,    um (cm/s),    vm,    spd,    dir\n');

jlength = length(num2str(timestop)) + 2;
timelength = 6;
jdaylength = length(num2str(floor(max(jday)))) + 4;
trklength = length(num2str(floor(max(bottrack)))) + 4;
eastlength = length(num2str(floor(max(avgEast)))) + 4;
northlength = length(num2str(floor(max(avgNorth)))) + 4;
maglength = length(num2str(floor(max(avgMag)))) + 4;
dirlength = length(num2str(floor(max(avgDIR)))) + 4;

for j = timestart:timestop
    fprintf(avgcsv, '%i,', j);
    if ~isnan(timestamps(j))
        fprintf(avgcsv, '%s,', datestr(timestamps(j)));
    else
        fprintf(avgcsv, 'NaN,');
    end
    fprintf(avgcsv, '%4.2f,', avgEast(j));
    fprintf(avgcsv, '%4.2f,', avgNorth(j));
    fprintf(avgcsv, '%4.2f,', avgVert(j));
    fprintf(avgcsv, '%5.2f,', avgMag(j));
    fprintf(avgcsv, '%5.2f\n', avgDIR(j));

    qpad = jlength - length(num2str(j));
    for pad = 1:qpad
        fprintf(avgdat, ' ');
    end
    clear *pad
    fprintf(avgdat, '%i', j);

    [yrj moj dayj hrj mnj ~] = datevec(timestamps(j));
    qpad = timelength - length(num2str(yrj));
    for pad = 1:qpad
        fprintf(avgdat, ' ');
    end
    clear *pad
    fprintf(avgdat, '%i', yrj);
    qpad = timelength - length(num2str(moj));
    for pad = 1:qpad
        fprintf(avgdat, ' ');
    end
    clear *pad
    fprintf(avgdat, '%i', moj);
    qpad = timelength - length(num2str(dayj));
    for pad = 1:qpad
        fprintf(avgdat, ' ');
    end
    clear *pad
    fprintf(avgdat, '%i', dayj);
    qpad = timelength - length(num2str(hrj));
    for pad = 1:qpad
        fprintf(avgdat, ' ');
    end
    clear *pad
    fprintf(avgdat, '%i', hrj);
    qpad = timelength - length(num2str(mnj));
    for pad = 1:qpad
        fprintf(avgdat, ' ');
    end
    clear *pad
    fprintf(avgdat, '%i.5', mnj);
    clear yrj moj dayj hrj mnj

    qpad = timelength - length(num2str(jday(j)));
    for pad = 1:qpad
        fprintf(avgdat, ' ');
    end
    clear *pad
    fprintf(avgdat, '%10.6f', jday(j));

    qpad = timelength - length(num2str(bottrack(j)));
    for pad = 1:qpad
        fprintf(avgdat, ' ');
    end
    clear *pad
    fprintf(avgdat, '%6.3f', bottrack(j));

    qpad = timelength - length(num2str(avgEast(j)));
    for pad = 1:qpad
        fprintf(avgdat, ' ');
    end
    clear *pad
    fprintf(avgdat, '%4.2f', avgEast(j));

    qpad = timelength - length(num2str(avgNorth(j)));
    for pad = 1:qpad
        fprintf(avgdat, ' ');
    end
    clear *pad
    fprintf(avgdat, '%4.2f', avgNorth(j));

    qpad = timelength - length(num2str(avgMag(j)));
    for pad = 1:qpad
        fprintf(avgdat, ' ');
    end
    clear *pad
    fprintf(avgdat, '%5.2f', avgMag(j));

    qpad = timelength - length(num2str(avgDIR(j)));
    for pad = 1:qpad
        fprintf(avgdat, ' ');
    end
    clear *pad
    fprintf(avgdat, '%5.2f\n', avgDIR(j));
end
fclose(avgdat);
fclose(avgcsv);
clear j avgdat avgcsv *length *pad

%% for the Word doc

doclog = fopen([datestr(timestamps(timestart), 'yyyymmdd'), '_', datestr(timestamps(timestop), 'yyyymmdd'), '_worddoc.log'], 'w');

fprintf(doclog, 'ADCP\n\t: Export from 1 - %i, bin numbers %i - %i\n', length(SerYear), SerBins(1), SerBins(end));
fprintf(doclog, '\t: Process %i (%s = %10.6f) - %i (%s = %10.6f)\n', timestart, datestr(timestamps(timestart), 'HH:MM.5 on yymmdd'), jday(timestart), timestop, datestr(timestamps(timestop), 'HH:MM.5 on yymmdd'), jday(timestop));

if exist('firstbot', 'var')
    fprintf(doclog, '\t: From averaged bottom CTD, use %i (%s) - %i (%s)\n\n', firstbot, datestr(inbotstamps(firstbot), 'HH:MM.5 on yymmdd'), lastbot, datestr(inbotstamps(lastbot), 'HH:MM.5 on yymmdd'));
else
    fprintf(doclog, '\t: bottom CTD had no data in the same timeframe.\n\n');
end
fprintf(doclog, 'ADCP missing values filled in manually using linear interpolation in time:\n');

s = 0;
for n = timestart:timestop
    for m = 1:lastbindex
        if (isnan(East_u(n, m)) && ~isnan(fixedEast(n, m))) || (isnan(North_v(n, m)) && ~isnan(fixedNorth(n, m))) || (isnan(Mag(n, m)) && ~isnan(fixedMag(n, m))) || (isnan(DIR(n, m)) && ~isnan(fixedDIR(n, m))) && ~isnan(timestamps(n))
            fprintf(doclog, '%s, bin %i: ', datestr(timestamps(n), 'mm/dd/yyyy HH:MM.5'), SerBins(m));
            if isnan(East_u(n, m)) && ~isnan(fixedEast(n, m))
                fprintf(doclog, 'East = %4.2f, ', fixedEast(n, m));
                s = s + 1;
            end
            if isnan(North_v(n, m)) && ~isnan(fixedNorth(n, m))
                fprintf(doclog, 'North = %4.2f, ', fixedNorth(n, m));
                s = s + 1;
            end
            if isnan(Mag(n, m)) && ~isnan(fixedMag(n, m))
                fprintf(doclog, 'Magnitude = %4.2f, ', fixedMag(n, m));
                s = s + 1;
            end
            if isnan(DIR(n, m)) && ~isnan(fixedDIR(n, m))
                fprintf(doclog, 'Direction = %4.2f\n', fixedDIR(n, m));
                s = s + 1;
            else
                fprintf(doclog, '\n');
            end
        end
    end
end
fprintf(doclog, '%i total missing values.\n', s);
fclose(doclog);
clear s n m doclog ans textdata colheaders data numfigs overall rangeD

save([datestr(timestamps(timestart), 'yyyymmdd'), '_', datestr(timestamps(timestop), 'yyyymmdd'), '_ADCP-processed.mat'])
