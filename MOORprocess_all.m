%% This script to process batches of mooring data using Matlab 2010b.
% import moor-timestamps.txt before starting.
%
% !!! check that CTD/YSI variable columns are correct !!! They are not
% always in the same order, and not every variable will exist in every
% dataset. 
% YSI: MOORprocess_all.m lines 77-95
% CTDSURF/CTDBOT: lines 251-259
% CTDcast: lines 1255-1275

stations = textdata(2:end,1);
startdates = textdata(2:end,2);
starttimes = textdata(2:end,3);
stationname = textdata(2:end,4);
scanstarts = data(:,1);
scanstops = data(:,2);

%This tells matlab which line of timestamps.txt refers to which instrument.
% Be sure the names in timestamp.txt are correct and exactly as below.
inputT = cell(1,5);
inputctd = cell(1,3);
for l = 1:length(stationname)
    if strcmp(stationname(l), 'T1') 
        inputT{1} = l;
    end
    if strcmp(stationname(l), 'T2') 
        inputT{2} = l;
    end
    if strcmp(stationname(l), 'T3') 
        inputT{3} = l;
    end
    if strcmp(stationname(l), 'T4') 
        inputT{4} = l;
    end
    if strcmp(stationname(l), 'T5') 
        inputT{5} = l;
    end
    if strcmp(stationname(l), 'TD')
        inputTD = l;
    end
    if strcmp(stationname(l), 'YSI') 
        inputctd{3} = l;
    end
    if strcmp(stationname(l), 'CTDBOT') 
        inputctd{1} = l;
    end
    if strcmp(stationname(l), 'CTDSURF') 
        inputctd{2} = l;
    end
    if strcmp(stationname(l), 'cast')
        inputcast = l;
    end
end

% These initialize the cell arrays that will hold the data from the two
% CTDs and the YSI.
ctdburst = cell(1,3);
ctdburststamp = cell(1,3);
ctdstart =  cell(1,3);
ctdstop =  cell(1,3);
ctdscan =  cell(1,3);
ctdtimestamp =  cell(1,3);
ctddepth =  cell(1,3);
ctdtemp =  cell(1,3);
ctdsal =  cell(1,3);
ctddens =  cell(1,3);
calcdens = cell(1,3);
ctdpress =  cell(1,3);
ctdcond =  cell(1,3);
ctdspcond =  cell(1,3);
ctdsound =  cell(1,3);
% This will hold replaced bad values later reported into the MS-Word doc. 
interpnums = [];

%% YSI data: = cell 3
if ~isempty(inputctd{3})
ctdstart{3} = scanstarts(inputctd{3});
ctdstop{3} = scanstops(inputctd{3});
ysistation = char(stations(inputctd{3}));
YSI = importdata([ysistation,'-matlab.txt']);
yr = YSI(:,3);
mo = YSI(:,1);
day = YSI(:,2);
hr = YSI(:,4);
mn = YSI(:,5);
ctdtimestamp{3}=datenum(yr+2000,mo,day,hr,mn,0) - datenum(0,0,0,0,2,30); %subtract 2.5 min to line up with everything else
ctddepth{3} = YSI(:,6);
ctdsal{3} = YSI(:,7);
ctdtemp{3} = YSI(:,8);
%ctdpress{3} = YSI(:,9);
ctdpress{3}(1:length(ctdtimestamp{3})) = NaN;
ctdcond{3} = YSI(:,9);
ctdspcond{3} = YSI(:,10);
ctdoxy{3} = YSI(:,12);
ctdoxysat{3} = YSI(:,13);
%ctdchla{3} = YSI(:,?);
ctdturb{3} = YSI(:,14);
%ctdturb{3}(1:length(ctdtimestamp{3})) = NaN; % a kludge for nonexistent variables in an existing package
ctdtds{3} = YSI(:,15);
ctdscan{3} = 1:length(ctdtimestamp{3});

%calculate day of year from timestamp
for i = 1:length(ctdtimestamp{3})
    excess = datenum(yr(i)+2000, 0, 1, 0, 0, 0); % Jan 1 is Day 0
    jdayYSI(i) = ctdtimestamp{3}(i) - excess;
    clear excess
end

%density (sigma-t) calculated from temperature and salinity
%equation from Jim Manning, http://globec.whoi.edu/globec-dir/sigmat-calc-matlab.html
t=ctdtemp{3};
s=ctdsal{3};
a=[6.536332e-9,-1.120083e-6,1.001685e-4,-9.095290e-3,6.793952e-2,999.842594];
b=[5.3875e-9,-8.2467e-7,7.6438e-5,-4.0899e-3,8.24493e-1];
c=[-1.6546e-6,1.0227e-4,-5.72466e-3];
d=4.8314e-4;
calcdens{3}=(polyval(a,t)+(polyval(b,t)+polyval(c,t).*sqrt(s)+d*s).*s) - 1000;
clear t s a b c d 

fiddat = fopen([ysistation, '.dat'], 'w');
fprintf(fiddat, 'n,      jday,    T(C),    s(ppt),    den,    dep(m),    O2(mg/l), O2s(%%), %s (17.00 m): t0 = %4.2f (%s)\n', ysistation, ctdstart{3}, datestr(ctdtimestamp{3}(ctdstart{3}), 'HH:MM mm/dd/yyyy'));
qlength = length(num2str(ctdstop{3})) + 2;
jdaylength = length(num2str(floor(max(jdayYSI)))) + 4;
templength = length(num2str(floor(max(ctdtemp{3})))) + 4;
sallength = length(num2str(floor(max(ctdsal{3})))) + 4;
denslength = length(num2str(floor(max(calcdens{3})))) + 4;
depthlength = length(num2str(floor(max(ctddepth{3})))) + 4;
oxylength = length(num2str(floor(max(ctdoxy{3})))) + 4;
oxysatlength = length(num2str(floor(max(ctdoxysat{3})))) + 4;
timelength = 6;

fidrange = fopen([ysistation, '_goodrange.csv'],'w');
fprintf(fidrange, 'Scan,Timestamp,Depth,Sal,Temp,Dens,Press,Cond,Sp.Cond,Oxy,Oxysat,Turbidity,TDS\n');
for q = ctdstart{3}:ctdstop{3}
    qpad = qlength - length(num2str(q));
    for pad = 1:qpad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%i', q);
    
    jdaypad = jdaylength - length(num2str(floor(jdayYSI(q))));
    for pad = 1:jdaypad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%10.6f', jdayYSI(q)); %jday

    temppad = templength - length(num2str(floor(ctdtemp{3}(q))));
    for pad = 1:temppad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%6.4f', ctdtemp{3}(q));

    salpad = sallength - length(num2str(floor(ctdsal{3}(q))));
    for pad = 1:salpad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%6.4f', ctdsal{3}(q));

    denspad = denslength - length(num2str(floor(calcdens{3}(q))));
    for pad = 1:denspad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%6.4f', calcdens{3}(q));

    depthpad = depthlength - length(num2str(floor(ctddepth{3}(q))));
    for pad = 1:depthpad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%6.4f', ctddepth{3}(q));

    oxypad = oxylength - length(num2str(floor(ctdoxy{3}(q))));
    for pad = 1:oxypad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%6.4f', ctdoxy{3}(q));

    oxysatpad = oxysatlength - length(num2str(floor(ctdoxysat{3}(q))));
    for pad = 1:oxysatpad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%6.4f', ctdoxysat{3}(q));

    fprintf(fiddat, '    %i', yr(q)+2000);

    timepad = timelength - length(num2str(mo(q)));
    for pad = 1:timepad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%i', mo(q));

    timepad = timelength - length(num2str(day(q)));
    for pad = 1:timepad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%i', day(q));

    timepad = timelength - length(num2str(hr(q)));
    for pad = 1:timepad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%i', hr(q));

    timepad = timelength - length(num2str(mn(q)));
    for pad = 1:timepad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%4.2f\n', mn(q)-2.5);

    fprintf(fidrange, '%i,', q);
    if ~isnan(ctdtimestamp{3}(q))
        fprintf(fidrange, '%s,', datestr(ctdtimestamp{3}(q)));
    else
        fprintf(fidrange, 'NaN,');
    end
    fprintf(fidrange, '%6.4f,', ctddepth{3}(q));
    fprintf(fidrange, '%6.4f,', ctdsal{3}(q));
    fprintf(fidrange, '%6.4f,', ctdtemp{3}(q));
    fprintf(fidrange, '%6.4f,', calcdens{3}(q));
    fprintf(fidrange, '%6.4f,', ctdpress{3}(q));
    fprintf(fidrange, '%6.4f,', ctdcond{3}(q));
    fprintf(fidrange, '%6.4f,', ctdspcond{3}(q));
    fprintf(fidrange, '%6.4f,', ctdoxy{3}(q));
    fprintf(fidrange, '%6.4f,', ctdoxysat{3}(q));
%    fprintf(fidrange, '%6.4f,', ctdchla{3}(q));
    fprintf(fidrange, '%6.4f,', ctdturb{3}(q));
    fprintf(fidrange, '%6.4f\n', ctdtds{3}(q));
end
fclose(fidrange);
fclose(fiddat);

clear yr mo day hr mn i q fidrange fiddat *length *pad
end %end YSI

%% CTDs: BOT = 1, SURF = 2
ctdmfromsurf = [20 4 17]; %depths below surface for bot, surf, ysi

for r = 1:2
if ~isempty(inputctd{r})
ctdstart{r} = scanstarts(inputctd{r});
ctdstop{r} = scanstops(inputctd{r});
ctdstation = char(stations(inputctd{r}));
CTD = importdata([ctdstation, '-matlab.txt']);
ctdtemp{r} = CTD(:,1);
ctdcond{r} = CTD(:,2);
ctdpress{r} = CTD(:,3);
ctdsal{r} = CTD(:,4);
ctdsound{r} = CTD(:,5);
ctdjday = CTD(:,6); %coming in, jday is Jan 1 = Day 1
ctddens{r} = CTD(:,8);
ctddepth{r} = CTD(:,9);
ctdspcond{r} = CTD(:,10);

% interpolate through any NaNs
fixed = cell(1,9);
fixed{1} = ctdsal{r};
fixed{2} = ctdtemp{r};
fixed{3} = ctddens{r};
fixed{4} = ctdpress{r};
fixed{5} = ctdcond{r};
fixed{6} = ctdspcond{r};
fixed{7} = ctdsound{r};
fixed{8} = ctddepth{r};
fixed{9} = ctdjday;

for w = 1:9
    blort = fixed{w};
    blanks = find(isnan(blort));
    if ~isempty(blanks)
        for e = 1:length(blanks)
            m = blanks(e);
            while isnan(blort(m))
                if m > 1
                    m = m - 1;
                else
                    break
                end
            end
            n = blanks(e);
            while isnan(blort(n))
                if n < length(blort)
                    n = n + 1;
                else
                    break
                end
            end
            blort(blanks(e)) = interp1([m n], [blort(m) blort(n)], blanks(e));
            interpnums = [interpnums; r w blanks(e) blort(blanks(e))];
        end
        clear e m n
    end
    fixed{w} = blort;
    clear blort blanks
end
clear w

ctdsal{r} = fixed{1};
ctdtemp{r} = fixed{2};
ctddens{r} = fixed{3};
ctdpress{r} = fixed{4};
ctdcond{r} = fixed{5};
ctdspcond{r} = fixed{6};
ctdsound{r} = fixed{7};
ctddepth{r} = fixed{8};
ctdjday = fixed{9};
clear fixed

% time
[yr, ~, ~, ~, ~, ~] = datevec(startdates(inputctd{r}));
ctdtimestamp{r} = ctdjday + datenum(yr, 0, 0); 
ctdscan{r} = 1:length(ctddepth{r});

%%equation from Jim Manning, http://globec.whoi.edu/globec-dir/sigmat-calc-matlab.html
t=ctdtemp{r};
s=ctdsal{r};
a=[6.536332e-9,-1.120083e-6,1.001685e-4,-9.095290e-3,6.793952e-2,999.842594];
b=[5.3875e-9,-8.2467e-7,7.6438e-5,-4.0899e-3,8.24493e-1];
c=[-1.6546e-6,1.0227e-4,-5.72466e-3];
d=4.8314e-4;
calcdens{r}=(polyval(a,t)+(polyval(b,t)+polyval(c,t).*sqrt(s)+d*s).*s) - 1000;
clear t s a b c d 

%bin averaging of the 20-minute intervals
[firstest, lastest] = moorburst(ctdstart{r}, ctdtimestamp{r}(ctdstart{r}), ctdstop{r}, ctdtimestamp{r}(ctdstop{r}));

k = 1;
p = firstest;
while p < lastest - 18
        ctdburst{r}(k,1) = p;
        ctdburst{r}(k,2) = mean(ctdtimestamp{r}(p:p+19));
        ctdburst{r}(k,3) = mean(ctddepth{r}(p:p+19));
        ctdburst{r}(k,4) = mean(ctdsal{r}(p:p+19));
        ctdburst{r}(k,5) = mean(ctdtemp{r}(p:p+19));
        ctdburst{r}(k,6) = mean(calcdens{r}(p:p+19));
        ctdburst{r}(k,7) = mean(ctdpress{r}(p:p+19));
        ctdburst{r}(k,8) = mean(ctdcond{r}(p:p+19));
        ctdburst{r}(k,9) = mean(ctdspcond{r}(p:p+19));
        ctdburst{r}(k,10) = mean(ctdsound{r}(p:p+19));
        ctdburst{r}(k,11) = mean(ctdjday(p:p+19)); % Jan 1 is still Day 1 
        ctdburst{r}(k,12) = mean(ctddens{r}(p:p+19));
    p = p + 20;
    k = k + 1;
end
ctdburststamp{r} = ctdburst{r}(:,2);
ctdstart{r} = firstest;
ctdstop{r} = lastest;

%text file outputs
fidrange = fopen([ctdstation, '_goodrange.csv'],'w');
fprintf(fidrange, 'Scan,Timestamp,Depth,Sal,Temp,Calc Dens,SBE Dens,Press,Cond,SpCond,Sound\n');
for q = ctdstart{r}:ctdstop{r}
    fprintf(fidrange, '%i,', ctdscan{r}(q));
    if ~isnan(ctdtimestamp{r}(q))
        fprintf(fidrange, '%s,', datestr(ctdtimestamp{r}(q)));
    else
        fprintf(fidrange, 'NaN,');
    end
    fprintf(fidrange, '%6.4f,', ctddepth{r}(q));
    fprintf(fidrange, '%6.4f,', ctdsal{r}(q));
    fprintf(fidrange, '%6.4f,', ctdtemp{r}(q));
    fprintf(fidrange, '%6.4f,', calcdens{r}(q));
    fprintf(fidrange, '%6.4f,', ctddens{r}(q));
    fprintf(fidrange, '%6.4f,', ctdpress{r}(q));
    fprintf(fidrange, '%6.4f,', ctdcond{r}(q));
    fprintf(fidrange, '%6.4f,', ctdspcond{r}(q));
    fprintf(fidrange, '%6.4f\n', ctdsound{r}(q));   
end
fclose(fidrange);

fiddat = fopen([ctdstation, '.dat'],'w');
fprintf(fiddat, '    n,    jday,            T(C),    s(ppt),    den,      dep(m),    %s (%i.00 m): t0 = %4.2f (%s)\n', ctdstation, ctdmfromsurf(r), floor(ctdburst{r}(1,11)-1), datestr(floor(ctdburst{r}(1,2)), 'HH:MM mm/dd/yyyy'));

qlength = length(num2str(length(ctdburst{r}))) + 2;
jdaylength = length(num2str(floor(max(ctdburst{r}(:,11))))) + 4;
depthlength = length(num2str(floor(max(ctdburst{r}(:,3))))) + 4;
sallength = length(num2str(floor(max(ctdburst{r}(:,4))))) + 4;
templength = length(num2str(floor(max(ctdburst{r}(:,5))))) + 4;
denslength = length(num2str(floor(max(ctdburst{r}(:,6))))) + 4;

fidest = fopen([ctdstation, '_est.csv'],'w');
fprintf(fidest, 'First Scan,Timestamp,Depth,Sal,Temp,Calc Dens,SBE Dens,Press,Cond,SpCond,Sound\n');
for s = 1:length(ctdburst{r})
        qpad = qlength - length(num2str(s));
    for pad = 1:qpad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%i', s);

    qpad = jdaylength - length(num2str(floor(ctdburst{r}(s, 11)-1)));
    for pad = 1:qpad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%10.6f', ctdburst{r}(s, 11)-1); %jday

    qpad = templength - length(num2str(floor(ctdburst{r}(s, 5))));
    for pad = 1:qpad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%6.4f', ctdburst{r}(s, 5)); %temp

    qpad = sallength - length(num2str(floor(ctdburst{r}(s, 4))));
    for pad = 1:qpad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%6.4f', ctdburst{r}(s, 4)); %sal

    qpad = denslength - length(num2str(floor(ctdburst{r}(s, 6))));
    for pad = 1:qpad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%6.4f', ctdburst{r}(s, 6)); %calcdens

    qpad = depthlength - length(num2str(floor(ctdburst{r}(s, 3))));
    for pad = 1:qpad
        fprintf(fiddat, ' ');
    end
    clear *pad
    fprintf(fiddat, '%6.4f', ctdburst{r}(s, 3)); %depth

    fprintf(fidest, '%i,', ctdburst{r}(s,1)); %scan
    if ~isnan(ctdburst{r}(s, 2))
        fprintf(fiddat, '    %s\n', datestr(ctdburst{r}(s, 2), 'mm dd HH MM SS')); %timestamp
        fprintf(fidest, '%s,', datestr(ctdburst{r}(s,2))); %timestamp
    else
        fprintf(fiddat, '    NaN\n');
        fprintf(fidest, 'NaN,');
    end
    fprintf(fidest, '%6.4f,', ctdburst{r}(s,3)); %depth
    fprintf(fidest, '%6.4f,', ctdburst{r}(s,4)); %sal
    fprintf(fidest, '%6.4f,', ctdburst{r}(s,5)); %temp
    fprintf(fidest, '%6.4f,', ctdburst{r}(s,6)); %calcdens
    fprintf(fidest, '%6.4f,', ctdburst{r}(s,12)); %ctddens
    fprintf(fidest, '%6.4f,', ctdburst{r}(s,7)); %press
    fprintf(fidest, '%6.4f,', ctdburst{r}(s,8)); %cond
    fprintf(fidest, '%6.4f,', ctdburst{r}(s,9)); %spcond
    fprintf(fidest, '%6.4f\n', ctdburst{r}(s,10)); %sound
end
fclose(fidest);
fclose(fiddat);

clear k p yr firstest lastest fidest fiddat s fidrange q ctdstation CTD *length *pad
end %end CTD if exist
end %end CTD loop
clear r

%% thermistors
Tmfromsurf = [6 9 11 13 18]; %depths below surface for each, from T1 to T5

timeT = cell(1,5);
tempT = cell(1,5);
startT = cell(1,5);
stopT = cell(1,5);
burstT = cell(1,5);
burstTstamp = cell(1,5);

for r = 1:length(inputT)
if ~isempty(inputT{r})
stationT = char(stations(inputT{r}));
Tname = [stationT, '-matlab.txt'];
thermistor = importdata(Tname);
temptextT = getfield(thermistor, 'textdata');
dataT = getfield(thermistor, 'data');
tempT{r} = dataT(:,1);
[~,b] = size(dataT);
if b > 1
    pressT = dataT(:,2);
    cnv = importdata([char(stations(inputTD)), 'cnv-matlab.txt']);
    depthT = cnv(:,2);
end

startT{r} = scanstarts(inputT{r});
stopT{r} = scanstops(inputT{r});

[yr, mo, day, ~, ~, ~] = datevec(temptextT(:,1));
[~, ~, ~, hr, mn, sec] = datevec(temptextT(:,2));
timeT{r} = datenum(yr, mo, day, hr, mn, sec);
starttimeT = timeT{r}(startT{r});
stoptimeT = timeT{r}(stopT{r});

for i = 1:length(timeT{r})
    excess = datenum(yr(i), 0, 1, 0, 0, 0); % Jan 1 is Day 0
    jdayT(i) = timeT{r}(i) - excess;
    clear excess
end

clear yr mo day hr mn sec i b

[firstest, lastest] = moorburst(startT{r}, starttimeT, stopT{r}, stoptimeT);

k = 1;
p = firstest;
while p < lastest-18
        burstT{r}(k,1) = p;
        burstT{r}(k,2) = mean(timeT{r}(p:p+19));
        burstT{r}(k,3) = mean(jdayT(p:p+19));
        burstT{r}(k,4) = mean(tempT{r}(p:p+19));
        if exist('pressT','var') ==1
            burstT{r}(k,5) = mean(pressT(p:p+19));
            burstT{r}(k,6) = mean(depthT(p:p+19));
        end
    p = p + 20;
    k = k +1;
end
clear p k
burstTstamp{r} = burstT{r}(:,2);
startT{r} = firstest;
stopT{r} = lastest;

fidT = fopen([stationT, '.dat'], 'w');
fprintf(fidT, '    n,     jday,        t(C),      T01 %s. ( %i.00 m): t0 = %4.2f (%s)\n', datestr(burstT{r}(end,2), 'mmddyy'), Tmfromsurf(r), floor(jdayT(1)), datestr(floor(burstT{r}(1,2)), 'HH:MM mm/dd/yyyy'));

qlength = length(num2str(length(burstT{r}))) + 2;
jdaylength = length(num2str(floor(max(burstT{r}(:,3))))) + 4;
templength = length(num2str(floor(max(burstT{r}(:,4))))) + 4;
if exist('pressT','var') ==1
    depthlength = length(num2str(floor(max(burstT{r}(:,6))))) + 4;
end

fidest = fopen([stationT, '_est.csv'], 'w');
fprintf(fidest, 'first scan,avg timestamp,jday, temp, press,depth\n');

for j = 1:length(burstT{r})
    qpad = qlength - length(num2str(j));
    for pad = 1:qpad
        fprintf(fidT, ' ');
    end
    clear *pad
    fprintf(fidT, '%i', j); %line number

    qpad = jdaylength - length(num2str(floor(burstT{r}(j,3))));
    for pad = 1:qpad
        fprintf(fidT, ' ');
    end
    clear *pad
    fprintf(fidT, '%10.6f', burstT{r}(j,3)); %jday

    qpad = templength - length(num2str(floor(burstT{r}(j,4))));
    for pad = 1:qpad
        fprintf(fidT, ' ');
    end
    clear *pad
    fprintf(fidT, '%6.4f', burstT{r}(j, 4)); %temp
    if exist('pressT','var') ==1
        qpad = depthlength - length(num2str(floor(burstT{r}(j,6))));
        for pad = 1:qpad
            fprintf(fidT, ' ');
        end
        clear *pad
        fprintf(fidT, '%6.4f', burstT{r}(j,6)); %depth
    end
    fprintf(fidest, '%i,', burstT{r}(j,1)); %first scan
    if ~isnan(burstT{r}(j,2))
        fprintf(fidT, '    %s\n', datestr(burstT{r}(j, 2), 'mm dd HH MM SS')); %timestamp
        fprintf(fidest, '%s,', datestr(burstT{r}(j,2))); %timestamp
    else
        fprintf(fidT, '    NaN\n');
        fprintf(fidest, 'NaN,');
    end
    fprintf(fidest, '%10.6f,', burstT{r}(j, 3)); %jday
    fprintf(fidest, '%6.4f', burstT{r}(j,4)); %temp
    if exist('pressT','var') ==1
        fprintf(fidest, ',%6.4f', burstT{r}(j,5)); %press
        fprintf(fidest, ',%6.4f\n', burstT{r}(j,6)); %depth
    else
        fprintf(fidest, '\n');
    end
end
fclose(fidT);
fclose(fidest);
clear j *pad *length

fidgood = fopen([stationT, '_goodrange.csv'], 'w');
fprintf(fidgood, 'scan, timestamp, temp, press\n');
for j = startT{r}:stopT{r}
    fprintf(fidgood, '%i,', j);
    if ~isnan(timeT{r}(j))
        fprintf(fidgood, '%s,', datestr(timeT{r}(j)));
    else
        fprintf(fidgood, 'NaN,');
    end
    fprintf(fidgood, '%6.4f', tempT{r}(j));
    if exist('pressT','var') ==1
        fprintf(fidgood, ',%6.4f', pressT(j));
        fprintf(fidgood, ',%6.4f\n', depthT(j));
    else
        fprintf(fidgood,'\n');
    end
end
fclose(fidgood);

clear stationT Tname thermistor dataT temptextT starttimeT stoptimeT jdayT firstest lastest fidT fidest fidgood j

end % end if exists
end % end thermistor for loop
clear r

%% Figures

%cellfun() does not work for min() or max()
%mintime = min([min(ctdtimestamp{1}), min(ctdtimestamp{2}), min(ctdtimestamp{3}), min(timeT{1}) min(timeT{2}) min(timeT{3}) min(timeT{4}) min(timeT{5})]);
%maxtime = max([max(ctdtimestamp{1}), max(ctdtimestamp{2}), max(ctdtimestamp{3}), max(timeT{1}) max(timeT{2}) max(timeT{3}) max(timeT{4}) max(timeT{5})]);

onehour = datenum(0,0,0,1,0,0);
%thermistor axes
ystart = floor(min([min(tempT{1}(startT{1}:stopT{1})) min(tempT{2}(startT{2}:stopT{2})) min(tempT{3}(startT{3}:stopT{3})) min(tempT{4}(startT{4}:stopT{4})) min(tempT{5}(startT{5}:stopT{5}))]));
yend = ceil(max([max(tempT{1}(startT{1}:stopT{1})) max(tempT{2}(startT{2}:stopT{2})) max(tempT{3}(startT{3}:stopT{3})) max(tempT{4}(startT{4}:stopT{4})) max(tempT{5}(startT{5}:stopT{5}))]));
minx = min([timeT{1}(startT{1}) timeT{2}(startT{2}) timeT{3}(startT{3}) timeT{4}(startT{4}) timeT{5}(startT{5})]);
maxx = max([timeT{1}(stopT{1}) timeT{2}(stopT{2}) timeT{3}(stopT{3}) timeT{4}(stopT{4}) timeT{5}(stopT{5})]);

%CTD x-axes
minstart = min([ctdtimestamp{1}(ctdstart{1}), ctdtimestamp{2}(ctdstart{2}), ctdtimestamp{3}(ctdstart{3})]);
maxstop = max([ctdtimestamp{1}(ctdstop{1}), ctdtimestamp{2}(ctdstop{2}), ctdtimestamp{3}(ctdstop{3})]);

mintime = min([minx minstart]);
maxtime = max([maxx maxstop]);

%total timeframe for all instruments
Tframestart = datestr(floor(mintime), 'yyyymmdd');
Tframeend = datestr(floor(maxtime),'yyyymmdd');

%% Thermistor Figures

%thermistors all
figure(1)
if ~isempty(inputT{1}) %black dots, gray line
    T1line = plot(timeT{1}, tempT{1});
    set(T1line, 'Color', [.8 .8 .8])
    hold on
    plot(burstT{1}(:,2), burstT{1}(:,4),'.k');
end
if ~isempty(inputT{2}) %dark green dots, green line
    plot(timeT{2}, tempT{2},'g')
    hold on
    T2dot = plot(burstT{2}(:,2),burstT{2}(:,4),'.');
    set(T2dot, 'color', [0 .75 0])
end
if ~isempty(inputT{3}) %magenta dots, pink line
    T3line = plot(timeT{3}, tempT{3});
    set(T3line, 'color', [1 .6 .78])
    hold on
    T3dot = plot(burstT{3}(:,2),burstT{3}(:,4),'.m');
end
if ~isempty(inputT{4}) %orange dots, yellow line
    plot(timeT{4}, tempT{4},'y')
    hold on
    T4dot = plot(burstT{4}(:,2), burstT{4}(:,4),'.');
    set(T4dot, 'color', [1 .69 .39]);
end
if ~isempty(inputT{5}) %blue dots, cyan line
    plot(timeT{5}, tempT{5},'c')
    hold on
    plot(burstT{5}(:,2), burstT{5}(:,4),'.b')
end
for r = 1:5
    if ~isempty(inputT{r})
        rcirc = plot(timeT{r}(startT{r}), tempT{r}(startT{r}), 'ro');
        set(rcirc, 'linewidth', 2);
        clear rcirc
        rcirc = plot(timeT{r}(stopT{r}), tempT{r}(stopT{r}), 'ro');
        set(rcirc, 'linewidth', 2);
        clear rcirc
    end
end
orient landscape
ylabel('Temperature ( \circC)')
title([Tframestart, '\_', Tframeend, ' good ranges (dates): T1 = ', datestr(timeT{1}(startT{1})), ' - ', datestr(timeT{1}(stopT{1})), ', T2 = ', datestr(timeT{2}(startT{2})), ' - ', datestr(timeT{2}(stopT{2})),sprintf(',\n'), 'T3 = ', datestr(timeT{3}(startT{3})), ' - ', datestr(timeT{3}(stopT{3})), ', T4 = ', datestr(timeT{4}(startT{4})), ' - ', datestr(timeT{4}(stopT{4})), ', T5 = ', datestr(timeT{5}(startT{5})), ' - ', datestr(timeT{5}(stopT{5}))])
xlabel('date')
axis([minx-onehour.*24 maxx+onehour.*24 ystart yend])

a = axis;
rangex = a(2)-a(1);
rangey = a(4)-a(3);
allticks = [minx, minx + rangex/4, minx + rangex/2, maxx - rangex/4, maxx];
set(gca,'xtick', allticks);
set(gca,'xticklabel', datestr(allticks));
set(gca,'xgrid','on')

spot1 = [a(2)+0.01*rangex a(4)-0.3*rangey];
leg1 = text(spot1(1), spot1(2), 'T1');
set(leg1, 'fontweight', 'bold');
spot2 = [a(2)+0.01*rangex a(4)-0.25*rangey];
leg2 = text(spot2(1), spot2(2), 'T2');
set(leg2, 'Color', [0 .75 0]);
set(leg2, 'fontweight', 'bold');
spot3 = [a(2)+0.01*rangex a(4)-0.2*rangey];
leg3 = text(spot3(1), spot3(2), 'T3');
set(leg3, 'Color', [1 0 1]);
set(leg3, 'fontweight', 'bold');
spot4 = [a(2)+0.01*rangex a(4)-0.15*rangey];
leg4 = text(spot4(1), spot4(2), 'T4');
set(leg4, 'Color', [1 .69 .39]);
set(leg4, 'fontweight', 'bold');
spot5 = [a(2)+0.01*rangex a(4)-0.1*rangey];
leg5 = text(spot5(1), spot5(2), 'T5');
set(leg5, 'Color', [0 0 1]);
set(leg5, 'fontweight', 'bold');
idspot = [a(2)-0.2*rangex a(3)-0.1*rangey];
text(idspot(1), idspot(2), datestr(now));

clear a range* *spot* leg* allticks

% thermistors beginning and end
figure(2)
subplot(2,1,1)
if ~isempty(inputT{1}) %black dots, gray line
    T1line = plot(timeT{1}, tempT{1});
    set(T1line, 'Color', [.8 .8 .8])
    hold on
    plot(burstT{1}(:,2), burstT{1}(:,4),'.k');
end
if ~isempty(inputT{2}) %dark green dots, green line
    plot(timeT{2}, tempT{2},'g')
    hold on
    T2dot = plot(burstT{2}(:,2),burstT{2}(:,4),'.');
    set(T2dot, 'color', [0 .75 0])
end
if ~isempty(inputT{3}) %magenta dots, pink line
    T3line = plot(timeT{3}, tempT{3});
    set(T3line, 'color', [1 .6 .78])
    hold on
    T3dot = plot(burstT{3}(:,2),burstT{3}(:,4),'.m');
end
if ~isempty(inputT{4}) %orange dots, yellow line
    plot(timeT{4}, tempT{4},'y')
    hold on
    T4dot = plot(burstT{4}(:,2), burstT{4}(:,4),'.');
    set(T4dot, 'color', [1 .69 .39]);
end
if ~isempty(inputT{5}) %blue dots, cyan line
    plot(timeT{5}, tempT{5},'c')
    hold on
    plot(burstT{5}(:,2), burstT{5}(:,4),'.b')
end
for r = 1:5
    if ~isempty(inputT{r})
        plot(timeT{r}(startT{r}), tempT{r}(startT{r}), 'ro');
    end
end
orient landscape
ylabel('Temperature ( \circC)')
axis([minx-onehour minx+onehour ystart yend])
title([Tframestart, '\_', Tframeend, ' good ranges (scan):', 'T1 = ', num2str(startT{1}), ' - ', num2str(stopT{1}), ', T2 = ', num2str(startT{2}), ' - ', num2str(stopT{2}), ', T3 = ',num2str(startT{3}), ' - ', num2str(stopT{3}), ', T4 = ', num2str(startT{4}), ' - ', num2str(stopT{4}), ', T5 = ', num2str(startT{5}), ' - ', num2str(stopT{5})])

a = axis;
rangex = a(2)-a(1);
rangey = a(4)-a(3);

startticks = [minx-onehour, minx-onehour + rangex/4, minx-onehour + rangex/2, minx+onehour - rangex/4, minx+onehour];
set(gca,'xtick', startticks);
set(gca,'xticklabel', datestr(startticks));
set(gca,'xgrid', 'on');

spot1 = [a(2)+0.01*rangex a(4)-0.3*rangey];
leg1 = text(spot1(1), spot1(2), 'T1');
set(leg1, 'fontweight', 'bold');
spot2 = [a(2)+0.01*rangex a(4)-0.25*rangey];
leg2 = text(spot2(1), spot2(2), 'T2');
set(leg2, 'Color', [0 .75 0]);
set(leg2, 'fontweight', 'bold');
spot3 = [a(2)+0.01*rangex a(4)-0.2*rangey];
leg3 = text(spot3(1), spot3(2), 'T3');
set(leg3, 'Color', [1 0 1]);
set(leg3, 'fontweight', 'bold');
spot4 = [a(2)+0.01*rangex a(4)-0.15*rangey];
leg4 = text(spot4(1), spot4(2), 'T4');
set(leg4, 'Color', [1 .69 .39]);
set(leg4, 'fontweight', 'bold');
spot5 = [a(2)+0.01*rangex a(4)-0.1*rangey];
leg5 = text(spot5(1), spot5(2), 'T5');
set(leg5, 'Color', [0 0 1]);
set(leg5, 'fontweight', 'bold');

clear range* a spot* leg*

subplot(2,1,2)
if ~isempty(inputT{1}) %black dots, gray line
    T1line = plot(timeT{1}, tempT{1});
    set(T1line, 'Color', [.8 .8 .8])
    hold on
    plot(burstT{1}(:,2), burstT{1}(:,4),'.k');
end
if ~isempty(inputT{2}) %dark green dots, green line
    plot(timeT{2}, tempT{2},'g')
    hold on
    T2dot = plot(burstT{2}(:,2),burstT{2}(:,4),'.');
    set(T2dot, 'color', [0 .75 0])
end
if ~isempty(inputT{3}) %magenta dots, pink line
    T3line = plot(timeT{3}, tempT{3});
    set(T3line, 'color', [1 .6 .78])
    hold on
    T3dot = plot(burstT{3}(:,2),burstT{3}(:,4),'.m');
end
if ~isempty(inputT{4}) %orange dots, yellow line
    plot(timeT{4}, tempT{4},'y')
    hold on
    T4dot = plot(burstT{4}(:,2), burstT{4}(:,4),'.');
    set(T4dot, 'color', [1 .69 .39]);
end
if ~isempty(inputT{5}) %blue dots, cyan line
    plot(timeT{5}, tempT{5},'c')
    hold on
    plot(burstT{5}(:,2), burstT{5}(:,4),'.b')
end
for r = 1:5
    if ~isempty(inputT{r})
        plot(timeT{r}(stopT{r}), tempT{r}(stopT{r}), 'ro');
    end
end
xlabel('date')
ylabel('Temperature ( \circC)')
axis([maxx-onehour maxx+onehour ystart yend])

a = axis;
rangex = a(2)-a(1);
rangey = a(4)-a(3);
idspot = [a(2)-0.2*rangex a(3)-0.2*rangey];
text(idspot(1), idspot(2), datestr(now));

endticks = [maxx-onehour, maxx-onehour + rangex/4, maxx-onehour + rangex/2, maxx+onehour - rangex/4, maxx+onehour];
set(gca,'xtick', endticks);
set(gca,'xticklabel', datestr(endticks));
set(gca,'xgrid','on');

saveas(1, [Tframestart, '_', Tframeend, '_Tall.pdf'],'pdf');
saveas(2, [Tframestart, '_', Tframeend, '_Tends.pdf'], 'pdf');

clear ystart yend minx maxx a range* *spot startticks endticks T*line T*dot

%% CTD Figures

ynames = char('Depth (m)', 'Salinity', 'Temperature ( \circC)', 'Density');
n = 3; %figure numbers
m = 3; %ctdburst columns, also ynames counter at -2
for t = {ctddepth, ctdsal, ctdtemp, calcdens}
blort = t{1};

ymin = floor(min([min(blort{1}(ctdstart{1}:ctdstop{1})) min(blort{2}(ctdstart{2}:ctdstop{2})) min(blort{3}(ctdstart{3}:ctdstop{3}))]));
ymax = ceil(max([max(blort{1}(ctdstart{1}:ctdstop{1})) max(blort{2}(ctdstart{2}:ctdstop{2})) max(blort{3}(ctdstart{3}:ctdstop{3}))]));

%CTD all
figure(n)
if ~isempty(blort{1}) %BOT: blue dots, cyan line
    plot(ctdtimestamp{1}, blort{1}, 'c');
    hold on
    plot(ctdburst{1}(:,2), ctdburst{1}(:,m),'.b');
end
if ~isempty(blort{2}) %SURF: black dots, gray line
    surfline = plot(ctdtimestamp{2}, blort{2});
    set(surfline, 'color', [.8 .8 .8])
    hold on
    plot(ctdburst{2}(:,2), ctdburst{2}(:,m),'.k');
end
if ~isempty(blort{3}) %YSI: dark green dots, green line
    plot(ctdtimestamp{3}, blort{3}, 'g');
    hold on
    ysidots = plot(ctdtimestamp{3}(ctdstart{3}:ctdstop{3}), blort{3}(ctdstart{3}:ctdstop{3}), '.');
    set(ysidots, 'color', [0 .75 0]);
end
for r = 1:3
    if ~isempty(blort{r})
        rcirc = plot(ctdtimestamp{r}(ctdstart{r}), blort{r}(ctdstart{r}), 'ro');
        set(rcirc, 'linewidth', 2);
        clear rcirc
        rcirc = plot(ctdtimestamp{r}(ctdstop{r}), blort{r}(ctdstop{r}), 'ro');
        set(rcirc, 'linewidth', 2);
        clear rcirc
    end
end
orient landscape
ylabel(ynames(m-2,:))
xlabel('date')
title([Tframestart, '\_', Tframeend, ' good ranges (dates): YSI = ', datestr(ctdtimestamp{3}(ctdstart{3})), ' - ', datestr(ctdtimestamp{3}(ctdstop{3})), sprintf(',\n'), 'Surface CTD = ', datestr(ctdtimestamp{2}(ctdstart{2})), ' - ', datestr(ctdtimestamp{2}(ctdstop{2})), ', Bottom CTD = ', datestr(ctdtimestamp{1}(ctdstart{1})), ' - ', datestr(ctdtimestamp{1}(ctdstop{1}))]);
axis([minstart-onehour.*24 maxstop+onehour.*24 ymin ymax])

a = axis;
rangex = a(2)-a(1);
rangey = a(4)-a(3);
allticks = [minstart, minstart + rangex/4, minstart + rangex/2, maxstop - rangex/4, maxstop];
set(gca,'xtick', allticks);
set(gca,'xticklabel', datestr(allticks));
set(gca,'xgrid','on')

spot1 = [a(2)+0.01*rangex a(4)-0.2*rangey];
leg1 = text(spot1(1), spot1(2), 'CTDBOT');
set(leg1, 'color', [0 0 1]);
set(leg1, 'fontweight', 'bold');
spot2 = [a(2)+0.01*rangex a(4)-0.15*rangey];
leg2 = text(spot2(1), spot2(2), 'YSI');
set(leg2, 'Color', [0 .75 0]);
set(leg2, 'fontweight', 'bold');
spot3 = [a(2)+0.01*rangex a(4)-0.1*rangey];
leg3 = text(spot3(1), spot3(2), 'CTDSURF');
set(leg3, 'Color', [0 0 0]);
set(leg3, 'fontweight', 'bold');
idspot = [a(2)-0.2*rangex a(3)-0.1*rangey];
text(idspot(1), idspot(2), datestr(now));

clear a range* leg* *spot* r

%CTD ends
figure(n+1)
subplot(2,1,1)
if ~isempty(blort{1}) %BOT: blue dots, cyan line
    plot(ctdtimestamp{1}, blort{1}, 'c');
    hold on
    plot(ctdburst{1}(:,2), ctdburst{1}(:,m),'.b');
end
if ~isempty(blort{2}) %SURF: black dots, gray line
    surfline = plot(ctdtimestamp{2}, blort{2});
    set(surfline, 'color', [.8 .8 .8])
    hold on
    plot(ctdburst{2}(:,2), ctdburst{2}(:,m),'.k');
end
if ~isempty(blort{3}) %YSI: dark green dots, green line
    plot(ctdtimestamp{3}, blort{3}, 'g');
    hold on
    ysidots = plot(ctdtimestamp{3}(ctdstart{3}:ctdstop{3}), blort{3}(ctdstart{3}:ctdstop{3}), '.');
    set(ysidots, 'color', [0 .75 0]);
end
for r = 1:3
    if ~isempty(blort{r})
        plot(ctdtimestamp{r}(ctdstart{r}), blort{r}(ctdstart{r}), 'ro');
    end
end
orient landscape
title([Tframestart, '\_', Tframeend, ' good ranges (scan): YSI = ', num2str(ctdstart{3}), ' - ', num2str(ctdstop{3}), sprintf(',\n'), 'Surface CTD = ', num2str(ctdstart{2}), ' - ', num2str(ctdstop{2}), ', Bottom CTD = ', num2str(ctdstart{1}), ' - ', num2str(ctdstop{1})]);
ylabel(ynames(m-2,:))
axis([minstart-onehour minstart+onehour ymin ymax])

a=axis;
rangex = a(2)-a(1);
rangey = a(4)-a(3);

startticks = [minstart-onehour, minstart-onehour + rangex/4, minstart-onehour + rangex/2, minstart+onehour - rangex/4, minstart+onehour];
set(gca,'xtick', startticks);
set(gca,'xticklabel', datestr(startticks));
set(gca,'xgrid', 'on');

spot1 = [a(2)+0.01*rangex a(4)-0.3*rangey];
leg1 = text(spot1(1), spot1(2), 'CTDBOT');
set(leg1, 'color', [0 0 1]);
set(leg1, 'fontweight', 'bold');
spot2 = [a(2)+0.01*rangex a(4)-0.2*rangey];
leg2 = text(spot2(1), spot2(2), 'YSI');
set(leg2, 'Color', [0 .75 0]);
set(leg2, 'fontweight', 'bold');
spot3 = [a(2)+0.01*rangex a(4)-0.1*rangey];
leg3 = text(spot3(1), spot3(2), 'CTDSURF');
set(leg3, 'Color', [0 0 0]);
set(leg3, 'fontweight', 'bold');

clear a range* *spot* leg*

subplot(2,1,2)

if ~isempty(blort{1}) %BOT: blue dots, cyan line
    plot(ctdtimestamp{1}, blort{1}, 'c');
    hold on
    plot(ctdburst{1}(:,2), ctdburst{1}(:,m),'.b');
end
if ~isempty(blort{2}) %SURF: black dots, gray line
    surfline = plot(ctdtimestamp{2}, blort{2});
    set(surfline, 'color', [.8 .8 .8])
    hold on
    plot(ctdburst{2}(:,2), ctdburst{2}(:,m),'.k');
end
if ~isempty(blort{3}) %YSI: dark green dots, green line
    plot(ctdtimestamp{3}, blort{3}, 'g');
    hold on
    ysidots = plot(ctdtimestamp{3}(ctdstart{3}:ctdstop{3}), blort{3}(ctdstart{3}:ctdstop{3}), '.');
    set(ysidots, 'color', [0 .75 0]);
end
for r = 1:3
    if ~isempty(blort{r})
        plot(ctdtimestamp{r}(ctdstop{r}), blort{r}(ctdstop{r}), 'ro');
    end
end
orient landscape
ylabel(ynames(m-2,:))
xlabel('date')
axis([maxstop-onehour maxstop+onehour ymin ymax])

a = axis;
rangex = a(2)-a(1);
rangey = a(4)-a(3);
endticks = [maxstop-onehour, maxstop-onehour + rangex/4, maxstop-onehour + rangex/2, maxstop+onehour - rangex/4, maxstop+onehour];
set(gca,'xtick', endticks);
set(gca,'xticklabel', datestr(endticks));
set(gca,'xgrid', 'on');

idspot = [a(2)-0.2*rangex a(3)-0.2*rangey];
text(idspot(1), idspot(2), datestr(now));

clear allticks startticks endticks r range* a idspot ymin ymax blort surfline ysidots rcirc

saveas(n, [Tframestart, '_', Tframeend, '_', ynames(m-2,1:3),'_all.pdf'],'pdf');
saveas(n+1, [Tframestart, '_', Tframeend, '_', ynames(m-2,1:3),'_ends.pdf'],'pdf');

n = n+2;
m = m+1;
end %end blort

clear l m n t minstart maxstop

%% all temperatures in one datafile
minest = min([min(burstTstamp{1}) min(burstTstamp{2}) min(burstTstamp{3}) min(burstTstamp{4}) min(burstTstamp{5}) ctdtimestamp{3}(ctdstart{3}) min(ctdburststamp{1}) min(ctdburststamp{2})]);
maxest = max([max(burstTstamp{1}) max(burstTstamp{2}) max(burstTstamp{3}) max(burstTstamp{4}) max(burstTstamp{5}) ctdtimestamp{3}(ctdstop{3}) max(ctdburststamp{1}) max(ctdburststamp{2})]);

for o = 1:floor((maxest - minest)/datenum(0,0,0,0,20,0))
    stampinterval(o) = minest + o*datenum(0,0,0,0,20,0);
end
stampdex = [minest stampinterval maxest];

temptable(1:length(stampdex), 1:10) = NaN;
temptable(:,1) = datenum(datestr(stampdex));

if ~isempty(ctdburst{2})
g = find(datenum(datestr(ctdburststamp{2}(1))) > (datenum(datestr(stampdex))-datenum(0,0,0,0,0,5)) & datenum(datestr(ctdburststamp{2}(1))) < (datenum(datestr(stampdex))+datenum(0,0,0,0,0,5)));
temptable(g:((g-1)+length(ctdburst{2})), 2) = ctdburst{2}(:, 5);
end
clear g
if ~isempty(burstT{1})
g = find(datenum(datestr(burstTstamp{1}(1))) > (datenum(datestr(stampdex))-datenum(0,0,0,0,0,5)) & datenum(datestr(burstTstamp{1}(1))) < (datenum(datestr(stampdex))+datenum(0,0,0,0,0,5)));
temptable(g:((g-1)+length(burstT{1})), 3) = burstT{1}(:,4);
end
clear g
if ~isempty(burstT{2})
g = find(datenum(datestr(burstTstamp{2}(1))) > (datenum(datestr(stampdex))-datenum(0,0,0,0,0,5)) & datenum(datestr(burstTstamp{2}(1))) < (datenum(datestr(stampdex))+datenum(0,0,0,0,0,5)));
temptable(g:((g-1)+length(burstT{2})), 4) = burstT{2}(:,4);
end
clear g
if ~isempty(burstT{3})
g = find(datenum(datestr(burstTstamp{3}(1))) > (datenum(datestr(stampdex))-datenum(0,0,0,0,0,5)) & datenum(datestr(burstTstamp{3}(1))) < (datenum(datestr(stampdex))+datenum(0,0,0,0,0,5)));
temptable(g:((g-1)+length(burstT{3})), 5) = burstT{3}(:,4);
end
clear g
if ~isempty(burstT{4})
g = find(datenum(datestr(burstTstamp{4}(1))) > (datenum(datestr(stampdex))-datenum(0,0,0,0,0,5)) & datenum(datestr(burstTstamp{4}(1))) < (datenum(datestr(stampdex))+datenum(0,0,0,0,0,5)));
temptable(g:((g-1)+length(burstT{4})), 6) = burstT{4}(:,4);
end
clear g
if ~isempty(ctdtemp{3})
g = find(datenum(datestr(ctdtimestamp{3}(ctdstart{3}))) == datenum(datestr(stampdex)));
temptable(g:((g-1)+length(ctdtemp{3}(ctdstart{3}:ctdstop{3}))), 7) = ctdtemp{3}(ctdstart{3}:ctdstop{3});
end
clear g
if ~isempty(burstT{5})
g = find(datenum(datestr(burstTstamp{5}(1))) > (datenum(datestr(stampdex))-datenum(0,0,0,0,0,5)) & datenum(datestr(burstTstamp{5}(1))) < (datenum(datestr(stampdex))+datenum(0,0,0,0,0,5)));
temptable(g:((g-1)+length(burstT{5})), 8) = burstT{5}(:,4);
end
clear g
if ~isempty(ctdburst{1})
g = find(datenum(datestr(ctdburststamp{1}(1))) > (datenum(datestr(stampdex))-datenum(0,0,0,0,0,5)) & datenum(datestr(ctdburststamp{1}(1))) < (datenum(datestr(stampdex))+datenum(0,0,0,0,0,5)));
temptable(g:(length(ctdburst{1}) + (g-1)), 9) = ctdburst{1}(:, 5);
end

%calculate jday
[yr, ~, ~, ~, ~, ~] = datevec(stampdex);
for i = 1:length(stampdex)
    excess = datenum(yr(i), 0, 1, 0, 0, 0); % Jan 1 is Day 0
    temptable(i, 10) = stampdex(i) - excess;
    clear excess
end
clear g o i yr excess

fidalltemp = fopen([Tframestart, '_', Tframeend, '_alltemp.dat'], 'w');

qlength = length(num2str(length(stampdex))) + 2;

for hpad = 1:qlength - 1
    fprintf(fidalltemp, ' ');
end
fprintf(fidalltemp, 'n       y  m  d  h mn  s    jday           sur       T01       T02       T03       T04       3ma       T05       bot\n');
for hpad = 1:qlength
    fprintf(fidalltemp, ' ');
end
fprintf(fidalltemp,'                                dep(m)    4.00      6.00      9.00      11.00     13.00     17.00     18.00     20.00\n');
clear hpad

for t = 1:length(stampdex)
    qpad = qlength - length(num2str(t));
    for pad = 1:qpad
        fprintf(fidalltemp, ' ');
    end
    clear *pad
    fprintf(fidalltemp, '%i', t);
    fprintf(fidalltemp, '    %s', datestr(temptable(t,1), 'yyyy mm dd HH MM SS'));
    fprintf(fidalltemp, '    %10.6f', temptable(t, 10));
    fprintf(fidalltemp, '    %6.2f', temptable(t, 2));
    fprintf(fidalltemp, '    %6.2f', temptable(t, 3));
    fprintf(fidalltemp, '    %6.2f', temptable(t, 4));
    fprintf(fidalltemp, '    %6.2f', temptable(t, 5));
    fprintf(fidalltemp, '    %6.2f', temptable(t, 6));
    fprintf(fidalltemp, '    %6.2f', temptable(t, 7));
    fprintf(fidalltemp, '    %6.2f', temptable(t, 8));
    fprintf(fidalltemp, '    %6.2f\n', temptable(t, 9));    
end
fclose(fidalltemp);

clear ans fidalltemp onehour stampinterval t minest mintime maxest maxtime stampdex qlength qpad
save([Tframestart, '_', Tframeend]);

%% Mooring_field_FOCAL.doc

doclog = fopen([Tframestart, '_', Tframeend, '_docsummary.txt'], 'w');

if ~isempty(inputctd{2})
    fprintf(doclog, 'Surface CTD\n\t: Use %i (%s) - %i (%s)\n\t: Removed are ', ctdstart{2}, datestr(ctdtimestamp{2}(ctdstart{2})), ctdstop{2}, datestr(ctdtimestamp{2}(ctdstop{2})));
    if ctdstart{2} - 1 == 0
        fprintf(doclog, 'none and ');
    elseif ctdstart{2} - 1 == 1
        fprintf(doclog, '1 (%s) and ', datestr(ctdtimestamp{2}(1)));
    else
        fprintf(doclog, '1 - %i (%s - %s) and ', ctdstart{2}-1, datestr(ctdtimestamp{2}(1)), datestr(ctdtimestamp{2}(ctdstart{2}-1)));
    end
    if ctdstop{2} == ctdscan{2}(end)
        fprintf(doclog, 'none\n');
    elseif ctdstop{2} +1 == ctdscan{2}(end)
        fprintf(doclog, '%i (%s)\n', ctdstop{2}+1, datestr(ctdtimestamp{2}(end)));
    else
        fprintf(doclog, '%i - end (%s - %s)\n', ctdstop{2}+1, datestr(ctdtimestamp{2}(ctdstop{2}+1)), datestr(ctdtimestamp{2}(end)));
    end
else
    fprintf(doclog, 'Surface CTD\n\t: no data\n');
end

if ~isempty(inputctd{3})
    fprintf(doclog, 'YSI\n\t: Use %i (%s) - %i (%s)\n\t: Removed are ', ctdstart{3}, datestr(ctdtimestamp{3}(ctdstart{3})), ctdstop{3}, datestr(ctdtimestamp{3}(ctdstop{3})));
    if ctdstart{3} - 1 == 0
        fprintf(doclog, 'none and ');
    elseif ctdstart{3} - 1 == 1
        fprintf(doclog, '1 (%s) and ', datestr(ctdtimestamp{3}(1)));
    else
        fprintf(doclog, '1 - %i (%s - %s) and ', ctdstart{3}-1, datestr(ctdtimestamp{3}(1)), datestr(ctdtimestamp{3}(ctdstart{3}-1)));
    end
    if ctdstop{3} == ctdscan{3}(end)
        fprintf(doclog, 'none\n');
    elseif ctdstop{3} +1 == ctdscan{3}(end)
        fprintf(doclog, '%i (%s)\n', ctdstop{3}+1, datestr(ctdtimestamp{3}(end)));
    else
        fprintf(doclog, '%i - end (%s - %s)\n', ctdstop{3}+1, datestr(ctdtimestamp{3}(ctdstop{3}+1)), datestr(ctdtimestamp{3}(end)));
    end
else
    fprintf(doclog, 'YSI\n\t: no data\n');
end

if ~isempty(inputctd{1})
    fprintf(doclog, 'Bottom CTD\n\t: Use %i (%s) - %i (%s)\n\t: Removed are ', ctdstart{1}, datestr(ctdtimestamp{1}(ctdstart{1})), ctdstop{1}, datestr(ctdtimestamp{1}(ctdstop{1})));
    if ctdstart{1} - 1 == 0
        fprintf(doclog, 'none and ');
    elseif ctdstart{1} - 1 == 1
        fprintf(doclog, '1 (%s) and ', datestr(ctdtimestamp{1}(1)));
    else
        fprintf(doclog, '1 - %i (%s - %s) and ', ctdstart{1}-1, datestr(ctdtimestamp{1}(1)), datestr(ctdtimestamp{1}(ctdstart{1}-1)));
    end
    if ctdstop{1} == ctdscan{1}(end)
        fprintf(doclog, 'none\n');
    elseif ctdstop{1} +1 == ctdscan{1}(end)
        fprintf(doclog, '%i (%s)\n', ctdstop{1}+1, datestr(ctdtimestamp{1}(end)));
    else
        fprintf(doclog, '%i - end (%s - %s)\n', ctdstop{1}+1, datestr(ctdtimestamp{1}(ctdstop{1}+1)), datestr(ctdtimestamp{1}(end)));
    end
else
    fprintf(doclog, 'Bottom CTD\n\t: no data\n');
end

for w = 1:5
    if ~isempty(inputT{w})
        fprintf(doclog, '%s\n\t: Use %i (%s) - %i (%s)\n\t: Removed are ', char(stationname(inputT{w})), startT{w}, datestr(timeT{w}(startT{w})), stopT{w}, datestr(timeT{w}(stopT{w})));
        if startT{w} - 1 == 0
            fprintf(doclog, 'none and ');
        elseif startT{w} - 1 == 1
            fprintf(doclog, '1 (%s) and ', datestr(timeT{w}(1)));
        else
            fprintf(doclog, '1 - %i (%s - %s) and ', startT{w}-1, datestr(timeT{w}(1)), datestr(timeT{w}(startT{w}-1)));
        end
        
        if stopT{w} == length(timeT{w})
            fprintf(doclog, 'none\n');
        elseif stopT{w} +1 == length(timeT{w})
            fprintf(doclog, '%i (%s)\n', stopT{w}+1, datestr(timeT{w}(end)));
        else
            fprintf(doclog, '%i - end (%s - %s)\n', stopT{w}+1, datestr(timeT{w}(stopT{w}+1)), datestr(timeT{w}(end)));
        end
    else
        fprintf(doclog, '%s\n\t: no data\n', char(stationname(inputT{w})));
    end
end

fprintf(doclog, '\nReplaced values:\n');
if ~isempty(interpnums)
    for J = 1:length(interpnums)
        switch interpnums(J,1)
            case 1
                fprintf(doclog, 'Bottom CTD\t');
            case 2
                fprintf(doclog, 'Surface CTD\t');
            case 3
                fprintf(doclog, 'YSI\t\t');
        end
        switch interpnums(J, 2)
            case 1
                fprintf(doclog, 'salinity\t\t');
            case 2
                fprintf(doclog, 'temperature\t\t');
            case 3
                fprintf(doclog, 'density\t\t');
            case 4
                fprintf(doclog, 'pressure\t\t');
            case 5
                fprintf(doclog, 'conductivity\t\t');
            case 6
                fprintf(doclog, 'specific conductance\t\t');
            case 7
                fprintf(doclog, 'sound velocity\t\t');
            case 8
                fprintf(doclog, 'depth\t\t');
            case 9
                fprintf(doclog, 'Julian Day\t\t');
        end
        fprintf(doclog, '%i\t%6.4f\n', interpnums(J,3), interpnums(J,4));
    end
else
    fprintf(doclog, 'none\n');
end

fclose(doclog);
clear J w

%% CTD cast

if exist('inputcast', 'var') == 1

station=char(stations(inputcast));
humandate=startdates(inputcast);
time=starttimes(inputcast);
inputmin = scanstarts(inputcast);
inputmax = scanstops(inputcast);
compdate=datestr(datevec(humandate), 'yyyymmdd');

datablock=importdata([station, '-matlab.txt']);
scan=datablock(:,1);
depth=datablock(:,3);
sal=datablock(:,4);
temp=datablock(:,5);
dens=datablock(:,6);
press=datablock(:,7);
cond=datablock(:,8);
spcond=datablock(:,9);
oxy=datablock(:,10);
oxysat=datablock(:,11);
ph=datablock(:,12);
chla=datablock(:,13);
%par=datablock(:,14);
par(1:length(scan)) = NaN;
cdom=datablock(:,14);
%batt=datablock(:,15);
batt(1:length(scan)) = NaN;
%xmiss=datablock(:,16);
xmiss(1:length(scan)) = NaN;

%xmiss(1:length(scan)) = NaN; %<-- use this if the variable does not exist

bin = .5;

if inputmax > 0
    firstmax = inputmax;
else
    firstmax = find(depth==max(depth),1);
end

% find the minimum of the last numbers less than 1 before the first max
if inputmin < 10
    mincap = find(depth(1:firstmax) < inputmin,1,'last');
    k = 0;
    while depth(mincap - k) - depth(mincap - (k+1)) > 0
        k = k+1;
    end
    lastmin = mincap - k;
else
    lastmin = inputmin;
end

%% depth check plot
fig1 = figure(11);
orient landscape
subplot(2,1,1)
[Hyy,Hdepth,Hsal] = plotyy(scan,depth,scan,sal);
title([Tframestart, '\_', Tframeend, ' CTD cast: ', char(humandate), ' ', char(time)]); 
set(get(Hyy(1),'Ylabel'),'String','Depth (m)')
set(Hyy(1),'Ycolor',[0 0 0])
set(Hyy(1),'Ylim',[0 40])
set(Hyy(1),'ytick',[0 10 20 30 40])
set(get(Hyy(2),'Ylabel'),'String','Salinity (psu)')
set(Hyy(2),'YColor', [0 0 0])
set(Hyy(2),'Ylim',[0 40])
set(Hyy(2),'ytick',[0 10 20 30 40])
set(Hdepth, 'color','k')
set(Hdepth,'marker','.')
set(Hdepth,'markersize',4)
set(Hsal, 'color','b')
set(Hsal,'marker','.')
set(Hsal,'markersize',4)
hold on
plot(scan(lastmin),depth(lastmin),'ro')
plot(scan(lastmin),sal(lastmin),'ro')
plot(scan(firstmax),depth(firstmax),'ro')
plot(scan(firstmax),sal(firstmax),'ro')

clear Hyy

subplot(2,1,2)
[Hyy, Htemp, Hdens] = plotyy(scan,temp,scan,dens);
xlabel('count number')
set(get(Hyy(1),'Ylabel'),'String','Temperature ( \circC)')
set(Hyy(1),'Ylim',[0 40])
set(Hyy(1),'YColor', [0 0 0])
set(Hyy(1),'ytick', [0 10 20 30 40])
set(get(Hyy(2),'Ylabel'),'String','Sigma-t (kg/m^3)')
set(Hyy(2),'Ylim',[0 40])
set(Hyy(2),'YColor',[0 0 0])
set(Hyy(2),'ytick',[0 10 20 30 40])
set(Htemp, 'color','r')
set(Htemp,'marker','.')
set(Htemp,'markersize',4)
set(Hdens, 'color','g')
set(Hdens,'marker','.')
set(Hdens,'markersize',4)
hold on
plot(scan(lastmin),temp(lastmin),'ro')
plot(scan(lastmin),dens(lastmin),'ro')
plot(scan(firstmax),temp(firstmax),'ro')
plot(scan(firstmax),dens(firstmax),'ro')

a = axis;
rangex = a(2)-a(1);
rangey = a(4)-a(3);
spot = [a(2)-0.25*rangex a(4)-0.05*rangey];
text(spot(1),spot(2), ['good data range: ', num2str(lastmin), ' to ', num2str(firstmax)]);
idspot = [a(2)-0.2*rangex a(3)-0.15*rangey];
text(idspot(1), idspot(2), datestr(now));
clear a range* *spot

%% bin averaging

% make bins
numbins = floor(max(depth)/bin);
binpoint = bin*(1:numbins);
for i = 1:length(binpoint)
    j = find(depth(lastmin:firstmax) >= binpoint(i)-.5*bin & depth(lastmin:firstmax) < binpoint(i)+.5*bin) + (lastmin-1);
    avgdepth(i) = mean(depth(j));
    avgsal(i) = mean(sal(j));
    avgtemp(i) = mean(temp(j));
    avgdens(i) = mean(dens(j));
    avgpress(i) = mean(press(j));
    avgcond(i) = mean(cond(j)*.0001);
    avgspcond(i) = mean(spcond(j)*.0001);
    avgoxy(i) = mean(oxy(j));
    avgoxysat(i) = mean(oxysat(j));
    avgph(i) = mean(ph(j));
    avgchla(i) = mean(chla(j));
    avgpar(i) = mean(par(j));
    avgcdom(i) = mean(cdom(j));
    avgbatt(i) = mean(batt(j));
    avgxmiss(i) = mean(xmiss(j));
    nbin(i) = length(j);
end

%% plots of binaveraged on non-binaveraged
fig2 = figure(12);
orient landscape
subplot(1,4,1)
plot(sal(lastmin:firstmax), depth(lastmin:firstmax),'.-k')
hold on 
plot(avgsal, binpoint,'ro')
xlabel('Salinity (psu)')
ylabel('Depth (m)')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')

subplot(1,4,2)
plot(temp(lastmin:firstmax), depth(lastmin:firstmax),'.-k')
hold on
plot(avgtemp, binpoint,'ro')
xlabel('Temperature ( \circC)')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')
set(gca,'yticklabel',{''})

subplot(1,4,3)
plot(dens(lastmin:firstmax), depth(lastmin:firstmax),'.-k')
hold on
plot(avgdens, binpoint, 'ro')
xlabel('Density (kg/m^3)')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')
set(gca,'yticklabel',{''})

subplot(1,4,4)
plot(press(lastmin:firstmax), depth(lastmin:firstmax),'.-k')
hold on
plot(avgpress, binpoint, 'ro')
xlabel('Pressure (db)')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')
set(gca,'yticklabel',{''})

a = axis;
rangex = a(2)-a(1);
rangey = a(4)-a(3);
spot = [a(2)-0.75*rangex a(3)+0.05*rangey];
text(spot(1),spot(2),[Tframestart, '\_', Tframeend, sprintf('\n'), char(humandate), ' ', char(time)]);
idspot = [a(1)+0.05*rangex a(4)+0.1*rangey];
text(idspot(1), idspot(2), datestr(now));
clear a range* *spot

fig3 = figure(13);
orient landscape
subplot(1,4,1)
plot(cond(lastmin:firstmax)*.0001, depth(lastmin:firstmax),'.-k')
hold on
plot(avgcond, binpoint, 'ro')
xlabel('Conductivity (10\circ4 uS/cm)')
ylabel('Depth (m)')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')

subplot(1,4,2)
plot(spcond(lastmin:firstmax)*.0001, depth(lastmin:firstmax),'.-k')
hold on
plot(avgspcond, binpoint,'ro')
xlabel('S Conductivity (10\circ4 uS/cm)')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')
set(gca,'yticklabel',{''})

subplot(1,4,3)
plot(oxy(lastmin:firstmax), depth(lastmin:firstmax),'.-k')
hold on
plot(avgoxy, binpoint, 'ro')
xlabel('Dissolved Oxygen (mg/l)')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')
set(gca,'yticklabel',{''})

subplot(1,4,4)
plot(oxysat(lastmin:firstmax), depth(lastmin:firstmax),'.-k')
hold on
plot(avgoxysat, binpoint,'ro')
xlabel('Oxygen Saturation (mg/l)')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')
set(gca,'yticklabel',{''})

a = axis;
rangex = a(2)-a(1);
rangey = a(4)-a(3);
idspot = [a(1)+0.05*rangex a(4)+0.1*rangey];
text(idspot(1), idspot(2), datestr(now));
clear a range* *spot

fig4 = figure(14);
orient landscape
subplot(1,4,1)
plot(ph(lastmin:firstmax), depth(lastmin:firstmax),'.-k')
hold on
plot(avgph, binpoint,'ro')
xlabel('pH')
ylabel('Depth (m)')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')

subplot(1,4,2)
plot(chla(lastmin:firstmax), depth(lastmin:firstmax),'.-k')
hold on
plot(avgchla, binpoint, 'ro')
xlabel('Fluorescence (mg/m^3)')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')
set(gca,'yticklabel',{''})

subplot(1,4,3)
plot(par(lastmin:firstmax), depth(lastmin:firstmax),'.-k')
%plot(0, depth(lastmin:firstmax),'.k') % <-- the kludge for if the variable doesn't exist
hold on
plot(avgpar, binpoint, 'ro')
xlabel('PAR')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')
set(gca,'yticklabel',{''})

subplot(1,4,4)
plot(cdom(lastmin:firstmax), depth(lastmin:firstmax),'.-k')
hold on
plot(avgcdom, binpoint,'ro')
xlabel('CDOM (mg/m^3)')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')
set(gca,'yticklabel',{''})

a = axis;
rangex = a(2)-a(1);
rangey = a(4)-a(3);
idspot = [a(1)+0.05*rangex a(4)+0.1*rangey];
text(idspot(1), idspot(2), datestr(now));
clear a range* *spot

fig5 = figure(15);
orient landscape
subplot(1,4,1)
plot(batt(lastmin:firstmax), depth(lastmin:firstmax),'.-k')
hold on
plot(avgbatt, binpoint, 'ro')
xlabel('Beam attenuation')
ylabel('Depth (m)')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')

subplot(1,4,3)
plot(xmiss(lastmin:firstmax), depth(lastmin:firstmax),'.-k')
hold on
plot(avgxmiss, binpoint, 'ro')
xlabel('Beam transmission')
set(gca,'YDir','reverse')
set(gca,'xaxislocation','top')
set(gca,'yticklabel',{''})

a = axis;
rangex = a(2)-a(1);
rangey = a(4)-a(3);
idspot = [a(1)+0.05*rangex a(4)+0.1*rangey];
text(idspot(1), idspot(2), datestr(now));
clear a range* *spot

%% output to text files and pdfs
fid = fopen([station, '_est.dat'],'w');
fprintf(fid, '%s   %s   %s %s: estimated data at %2.2f m intervals\n', compdate, station, char(humandate), char(time), bin);
fprintf(fid, '    i,    d(m),    S(psu),    T(C),    sigma-t,    p(db),    cd,    sc(10**4 uS/cm),    DO(mg/l), DOS(.), pH,    FL(mg/m3), PAR,      CDOM(mg/m3), b.att,    b.trans,  N,    meandepth(m)\n');

fid3 = fopen([station, '_binned.csv'],'w');
fprintf(fid3, 'bindepth, sal, temp, dens, press, cond, spcond, DO, DOsat, pH, fluor, PAR, CDOM, batt, xmiss, N bin, avgdepth\n');

qlength = length(num2str(length(avgdepth))) + 2;
bplength = length(num2str(floor(max(binpoint)))) + 4;
sallength = length(num2str(floor(max(avgsal)))) + 4;
templength = length(num2str(floor(max(avgtemp)))) + 4;
denslength = length(num2str(floor(max(avgdens)))) + 4;
presslength = length(num2str(floor(max(avgpress)))) + 4;
condlength = length(num2str(floor(max(avgcond)))) + 4;
spcondlength = length(num2str(floor(max(avgspcond)))) + 4;
oxylength = length(num2str(floor(max(avgoxy)))) + 6;
oxysatlength = length(num2str(floor(max(avgoxysat)))) + 4;
phlength = length(num2str(floor(max(avgph)))) + 4;
chlalength = length(num2str(floor(max(avgchla)))) + 4;
parlength = length(num2str(floor(max(avgpar)))) + 4;
cdomlength = length(num2str(floor(max(avgcdom)))) + 4;
battlength = length(num2str(floor(max(avgbatt)))) + 4;
xmisslength = length(num2str(floor(max(avgxmiss)))) + 4;
nbinlength = length(num2str(floor(max(nbin)))) + 4;
depthlength = length(num2str(floor(max(avgdepth)))) + 4;

for q=1:length(avgdepth)
    qpad = qlength - length(num2str(q));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%i', q);

    qpad = bplength - length(num2str(floor(binpoint(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid,'%2.2f', binpoint(q));
    fprintf(fid3, '%2.2f,', binpoint(q));

    qpad = sallength - length(num2str(floor(avgsal(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgsal(q));
    fprintf(fid3, '%6.4f,', avgsal(q));

    qpad = templength - length(num2str(floor(avgtemp(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgtemp(q));
    fprintf(fid3, '%6.4f,', avgtemp(q));
    
    qpad = denslength - length(num2str(floor(avgdens(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgdens(q));
    fprintf(fid3, '%6.4f,', avgdens(q));

    qpad = presslength - length(num2str(floor(avgpress(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgpress(q));
    fprintf(fid3, '%6.4f,', avgpress(q));

    qpad = condlength - length(num2str(floor(avgcond(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgcond(q));
    fprintf(fid3, '%6.4f,', avgcond(q));

    qpad = spcondlength - length(num2str(floor(avgspcond(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgspcond(q));
    fprintf(fid3, '%6.4f,', avgspcond(q));

    qpad = oxylength - length(num2str(floor(avgoxy(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgoxy(q));
    fprintf(fid3, '%6.4f,', avgoxy(q));

    qpad = oxysatlength - length(num2str(floor(avgoxysat(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgoxysat(q));
    fprintf(fid3, '%6.4f,', avgoxysat(q));

    qpad = phlength - length(num2str(floor(avgph(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgph(q));
    fprintf(fid3, '%6.4f,', avgph(q));

    qpad = chlalength - length(num2str(floor(avgchla(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgchla(q));
    fprintf(fid3, '%6.4f,', avgchla(q));

     qpad = parlength - length(num2str(floor(avgpar(q))));
     for pad = 1:qpad
         fprintf(fid, ' ');
     end
     clear *pad
    fprintf(fid, '%6.4f', avgpar(q));
    fprintf(fid3, '%6.4f,', avgpar(q));
%fprintf(fid, '    NaN');
%fprintf(fid3, 'NaN,');

    qpad = cdomlength - length(num2str(floor(avgcdom(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgcdom(q));
    fprintf(fid3, '%6.4f,', avgcdom(q));

    qpad = battlength - length(num2str(floor(avgbatt(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgbatt(q));
    fprintf(fid3, '%6.4f,', avgbatt(q));
    
    qpad = xmisslength - length(num2str(floor(avgxmiss(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f', avgxmiss(q));
    fprintf(fid3, '%6.4f,', avgxmiss(q));

    qpad = nbinlength - length(num2str(nbin(q)));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%i', nbin(q));
    fprintf(fid3, '%i,', nbin(q));

    qpad = depthlength - length(num2str(floor(avgdepth(q))));
    for pad = 1:qpad
        fprintf(fid, ' ');
    end
    clear *pad
    fprintf(fid, '%6.4f\n', avgdepth(q));
    fprintf(fid3, '%6.4f\n', avgdepth(q));
end
fclose(fid);
fclose(fid3);
clear *length *pad

fid2 = fopen([station, '_goodrange.csv'], 'w');
fprintf(fid2, 'depth, sal, temp, dens, press, cond, spcond, DO, DOsat, pH, fluor, PAR, CDOM, batt, xmiss, scan\n');
for t = 0:(length(lastmin:firstmax) - 1)
    fprintf(fid2, '%6.4f,', depth(lastmin + t));
    fprintf(fid2, '%6.4f,', sal(lastmin + t));
    fprintf(fid2, '%6.4f,', temp(lastmin + t));
    fprintf(fid2, '%6.4f,', dens(lastmin + t));
    fprintf(fid2, '%6.4f,', press(lastmin + t));
    fprintf(fid2, '%6.4f,', cond(lastmin + t));
    fprintf(fid2, '%6.4f,', spcond(lastmin + t));
    fprintf(fid2, '%6.4f,', oxy(lastmin + t));
    fprintf(fid2, '%6.4f,', oxysat(lastmin+t));
    fprintf(fid2, '%6.4f,', ph(lastmin + t));
    fprintf(fid2, '%6.4f,', chla(lastmin + t));
    fprintf(fid2, '%6.4f,', par(lastmin + t));
%fprintf(fid2, 'NaN,');  % <-- kludge for missing variables
    fprintf(fid2, '%6.4f,', cdom(lastmin + t));
    fprintf(fid2, '%6.4f,', batt(lastmin + t));
    fprintf(fid2, '%6.4f,', xmiss(lastmin + t));
    fprintf(fid2, '%i\n', scan(lastmin + t));
end
fclose(fid2);

saveas(fig1, [station,'_checkdep.pdf'], 'pdf');
saveas(fig2, [station, '_est1.pdf'], 'pdf');
saveas(fig3, [station, '_est2.pdf'], 'pdf');
saveas(fig4, [station, '_est3.pdf'], 'pdf');
saveas(fig5, [station, '_est4.pdf'], 'pdf');

fidlog = fopen([station, '.log'], 'w');
fprintf(fidlog, 'Diagnostic output\n\nCruise ID = %s\n\nfile: %s\n\n', station, [station,'.cnv']);
sum2 = [num2str(scan(firstmax) - scan(lastmin)), '     ', num2str(lastmin), '     ', num2str(depth(lastmin)), '     ', num2str(firstmax), '     ', num2str(depth(firstmax)) ];
fprintf(fidlog,'nvalid nrd_min   dep_min  nrd_max   dep_max\n%s\n\nmin, avg, max\n', sum2);
fprintf(fidlog, 'T (C)               ');
fprintf(fidlog, '%6.2f        ', min(temp(lastmin:firstmax)));
fprintf(fidlog, '%6.2f        ', mean(temp(lastmin:firstmax)));
fprintf(fidlog, '%6.2f\n', max(temp(lastmin:firstmax)));
fprintf(fidlog, 'S (psu)             ');
fprintf(fidlog, '%6.2f        ', min(sal(lastmin:firstmax)));
fprintf(fidlog, '%6.2f        ', mean(sal(lastmin:firstmax)));
fprintf(fidlog, '%6.2f\n', max(sal(lastmin:firstmax)));
fprintf(fidlog, 'Sg_t (kg/m3)        ');
fprintf(fidlog, '%6.3f        ', min(dens(lastmin:firstmax)));
fprintf(fidlog, '%6.3f        ', mean(dens(lastmin:firstmax)));
fprintf(fidlog, '%6.3f\n', max(dens(lastmin:firstmax)));
fprintf(fidlog, 'Pressure (db)       ');
fprintf(fidlog, '%6.2f        ', min(press(lastmin:firstmax)));
fprintf(fidlog, '%6.2f        ', mean(press(lastmin:firstmax)));
fprintf(fidlog, '%6.2f\n', max(press(lastmin:firstmax)));
fprintf(fidlog, 'Cond (S/cm)         ');
fprintf(fidlog, '%6.4f    ', min(cond(lastmin:firstmax)));
fprintf(fidlog, '%6.4f    ', mean(cond(lastmin:firstmax)));
fprintf(fidlog, '%6.4f\n', max(cond(lastmin:firstmax)));
fprintf(fidlog, 'Sp. Cond (S/cm)     ');
fprintf(fidlog, '%6.4f    ', min(spcond(lastmin:firstmax)));
fprintf(fidlog, '%6.4f    ', mean(spcond(lastmin:firstmax)));
fprintf(fidlog, '%6.4f\n', max(spcond(lastmin:firstmax)));
fprintf(fidlog, 'DO (mg/l)           ');
fprintf(fidlog, '%6.2f        ', min(oxy(lastmin:firstmax)));
fprintf(fidlog, '%6.2f        ', mean(oxy(lastmin:firstmax)));
fprintf(fidlog, '%6.2f\n', max(oxy(lastmin:firstmax)));
fprintf(fidlog, 'DO_sat(%%)          ');
fprintf(fidlog, '%6.2f        ', min(oxysat(lastmin:firstmax)));
fprintf(fidlog, '%6.2f        ', mean(oxysat(lastmin:firstmax)));
fprintf(fidlog, '%6.2f\n', max(oxysat(lastmin:firstmax)));
fprintf(fidlog, 'pH                  ');
fprintf(fidlog, '%6.3f        ', min(ph(lastmin:firstmax)));
fprintf(fidlog, '%6.3f        ', mean(ph(lastmin:firstmax)));
fprintf(fidlog, '%6.3f\n', max(ph(lastmin:firstmax)));
fprintf(fidlog, 'Fluor (mg/m3)       ');
fprintf(fidlog, '%6.3f        ', min(chla(lastmin:firstmax)));
fprintf(fidlog, '%6.3f        ', mean(chla(lastmin:firstmax)));
fprintf(fidlog, '%6.3f\n', max(chla(lastmin:firstmax)));
fprintf(fidlog, 'CDOM (mg/m3)        ');
fprintf(fidlog, '%6.3f        ', min(cdom(lastmin:firstmax)));
fprintf(fidlog, '%6.3f        ', mean(cdom(lastmin:firstmax)));
fprintf(fidlog, '%6.3f\n', max(cdom(lastmin:firstmax)));
fprintf(fidlog, 'B. Atten.           ');
fprintf(fidlog, '%6.3f        ', min(batt(lastmin:firstmax)));
fprintf(fidlog, '%6.3f        ', mean(batt(lastmin:firstmax)));
fprintf(fidlog, '%6.3f\n', max(batt(lastmin:firstmax)));
fprintf(fidlog, 'B. Trans.           ');
fprintf(fidlog, '%6.3f        ', min(xmiss(lastmin:firstmax)));
fprintf(fidlog, '%6.3f        ', mean(xmiss(lastmin:firstmax)));
fprintf(fidlog, '%6.3f\n', max(xmiss(lastmin:firstmax)));

fprintf(fidlog, '\ntotdep for dz2 estimation = %s        %i\n', num2str(binpoint(end)), length(binpoint));
fclose(fidlog);

end %end if inputcast exists
