function [firstest, lastest] = moorburst(startscan, starttime, stopscan, stoptime)
%this function finds the beginning and ending points of doing the
%5-minute averages per 20 minutes, for mooring thermistors and CTDs that
%take readings every 1 minute.

% 1. from beginning of good range, find first next :13, :33, or :53 by
% assuming 1 minute per scan.
[~, ~, ~, ~, mn, ~] = datevec(starttime);

i = startscan;
if mn < 13
    i = i + (13 - mn);
elseif mn > 13 && mn < 33
    i = i + (33 - mn);
elseif mn > 33 && mn < 53
    i = i + (53 - mn);
elseif mn > 53
    i = i + (60 - mn) + 13;
end
firstest = i;

clear mn
% 2. from good-range end, find nearest previous :32, :52, or :12 by
% assuming 1 minute per scan.
[~, ~, ~, ~, mn, ~] = datevec(stoptime);

j = stopscan;
if mn < 12
    j = j - mn - 8;
elseif mn > 12 && mn < 32
    j = j - (mn - 12);
elseif mn > 32 && mn < 52
    j = j - (mn - 32);
elseif mn > 52
    j = j - (mn - 52);
end
lastest = j;
