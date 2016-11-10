%%%%
%%
%% [04, 11, 2016] Dirk Slawinski
%% - time to give Ken J's code a go. One issue right away is  that .cal 
%% file has different name-value pairs and also no "," between some entries
%% For example Ken's code looks for "CalTemp" and our file has "T_CAL" which
%% I am sure is the same. The AGRO processing document says a Satlantic float
%% will have "T_CAL_SWA"
%%
%% - OK so a quick plot shows there seems to be a shift to hight NO3 values
%% with this when compared with the discrete values ... ooo boy ... best to
%% properly read some more of the details on the scripts.
%%
%%%%

clear all

fltNum = 0391; %0390;
sunaId = 0417; %0416;
prflNum = 10; %1; %3;
infile =  sprintf('../f%04d/SNA%04dA.cal', fltNum, sunaId);
cal = parseNO3cal(infile);
infile = sprintf('../f%04d/%04d.%03d.isus', fltNum, fltNum, prflNum);
dat = parse_NO3msg(infile);

webDir = '/OSM/HOME-PER/sla083/htdocs/bioargo';
fltDat = load(sprintf('%s/%04dsbedat.mat', webDir, fltNum));
fltPrsdx = 5;
fltTmpdx = 4;
fltSltdx = 6;
fltNo3dx = 11;
prflDat = fltDat.dat(prflNum).profile;

% presure correction flag
%prsFlg = 0;
prsFlg = 1;
cal.pres_coef = 0.02;
NO3 = calc_FLOAT_NO3(dat, cal, prsFlg);

% now lets make some test plots
prsdx = 3;
tmpdx = 4;
sltdx = 5;
no3dx = 6;

ds_clear
orient landscape;

nr = 1;
nc = 3;
subplot(nr, nc, 1)
plot(prflDat(:, fltTmpdx), -prflDat(:, fltPrsdx), '.r');
hold on
plot(NO3(:, tmpdx), -NO3(:, prsdx), '.b');
title('Temp [C^o]');

subplot(nr, nc, 2)
plot(prflDat(:, fltSltdx), -prflDat(:, fltPrsdx), '.r');
hold on
plot(NO3(:, sltdx), -NO3(:, prsdx), '.b');
title('Salt [PSU]');

subplot(nr, nc, 3)
plot(prflDat(:, fltNo3dx), -prflDat(:, fltPrsdx), '.r');
hold on
nShft = 0;
% hmmm, looks like the MBARI calc is ~2.5 greater than the float cals
%nShft = 2.5;
% Nick H-M suggested it might be sqrt of 5 it was not
%nShft = sqrt(5);
% alright lets calc the mean shift
ndx = find(~isnan(prflDat(:, fltNo3dx)));
%nShft = mean(NO3(:, no3dx) - prflDat(ndx, fltNo3dx));

plot(NO3(:, no3dx)-nShft, -NO3(:, prsdx), '.b');
title('NO_3 [umol/l]');
if(nShft == 0)
   legend('.msg', '.isus', 'location', 'southwest')
else
   legend('.msg', sprintf('.isus -% 0.1f', nShft)', 'location', 'southwest')
end % nShft ~=0

oname = sprintf('../plots/f%04d-%03d.png', fltNum, prflNum);
ds_print(oname, ds_islandscape, 0);


