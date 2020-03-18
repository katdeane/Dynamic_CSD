function BrownScyth(homedir)

% folder:   D:\MyCode\Dynamic_CSD_Analysis\DATA\avrec_compare
% Input:    stats from single trial avrec and relres 
% Output:   BrownScyth statistical comparison of variances; box plots and
%           stats saved as pictures

warning('OFF');
dbstop if error

% Change directory to your working folder
if ~exist('homedir','var')
    if exist('D:\MyCode\Dynamic_CSD','dir') == 7
        cd('D:\MyCode\Dynamic_CSD');
    elseif exist('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD','dir') == 7
        cd('C:\Users\kedea\Documents\Work Stuff\Dynamic_CSD')
    end
    
    homedir = pwd;
    addpath(genpath(homedir));
end
cd(homedir);
cd Data; cd avrec_compare;

figlist = {'vMusc_Box','vMusc_Stat','vAw_Box','vAw_Stat'};

%% Load and sort the data 
load('AvrecPlotData_single.mat','plottune'); AVRECdata = plottune;
load('RelresPlotData_single.mat','plottune'); RELRESdata = plottune; 
load('AbsresPlotData_single.mat','plottune'); ABSRESdata = plottune; clear plotdata

% pull out the needed structure 
AVREC.peakamp = AVRECdata.peakamp';
AVREC.peaklate = AVRECdata.peaklate';
AVREC.tgroup = AVRECdata.tgroup';
AVREC.frqz = AVRECdata.frqz';
RELRES.peakamp = RELRESdata.peakamp';
RELRES.peaklate = RELRESdata.peaklate';
RELRES.tgroup = RELRESdata.tgroup';
RELRES.frqz = RELRESdata.frqz';
ABSRES.peakamp = ABSRESdata.peakamp';
ABSRES.peaklate = ABSRESdata.peaklate';
ABSRES.tgroup = ABSRESdata.tgroup';
ABSRES.frqz = ABSRESdata.frqz';

% convert to table, and find the groups to compare at specific frequencies
% AV_BF is avrec BF, AV_m2 is avrec BF-2, RE is relres, AB is absres
AVREC_T = struct2table(AVREC);
AV_BFT = AVREC_T(strcmp(AVREC_T.frqz,'d BF'),1:3);
[AV_BF,AV_BFname]=findgroups(AV_BFT.tgroup);
AV_m2T = AVREC_T(strcmp(AVREC_T.frqz,'b -2'),1:3);
[AV_m2,AV_m2name]=findgroups(AV_m2T.tgroup);

RELRES_T = struct2table(RELRES);
RE_BFT = RELRES_T(strcmp(RELRES_T.frqz,'d BF'),1:3);
[RE_BF,RE_BFname]=findgroups(RE_BFT.tgroup);
RE_m2T = RELRES_T(strcmp(RELRES_T.frqz,'b -2'),1:3);
[RE_m2,RE_m2name]=findgroups(RE_m2T.tgroup);

ABSRES_T = struct2table(ABSRES);
AB_BFT = ABSRES_T(strcmp(ABSRES_T.frqz,'d BF'),1:3);
[AB_BF,AB_BFname]=findgroups(AB_BFT.tgroup);
AB_m2T = ABSRES_T(strcmp(ABSRES_T.frqz,'b -2'),1:3);
[AB_m2,AB_m2name]=findgroups(AB_m2T.tgroup);


%% AVREC PEAKAMP

% BF
clear peakamp
for ii =1:4
    %cut out groups and stick them side by side in a struct
    peakamp.(AV_BFname{ii}) = AV_BFT.peakamp(AV_BF == ii);
end
%stick the struct in a table instead to spread it out properly
peakAmpT = struct2table(peakamp); 
%tables can't go into vartestn so we need a matrix, no table2mat so 2 steps:
peakAmpC = table2cell(peakAmpT); PEAKamp = cell2mat(peakAmpC); 
%2 comparisons, ket vs awake and ket vs mus
PEAKamp1 = PEAKamp(:,2:3); %ketamine and awake
PEAKamp2 = PEAKamp(:,2:2:4); %ketamine and muscimol

[AV_Amp_BF_P1,AV_Amp_BF_STATS1] = vartestn(PEAKamp1, 'testtype','BrownForsythe');
[AV_Amp_BF_P2,AV_Amp_BF_STATS2] = vartestn(PEAKamp2, 'testtype','BrownForsythe');

%get each figure, save them and close them
for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['AVREC_' figlist{ifig} '_BFamp'])
saveas(gcf, ['AVREC_' figlist{ifig} '_BFamp.pdf'])
close (h)
end

%BF - 2
clear peakamp
for ii =1:4
    %cut out groups and stick them side by side in a struct
    peakamp.(AV_m2name{ii}) = AV_m2T.peakamp(AV_m2 == ii);
end
%stick the struct in a table instead to spread it out properly
peakAmpT = struct2table(peakamp); 
%tables can't go into vartestn so we need a matrix, no table2mat so 2 steps:
peakAmpC = table2cell(peakAmpT); PEAKamp = cell2mat(peakAmpC); 
%2 comparisons, ket vs awake and ket vs mus
PEAKamp1 = PEAKamp(:,2:3); %ketamine and awake
PEAKamp2 = PEAKamp(:,2:2:4); %ketamine and muscimol

[AV_Amp_m2_P1,AV_Amp_m2_STATS1] = vartestn(PEAKamp1, 'testtype','BrownForsythe');
[AV_Amp_m2_P2,AV_Amp_m2_STATS2] = vartestn(PEAKamp2, 'testtype','BrownForsythe');

for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['AVREC_' figlist{ifig} '_m2amp'])
saveas(gcf, ['AVREC_' figlist{ifig} '_m2amp.pdf'])
close (h)
end

%% AVREC PEAKLAT

% BF
clear peakLate
for ii =1:4
    %cut out groups and stick them side by side in a struct
    peakLate.(AV_BFname{ii}) = AV_BFT.peaklate(AV_BF == ii);
end
%stick the struct in a table instead to spread it out properly
peakLateT = struct2table(peakLate); 
%tables can't go into vartestn so we need a matrix, no table2mat so 2 steps:
peakLateC = table2cell(peakLateT); PEAKLate = cell2mat(peakLateC); 
%2 comparisons, ket vs awake and ket vs mus
PEAKLate1 = PEAKLate(:,2:3); %ketamine and awake
PEAKLate2 = PEAKLate(:,2:2:4); %ketamine and muscimol

[AV_Late_BF_P1,AV_Late_BF_STATS1] = vartestn(PEAKLate1, 'testtype','BrownForsythe');
[AV_Late_BF_P2,AV_Late_BF_STATS2] = vartestn(PEAKLate2, 'testtype','BrownForsythe');

%get each figure, save them and close them
for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['AVREC_' figlist{ifig} '_BFlate'])
saveas(gcf, ['AVREC_' figlist{ifig} '_BFlate.pdf'])
close (h)
end

%BF - 2
clear peakLate
for ii =1:4
    %cut out groups and stick them side by side in a struct
    peakLate.(AV_m2name{ii}) = AV_m2T.peaklate(AV_m2 == ii);
end
%stick the struct in a table instead to spread it out properly
peakLateT = struct2table(peakLate); 
%tables can't go into vartestn so we need a matrix, no table2mat so 2 steps:
peakLateC = table2cell(peakLateT); PEAKLate = cell2mat(peakLateC); 
%2 comparisons, ket vs awake and ket vs mus
PEAKLate1 = PEAKLate(:,2:3); %ketamine and awake
PEAKLate2 = PEAKLate(:,2:2:4); %ketamine and muscimol

[AV_Late_m2_P1,AV_Late_m2_STATS1] = vartestn(PEAKLate1, 'testtype','BrownForsythe');
[AV_Late_m2_P2,AV_Late_m2_STATS2] = vartestn(PEAKLate2, 'testtype','BrownForsythe');

for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['AVREC_' figlist{ifig} '_m2late'])
saveas(gcf, ['AVREC_' figlist{ifig} '_m2late.pdf'])
close (h)
end

%% RELRES PEAKAMP

% BF
clear peakamp
for ii =1:3
    %cut out groups and stick them side by side in a struct
    peakamp.(RE_BFname{ii}) = RE_BFT.peakamp(RE_BF==ii);
end
%stick the struct in a table instead to spread it out properly
peakAmpT = struct2table(peakamp); 
%tables can't go into vartestn so we need a matrix, no table2mat so 2 steps:
peakAmpC = table2cell(peakAmpT); PEAKamp = cell2mat(peakAmpC); 
%2 comparisons, ket vs awake and ket vs mus
PEAKamp1 = PEAKamp(:,1:2); %ketamine and awake
PEAKamp2 = PEAKamp(:,1:2:3); %ketamine and muscimol

[RE_Amp_BF_P1,RE_Amp_BF_STATS1] = vartestn(PEAKamp1, 'testtype','BrownForsythe');
[RE_Amp_BF_P2,RE_Amp_BF_STATS2] = vartestn(PEAKamp2, 'testtype','BrownForsythe');

%get each figure, save them and close them
for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['RELRES_' figlist{ifig} '_BFamp'])
saveas(gcf, ['RELRES_' figlist{ifig} '_BFamp.pdf'])
close (h)
end

%BF - 2
clear peakamp
for ii =1:3
    %cut out groups and stick them side by side in a struct
    peakamp.(RE_m2name{ii}) = RE_m2T.peakamp(RE_m2 == ii);
end
%stick the struct in a table instead to spread it out properly
peakAmpT = struct2table(peakamp); 
%tables can't go into vartestn so we need a matrix, no table2mat so 2 steps:
peakAmpC = table2cell(peakAmpT); PEAKamp = cell2mat(peakAmpC); 
%2 comparisons, ket vs awake and ket vs mus
PEAKamp1 = PEAKamp(:,1:2); %ketamine and awake
PEAKamp2 = PEAKamp(:,1:2:3); %ketamine and muscimol

[RE_Amp_m2_P1,RE_Amp_m2_STATS1] = vartestn(PEAKamp1, 'testtype','BrownForsythe');
[RE_Amp_m2_P2,RE_Amp_m2_STATS2] = vartestn(PEAKamp2, 'testtype','BrownForsythe');

for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['RELRES_' figlist{ifig} '_m2amp'])
saveas(gcf, ['RELRES_' figlist{ifig} '_m2amp.pdf'])
close (h)
end

%% RELRES PEAKLAT

% BF
clear peakLate
for ii =1:3
    %cut out groups and stick them side by side in a struct
    peakLate.(RE_BFname{ii}) = RE_BFT.peaklate(RE_BF==ii);
end
%stick the struct in a table instead to spread it out properly
peakLateT = struct2table(peakLate); 
%tables can't go into vartestn so we need a matrix, no table2mat so 2 steps:
peakLateC = table2cell(peakLateT); PEAKLate = cell2mat(peakLateC); 
%2 comparisons, ket vs awake and ket vs mus
PEAKLate1 = PEAKLate(:,1:2); %ketamine and awake
PEAKLate2 = PEAKLate(:,1:2:3); %ketamine and muscimol

[RE_Late_BF_P1,RE_Late_BF_STATS1] = vartestn(PEAKLate1, 'testtype','BrownForsythe');
[RE_Late_BF_P2,RE_Late_BF_STATS2] = vartestn(PEAKLate2, 'testtype','BrownForsythe');

%get each figure, save them and close them
for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['RELRES_' figlist{ifig} '_BFlate'])
saveas(gcf, ['RELRES_' figlist{ifig} '_BFlate.pdf'])
close (h)
end

%BF - 2
clear peakLate
for ii =1:3
    %cut out groups and stick them side by side in a struct
    peakLate.(RE_m2name{ii}) = RE_m2T.peaklate(RE_m2 == ii);
end
%stick the struct in a table instead to spread it out properly
peakLateT = struct2table(peakLate); 
%tables can't go into vartestn so we need a matrix, no table2mat so 2 steps:
peakLateC = table2cell(peakLateT); PEAKLate = cell2mat(peakLateC); 
%2 comparisons, ket vs awake and ket vs mus
PEAKLate1 = PEAKLate(:,1:2); %ketamine and awake
PEAKLate2 = PEAKLate(:,1:2:3); %ketamine and muscimol

[RE_Late_m2_P1,RE_Late_m2_STATS1] = vartestn(PEAKLate1, 'testtype','BrownForsythe');
[RE_Late_m2_P2,RE_Late_m2_STATS2] = vartestn(PEAKLate2, 'testtype','BrownForsythe');

for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['RELRES_' figlist{ifig} '_m2late'])
saveas(gcf, ['RELRES_' figlist{ifig} '_m2late.pdf'])
close (h)
end

%% ABSRES PEAKAMP

% BF
clear peakamp
for ii =1:3
    %cut out groups and stick them side by side in a struct
    peakamp.(AB_BFname{ii}) = AB_BFT.peakamp(AB_BF==ii);
end
%stick the struct in a table instead to spread it out properly
peakAmpT = struct2table(peakamp); 
%tables can't go into vartestn so we need a matrix, no table2mat so 2 steps:
peakAmpC = table2cell(peakAmpT); PEAKamp = cell2mat(peakAmpC); 
%2 comparisons, ket vs awake and ket vs mus
PEAKamp1 = PEAKamp(:,1:2); %ketamine and awake
PEAKamp2 = PEAKamp(:,1:2:3); %ketamine and muscimol

[AB_Amp_BF_P1,AB_Amp_BF_STATS1] = vartestn(PEAKamp1, 'testtype','BrownForsythe');
[AB_Amp_BF_P2,AB_Amp_BF_STATS2] = vartestn(PEAKamp2, 'testtype','BrownForsythe');

%get each figure, save them and close them
for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['ABSRES_' figlist{ifig} '_BFamp'])
saveas(gcf, ['ABSRES_' figlist{ifig} '_BFamp.pdf'])
close (h)
end

%BF - 2
clear peakamp
for ii =1:3
    %cut out groups and stick them side by side in a struct
    peakamp.(AB_m2name{ii}) = AB_m2T.peakamp(AB_m2 == ii);
end
%stick the struct in a table instead to spread it out properly
peakAmpT = struct2table(peakamp); 
%tables can't go into vartestn so we need a matrix, no table2mat so 2 steps:
peakAmpC = table2cell(peakAmpT); PEAKamp = cell2mat(peakAmpC); 
%2 comparisons, ket vs awake and ket vs mus
PEAKamp1 = PEAKamp(:,1:2); %ketamine and awake
PEAKamp2 = PEAKamp(:,1:2:3); %ketamine and muscimol

[AB_Amp_m2_P1,AB_Amp_m2_STATS1] = vartestn(PEAKamp1, 'testtype','BrownForsythe');
[AB_Amp_m2_P2,AB_Amp_m2_STATS2] = vartestn(PEAKamp2, 'testtype','BrownForsythe');

for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['ABSRES_' figlist{ifig} '_m2amp'])
saveas(gcf, ['ABSRES_' figlist{ifig} '_m2amp.pdf'])
close (h)
end

%% ABSRES PEAKLAT

% BF
clear peakLate
for ii =1:3
    %cut out groups and stick them side by side in a struct
    peakLate.(AB_BFname{ii}) = AB_BFT.peaklate(AB_BF==ii);
end
%stick the struct in a table instead to spread it out properly
peakLateT = struct2table(peakLate); 
%tables can't go into vartestn so we need a matrix, no table2mat so 2 steps:
peakLateC = table2cell(peakLateT); PEAKLate = cell2mat(peakLateC); 
%2 comparisons, ket vs awake and ket vs mus
PEAKLate1 = PEAKLate(:,1:2); %ketamine and awake
PEAKLate2 = PEAKLate(:,1:2:3); %ketamine and muscimol

[AB_Late_BF_P1,AB_Late_BF_STATS1] = vartestn(PEAKLate1, 'testtype','BrownForsythe');
[AB_Late_BF_P2,AB_Late_BF_STATS2] = vartestn(PEAKLate2, 'testtype','BrownForsythe');

%get each figure, save them and close them
for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['ABSRES_' figlist{ifig} '_BFlate'])
saveas(gcf, ['ABSRES_' figlist{ifig} '_BFlate.pdf'])
close (h)
end

%BF - 2
clear peakLate
for ii =1:3
    %cut out groups and stick them side by side in a struct
    peakLate.(AB_m2name{ii}) = AB_m2T.peaklate(AB_m2 == ii);
end
%stick the struct in a table instead to spread it out properly
peakLateT = struct2table(peakLate); 
%tables can't go into vartestn so we need a matrix, no table2mat so 2 steps:
peakLateC = table2cell(peakLateT); PEAKLate = cell2mat(peakLateC); 
%2 comparisons, ket vs awake and ket vs mus
PEAKLate1 = PEAKLate(:,1:2); %ketamine and awake
PEAKLate2 = PEAKLate(:,1:2:3); %ketamine and muscimol

[AB_Late_m2_P1,AB_Late_m2_STATS1] = vartestn(PEAKLate1, 'testtype','BrownForsythe');
[AB_Late_m2_P2,AB_Late_m2_STATS2] = vartestn(PEAKLate2, 'testtype','BrownForsythe');

for ifig = 1:4
h = gcf;
set(h, 'PaperType', 'A4');
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperUnits', 'centimeters');
savefig(h, ['ABSRES_' figlist{ifig} '_m2late'])
saveas(gcf, ['ABSRES_' figlist{ifig} '_m2late.pdf'])
close (h)
end