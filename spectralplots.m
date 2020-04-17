clear

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

cd(homedir)

disp('loading mat file WT_CL')
tic
load('WT_CL.mat','WT_CL')
toc

% frequencies can be found in wtTable.freq{1} to clarify the following
% rows choices; actual intended rows commented
theta = (49:54);        %(4:7);
alpha = (44:48);        %(8:12);
beta_low = (39:43);     %(13:18);
beta_high = (34:38);    %(19:30);
gamma_low = (26:33);    %(31:60);
gamma_high = (19:25);   %(61:100);

osciName = {'theta' 'alpha' 'beta_low' 'beta_high' 'gamma_low' 'gamma_high'};
osciRows = {theta alpha beta_low beta_high gamma_low gamma_high};
Layer = 'IV';

% loop through stimulation frequencies
for iStim
% loop through layers
curWT = WT_CL(contains(WT_CL.layer,Layer),:);

% loop through oscillation frequencies
curOsci = curWT.scalogram{1,1}(gamma_low,:);
% get power 
curOsci = abs(curOsci).^2;
peak_power = max(max(curOsci));
[row,col] = find(curOsci == peak_power);


    
    
    
    
    
    
    
    
    
    
    
    
    

