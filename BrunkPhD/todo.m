%% HeadersRus
% would be good practice to document these steps (input/output, any
% important points to note, such as changes needed in certain conditions or
% the need to have a folder empty of all except target .mat files etc.

%% I love headers

try cd('D:\MyCode\Dynamic_CSD_Analysis');
catch
    cd('D:\Dynamic_CSD_Analysis');
end
home = pwd; 
addpath(genpath(home));

%% now I can announce when I do specific things to make them easier to find!

run('dynamic_CSD_Inf_single'); 
%   Input: groups and raw folders. It calculates 
%   and stores CSD information in Data struct which is
%   saved in the DATA folder. 
% 
%   6a and 6b and Inf are detected here and there is single trial data
%   IMPORTANT: DO NOT change sink list here. If you need another set of
%   sinks then create a SEPARATE and UNIQUELY NAMED script. 
%   Note: Currently the graphs being generated are limited to select sinks
%   although all sinks are being stored in the data structure
% 
%   CHANGE if needed: add your working directory to the try/catch; change
%   condition to run 

%% Which group analysis conditions do we want?

 GroupAnalysis_fnc_temoporal (1, 0, 0);
% GroupAnalysis_fnc_temoporal (1, 1, 0);
 GroupAnalysis_fnc_temoporal (1, 0, 1);
% GroupAnalysis_fnc_temoporal (1, 1, 1);
% GroupAnalysis_fnc_temoporal (1, 0);
% GroupAnalysis_fnc_temoporal (0, 1)
% GroupAnalysis_fnc_temoporal (0, 0)

%% Comparisons

cd MichaPhD
% Compare_Groups_Correlation;