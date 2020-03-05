function [DATA]=attenuation10(DATA)

Stimlist = [DATA.EPOCHS.code]; %all codes in order of presentation
if length(Stimlist) == 996
    Stimlist = Stimlist(1:995);
end

TempLFP = []; %to fill with each leaf of the LFP needed
TempEPOCHS = DATA.EPOCHS(1); %to copy all fields
i2 = 0; %to build the LFP from 1 to 350 or however many there end up to be

%NOTE: code numbers 1 through 7 are attenuated at 10dB
for i1 = 1:length(Stimlist)%cycle through all codes in order of presentation
    if Stimlist(i1) > 0 && Stimlist(i1) < 8 %will check for each number to match and then add matrix to TempLFP
        i2=i2+1;
        TempLFP(:,:,i2) = DATA.LFP(:,:,i1);
        TempEPOCHS(i2) = DATA.EPOCHS(i1);
    end
end

DATA.LFP = TempLFP;
DATA.EPOCHS = TempEPOCHS;