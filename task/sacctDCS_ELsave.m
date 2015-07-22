function sacctDCS_ELsave(xp,EDFname,iBlock,iLeg)

% save ET-file
fullEDFname = sprintf('%s%s_%s_%s_%s_block%i_%s.edf', xp.dataFolder, xp.codename, xp.subject, xp.tDCS, xp.legNames{iLeg}, iBlock, datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
Eyelink('CloseFile');
Eyelink('WaitForModeReady', 500); % helps to avoid errors in retrieving files
try
    status = Eyelink('ReceiveFile', EDFname, fullEDFname); %this collects the file from the eyelink
    disp(status);
    disp(['File ' fullEDFname ' saved to disk']);
catch
    warning(['File ' fullEDFname ' not saved to disk']);
end