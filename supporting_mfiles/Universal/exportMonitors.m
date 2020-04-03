function exportMonitors(MonitorList)
global DSSText;

for i = 1:length(MonitorList)
    DSSText.Command = ['export mon ',char(MonitorList(i))];
end


