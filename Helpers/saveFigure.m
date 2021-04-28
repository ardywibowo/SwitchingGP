function saveFigure(directory, fileName)
%SAVEFIGURE Summary of this function goes here
%   Detailed explanation goes here

createFolder(directory);

set(gcf, 'visible','off', 'CreateFcn', 'set(gcf,''visible'',''on'')');
saveas(gcf, [directory fileName '.fig']);
saveas(gcf, [directory fileName '.png']);
close;

end

