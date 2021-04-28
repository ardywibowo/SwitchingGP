function createFolder(baseDirectory, folderNames)
%CREATEFOLDER Summary of this function goes here
%   Detailed explanation goes here

folderExists = 7;

if nargin < 2
	if ~(exist(baseDirectory, 'dir') == folderExists)
		mkdir(pwd, baseDirectory);
	end	
	return;
end

for folderName = folderNames
	folderChar = char(folderName);
	if ~(exist([baseDirectory, folderChar], 'dir') == folderExists)
		mkdir(pwd, [baseDirectory, folderChar]);
	end	
end

end

