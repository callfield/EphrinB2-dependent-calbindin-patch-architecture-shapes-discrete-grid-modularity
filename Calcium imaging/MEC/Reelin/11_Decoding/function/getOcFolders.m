function folderPaths = getOcFolders(baseDir)
% Return paths to subfolders inside folders whose names start with 'Oc'.

ocFolders = dir(fullfile(baseDir, "Oc*"));
folderPaths = {};

for i = 1:numel(ocFolders)
    if ocFolders(i).isdir
        ocPath = fullfile(baseDir, ocFolders(i).name);

        subFolders = dir(ocPath);
        subFolders = subFolders([subFolders.isdir]);
        subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));

        for j = 1:numel(subFolders)
            folderPaths{end + 1} = fullfile(ocPath, subFolders(j).name);
        end
    end
end

end
