function [fid, cleaner] = open_text_file_for_write(filename)
%OPEN_TEXT_FILE_FOR_WRITE Open a text file for writing with clearer errors.

    [parent_dir, ~, ~] = fileparts(filename);
    if ~isempty(parent_dir) && exist(parent_dir, 'dir') ~= 7
        mkdir(parent_dir);
    end

    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open "%s" for writing.', filename);
    end

    cleaner = onCleanup(@() fclose(fid));
end
