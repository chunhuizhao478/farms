function areEqual = compareFiles(file1, file2)
    % Initialize the return value
    areEqual = false;

    % Try to open the first file
    fid1 = fopen(file1, 'r');
    if fid1 == -1
        error('Error opening file: %s', file1);
    end

    % Try to open the second file
    fid2 = fopen(file2, 'r');
    if fid2 == -1
        fclose(fid1);  % Close the first file before exiting
        error('Error opening file: %s', file2);
    end

    % Read the contents of the first file
    content1 = fread(fid1, '*char')';
    fclose(fid1);

    % Read the contents of the second file
    content2 = fread(fid2, '*char')';
    fclose(fid2);

    % Compare the contents of the two files
    areEqual = isequal(content1, content2);
    
    if areEqual
        fprintf('The files are identical.\n');
    else
        fprintf('The files are different.\n');
    end
end
