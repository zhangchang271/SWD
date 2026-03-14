function restoreMatlabFile(fileContent, filePath, overwrite)
    % 将存储的文件内容恢复为MATLAB .m 文件
    % 输入参数：
    %   fileContent - 存储文件内容的字符串变量
    %   filePath - 要保存的.m文件的完整路径（字符串）
    %   overwrite - 可选参数，如果设为true则覆盖现有文件（默认为false）
    
    % 设置默认参数
    if nargin < 3
        overwrite = false;
    end
    
    % 检查文件路径是否以.m结尾
    [fileDir, fileName, ext] = fileparts(filePath);
    if isempty(ext)
        filePath = fullfile(fileDir, [fileName, '.m']);
    elseif ~strcmpi(ext, '.m')
        warning('文件扩展名不是.m，将强制改为.m扩展名');
        filePath = fullfile(fileDir, [fileName, '.m']);
    end
    
    % 检查文件是否已存在
    if isfile(filePath) && ~overwrite
        error('文件已存在：%s。若要覆盖，请设置overwrite参数为true。', filePath);
    end
    
    % 确保目录存在
    if ~isempty(fileDir) && ~isfolder(fileDir)
        mkdir(fileDir);
    end
    
    % 以写入模式打开文件
    fid = fopen(filePath, 'w');
    if fid == -1
        error('无法创建文件：%s', filePath);
    end
    
    % 将内容写入文件
    fwrite(fid, fileContent, 'char');
    fclose(fid);
    
    fprintf('文件已成功恢复至：%s\n', filePath);
end