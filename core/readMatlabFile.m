function fileContent = readMatlabFile(filePath)
    % 读取MATLAB .m 文件内容并存储为字符串变量
    % 输入参数：
    %   filePath - .m 文件的完整路径（字符串）
    % 输出参数：
    %   fileContent - 存储文件内容的字符串变量
    
    % 检查文件是否存在
    if ~isfile(filePath)
        error('文件不存在：%s', filePath);
    end
    
    % 检查文件扩展名是否为.m
    [~, ~, ext] = fileparts(filePath);
    if ~strcmpi(ext, '.m')
        warning('文件类型可能不是MATLAB脚本（.m文件）');
    end
    
    % 以文本模式读取文件内容
    fid = fopen(filePath, 'r');
    if fid == -1
        error('无法打开文件：%s', filePath);
    end
    
    % 将文件内容读取为字符数组
    fileContent = fread(fid, '*char')';
    fclose(fid);
end