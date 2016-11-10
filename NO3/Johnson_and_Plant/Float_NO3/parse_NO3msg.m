function spec = parse_NO3msg(isus_file)

% PURPOSE: 
%   This function parses an APEX or NAVIS biochemical float nitrate *.isus
%   file. A structure is returned containing the seawater UV intensity
%   spectrum and other data used to determine the nitrate concentration
%
% USAGE:
%	data = parse_NO3msg(isus_file)
%
% INPUTS:
%	isus_file  = Path\calibration file name as a string
%
% OUTPUTS:
%   spec =   a structure of data and calibration info
%       spec.SDN        = sample time, Matlab sdn        (array)
%      	spec.P          = sample pressure, dbar          (array)
%     	spec.T          = sample temperature, C          (array)
%     	spec.S          = sample salinity, pss           (array)
%      	spec.DC         = sample DC intensity, counts    (scalar)
%     	spec.UV_INTEN   = Measured UV intensities        (matrix)
%      	spec.SWDC       = sea water DC intensity, counts (scalar)
%                         NaN if doesn't exist (earlier floats)
%
%      	spec.spectra_pix_range = pixel registration for the wavelengths in
%                                the calibration file
%       spec.WL_Fit_Range      = wavelength bounds of fit window - used in
%                                SUNA floats to determine spectra_pix_range
%    	spec.pix_fit_win       = pixel range for Multiple linear regression.
%                                Used to subset sample spectra to fit window
%
%     	spec.CalTemp    = calibration temperature this is the temperature
%                         the instrument was calibrated at in the lab.
%    	spec.CalDateStr = lab calibration date string
%     	spec.CalDate    = serial date number for the lab cal date
%
%
% EXAMPLE:
%   jp = parse_NO3msg('\\atlas\chemwebdata\floats\f5143\5143.005.isus')
%
% Created 01/13/2016 by jp

% isus_file ='\\atlas\chemwebdata\floats\f5143\5143.005.isus'; % for testing
% isus_file ='\\atlas\Chem\ISUS\Argo\8497Hawaii\8497.042.isus';
% isus_file ='C:\temp\9095.068.isus';
% isus_file ='C:\temp\6403.113.isus';
% isus_file ='C:\temp\6967.225.isus';
% isus_file ='C:\temp\6976.152.isus';
%isus_file = 'C:\temp\7622.002.isus';
% ************************************************************************
% PARSE *.isus FILE
% ************************************************************************
% FORMATS TO PARSE BASE 10 PART OF *.ISUS MESSAGE FILE
% HEX WILL BE DONE ON THE FLY BELOW AS HEX BLOCK CAN CHANGE SIZE
% FORMATS TO PARSE BASE 10 PART OF *.ISUS MESSAGE FILE
% HEX WILL BE DONE ON THE FLY BELOW AS HEX BLOCK CAN CHANGE SIZE
d_format = '%s %s %s %f %f %f %f %f %f %f'; % six cols, 1st is string
d_format = [d_format,'%f %f %f %f %f %f %f %f %f %f'];
d_format = [d_format,'%f %f %f %f %s %f'];  

%Predim output if no spectra present
spec.SDN        = [];
spec.P          = [];
spec.T          = [];
spec.S          = [];
spec.DC         = [];
spec.UV_INTEN   = [];
spec.SWDC       = [];

spec.CalTemp           = NaN;
spec.CalDate           = NaN;
spec.CalDateStr        = '';
spec.pix_fit_win       = NaN;
spec.spectra_pix_range = NaN;
spec.file_name         = '';
spec.float             = '';
spec.cast              = '';


% check for valid target
if exist(isus_file,'file')
    fid   = fopen(isus_file);
    spec.file_name = regexp(isus_file,'\d+\.\d+\.\D+','once','match');
    if ~isempty(spec.file_name)
        flt_info = regexp(spec.file_name,'\.','split');
        spec.float = flt_info{1};
        spec.cast  = flt_info{2};
        clear flt_info
    end
else
    disp(['Can not find *.isus file at: ',isus_file])
    return
end

% ************************************************************************
% FIND # OF SAMPLE SPECTRA DATA ROWS & HEX BLOCK LENGTHS
% FOR SAMPLE COMPLETENESS
% ************************************************************************
tline   = ' ';
data_ct = 0;
hdr_ct  = 0;
hdr_chk  = 0;
hex_len = [];
while ischar(tline) 
    if  regexp(tline,'^0x','once') % DATA LINES
        hdr_chk =1;
        data_ct = data_ct + 1;
        % GET HEX STRING FROM TLINE
        hex_str = regexp(tline,'\w{150,}','once','match'); % >= 150 char
        if ~isempty(hex_str)
            hex_len(data_ct) = length(hex_str);
        else
            hex_len(data_ct) = 0;
        end
    end
    if regexp(tline,'^H,','once') & hdr_chk == 0% DATA LINES
        hdr_ct = hdr_ct+1;
    end
    tline = fgetl(fid); % get 1st line
end

% Check for header lines
if hdr_ct ==  0
    disp(['File exists but no header lines found in file for: ', ...
        isus_file]);
    fclose(fid);
    return
end
        
% check for data lines
if data_ct == 0
    disp(['File exists but no data lines found for : ',isus_file]);
    fclose(fid);
    return
else
    hex_length      = mode(hex_len); % most common value
    spec_length     = hex_length /4;
    line_test       = hex_len ~= hex_length; % 1's are bad lines
    data_ct         = sum(~line_test); % Sum good lines
    
    spec.SDN        = ones(data_ct,1)*NaN; % predim
    spec.P          = spec.SDN;
    spec.T          = spec.SDN;
    spec.S          = spec.SDN;
    spec.DC         = spec.SDN;
    spec.SWDC       = spec.SDN;
    spec.UV_INTEN   = ones(data_ct,spec_length)*NaN; % predim

    %hex_length   = length(dline{1,25}{1})/4;
    hex_template = '%04x';
    hex_format   = hex_template;
    % BUILD HEX FORMAT STRING % all hex #'s = 4 char
    for i = 1:spec_length  - 1 % build hex format string
        hex_format = [hex_format,hex_template];
    end
end

% ************************************************************************
% GO TO BEGINING OF FILE AND START PARSING THE DATA FILE
% ************************************************************************
frewind(fid); 
tline = fgetl(fid); % prime the engine - get 1st line 
line_ct = 0;
hdr_chk = 0;
while ischar(tline)
    % GET HEADER LINE INFO
    if  hdr_chk == 0 && ~isempty(regexp(tline,'^H','once'))
        hdr_ct = hdr_ct +1;
        if regexp(tline,'Calibration\s+Date','once') % calibration date
            % use to match to cal file calibration date
            d_str = textscan(tline, '%*s %*s %s', 'Delimiter',',',...
                'CollectOutput',1); 
            spec.CalDateStr = d_str{1,1}{1};
            spec.CalDate = datenum(d_str{1,1},'mm/dd/yyyy');
        end        
        if regexp(tline,'Sw\s+Calibration\s+Temp','once') % cal temp
            % use to match to cal file calibration temp
            spec.CalTemp = str2double(regexp(tline,'\d+\.\d+','match')); 
        end
        if regexp(tline,'Wavelength Fit','once') % Get WL MLR (SUNA)
            % This can be used to create pix_fit_win for SUNA files
            spec.WL_fit_win = str2double(regexp(tline,'\d+\.*\d*','match')); 
        end
        if regexp(tline,'Pixel Fit','once') % Get pixel range for MLR 
            % This will be used to subset sample spectra to fit window
            % for MLR - eventually have user control to adjust
            spec.pix_fit_win = str2double(regexp(tline,'\d+','match')); 
        end
    end
    


    % PARSE DATA LINES
    if  ~isempty(regexp(tline,'^0x','once')) ; % data line
        dline = textscan(tline,d_format,'Delimiter',',');
        if ~isempty(dline{1,25})
            % CHECK SIZE OF HEX STRING, IF PARTIAL SPECTRA MOVE TO NEXT LINE
            if length(dline{1,25}{1}) == hex_length % right size
                line_ct = line_ct +1;
                
                spec.SDN(line_ct)  = datenum(dline{1,3},'mm/dd/yyyy HH:MM:SS');
                spec.P(line_ct)    = dline{1,5};
                spec.T(line_ct)    = dline{1,6};
                spec.S(line_ct)    = dline{1,7};
                spec.DC(line_ct)   = dline{1,18}; % dark current
                
                % EARLIER FLOATS DID NOT RETURN SW DC AFTER HEX BLOCK so SET
                % SW DC TO DC IF NEED BE
                if isempty(dline{1,26})
                    spec.SWDC(line_ct) = spec.DC(line_ct); % Older model
                else % used if shutter stuck open
                    spec.SWDC(line_ct) = dline{1,26};
                end
                
                hex_UV_INTEN       = dline{1,25}{1}; % hex string still
                UV_INTEN = (sscanf(hex_UV_INTEN, hex_format))'; %col to row too
                spec.UV_INTEN(line_ct,1:spec_length) = UV_INTEN;
                if line_ct == 1 % get some one time info
                    hdr_chk = 1; % header before and after only need once
                    
                    % this range is the pixel registration for WL's in cal file
                    spec.spectra_pix_range = [dline{23},dline{24}];
                end
            end
        end
    end
    tline = fgetl(fid);
end
fclose(fid);
clearvars -except spec









