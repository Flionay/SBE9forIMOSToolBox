function sample_data = SBE9Parse( filename, mode )
%SBE26PARSE Parses a .tid data file from a Seabird SBE26
% TP logger.
%
% This function is able to read in a .tid data file retrieved
% from a Seabird SBE26 Temperature and Pressure Logger. It is 
% assumed the file consists in the following columns:
%
%   - measurement number
%   - date and time (mm/dd/yyyy HH:MM:SS) of beginning of measurement
%   - pressure in psia
%   - temperature in degrees Celsius
%
% Inputs:
%   filename    - cell array of files to import (only one supported).
%   mode        - Toolbox data type mode.
%
% Outputs:
%   sample_data - Struct containing sample data.
%
% Author:       Guillaume Galibert <guillaume.galibert@utas.edu.au>
%

%
% Copyright (C) 2017, Australian Ocean Data Network (AODN) and Integrated 
% Marine Observing System (IMOS).
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.
% If not, see <https://www.gnu.org/licenses/gpl-3.0.en.html>.
%
narginchk(1,2);

if ~iscellstr(filename)
    error('filename must be a cell array of strings');
end

% only one file supported currently
filename = filename{1};

  [dataLines, instHeaderLines, procHeaderLines] = readSBEcnv( filename, mode );
  
  % read in the raw instrument header
  instHeader = parseInstrumentHeader(instHeaderLines, mode);
  
  % replace procHeaderLines from .XMLCON file --- aim to CTD SBE 9
  file_path = split(filename,'.');
  path_xml = [file_path{1},'.XMLCON'];
    
  try
    fid = fopen(path_xml, 'rt');
    procHearder = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
  catch 
      error(['The file ',path_xml,'if not exist']);
  end
  procHeaderLines = procHearder{1};
  procHeader = parseProcessedHeader(procHeaderLines);
  
  % use the appropriate subfunction to read in the data
  % assume that anything with a suffix not equal to .hex
  % is a .cnv file
  [~, ~, ext] = fileparts(filename);
  
  if strcmpi(ext, '.hex')
      [data, comment] = readSBE9hex(dataLines, instHeaderLines, procHearder);
  else
      [data, comment] = readSBEcnvData(dataLines, instHeader, procHeader, mode);
  end
  
  % create sample data struct,
  % and copy all the data in
  sample_data = struct;
  
  sample_data.toolbox_input_file    = filename;
  sample_data.meta.featureType      = mode;
  sample_data.meta.instHeader       = instHeader;
  sample_data.meta.procHeader       = procHeader;
  



time        =  genTimestamps(instHeader, procHearder,data); % we assume a measurement is averaged over 4 minutes
pressure    = data.PRES_REL; % 1psi = 0.6894757dbar
temperature = data.TEMP;
cndc = data.CNDC;

% create sample data struct,
% and copy all the data in




sample_data.meta.instrument_make = 'Seabird';
sample_data.meta.instrument_model = 'SBE9';

sample_data.meta.instrument_firmware = '';

sample_data.meta.instrument_serial_no = '';

% sample_data.meta.instrument_sample_interval = median(diff(time*24*3600));

sample_data.dimensions = {};
sample_data.variables  = {};

% generate time data from header information
sample_data.dimensions{1}.name              = 'TIME';
sample_data.dimensions{1}.typeCastFunc      = str2func(netcdf3ToMatlabType(imosParameters(sample_data.dimensions{1}.name, 'type')));
sample_data.dimensions{1}.data              = sample_data.dimensions{1}.typeCastFunc(time);
sample_data.dimensions{1}.comment           = 'Time stamp corresponds to the centre of the measurement which lasts 4 minutes.';

sample_data.variables{end+1}.name           = 'TIMESERIES';
sample_data.variables{end}.typeCastFunc     = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end}.data             = sample_data.variables{end}.typeCastFunc(1);
sample_data.variables{end}.dimensions       = [];
sample_data.variables{end+1}.name           = 'LATITUDE';
sample_data.variables{end}.typeCastFunc     = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end}.data             = sample_data.variables{end}.typeCastFunc(data.Latitude);
sample_data.variables{end}.dimensions       = [];
sample_data.variables{end+1}.name           = 'LONGITUDE';
sample_data.variables{end}.typeCastFunc     = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end}.data             = sample_data.variables{end}.typeCastFunc(data.Longitude);
sample_data.variables{end}.dimensions       = [];
sample_data.variables{end+1}.name           = 'NOMINAL_DEPTH';
sample_data.variables{end}.typeCastFunc     = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end}.data             = sample_data.variables{end}.typeCastFunc(NaN);
sample_data.variables{end}.dimensions       = [];

% create a variable for each parameter
coordinates = 'TIME LATITUDE LONGITUDE NOMINAL_DEPTH';
    
% dimensions definition must stay in this order : T, Z, Y, X, others;
% to be CF compliant
sample_data.variables{end+1}.dimensions     = 1;
sample_data.variables{end  }.name           = 'DEPTH';
sample_data.variables{end  }.typeCastFunc   = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end  }.data           = sample_data.variables{end}.typeCastFunc(pressure);
sample_data.variables{end  }.coordinates    = coordinates;
% let's document the constant pressure atmosphere offset previously
% applied by SeaBird software on the absolute presure measurement
% sample_data.variables{end}.applied_offset   = sample_data.variables{end}.typeCastFunc(-14.7*0.689476);

sample_data.variables{end+1}.dimensions     = 1;
sample_data.variables{end  }.name           = 'TEMP';
sample_data.variables{end  }.typeCastFunc   = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end  }.data           = sample_data.variables{end}.typeCastFunc(temperature);
sample_data.variables{end  }.coordinates    = coordinates;

sample_data.variables{end+1}.dimensions     = 1;
sample_data.variables{end  }.name           = 'CNDC';
sample_data.variables{end  }.typeCastFunc   = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end  }.data           = sample_data.variables{end}.typeCastFunc(cndc);
sample_data.variables{end  }.coordinates    = coordinates;

end

function header = parseInstrumentHeader(headerLines, mode)
%PARSEINSTRUMENTHEADER Parses the header lines from a SBE19/37 .cnv file.
% Returns the header information in a struct.
%
% Inputs:
%   headerLines - cell array of strings, the lines of the header section.
%   mode        - Toolbox data type mode.
%
% Outputs:
%   header      - struct containing information that was in the header
%                 section.
%
header = struct;

% there's no real structure to the header information, which
% is annoying. my approach is to use various regexes to search
% for info we want, and to ignore everything else. inefficient,
% but it's the nicest way i can think of

headerExpr   = '^\*\s*(SBE \S+|SeacatPlus)\s+V\s+(\S+)\s+SERIAL NO.\s+(\d+)';
%BDM (18/2/2011) - new header expressions to reflect newer SBE header info
headerExpr2  = '<HardwareData DeviceType=''(\S+)'' SerialNumber=''(\S+)''>';
headerExpr3  = 'Sea-Bird (.*?) *?Data File\:';
scanExpr     = 'number of scans to average = (\d+)';
scanExpr2    = '*\s+ <ScansToAverage>(\d+)</ScansToAverage>';
memExpr      = 'samples = (\d+), free = (\d+), casts = (\d+)';
sampleExpr   = ['sample interval = (\d+) (\w+), ' ...
    'number of measurements per sample = (\d+)'];
sampleExpr2  ='*\s+ <Samples>(\d+)</Samples>';
profExpr     = '*\s+ <Profiles>(\d+)</Profiles>';
modeExpr     = 'mode = (\w+), minimum cond freq = (\d*), pump delay = (\d*)';
pressureExpr = 'pressure sensor = (strain gauge|quartz)';
voltExpr     = 'Ext Volt ?(\d+) = (yes|no)';
outputExpr   = 'output format = (.*)$';
castExpr     = ['(?:cast|hdr)\s+(\d+)\s+' ...
    '(\d+ \w+ \d+ \d+:\d+:\d+)\s+'...
    'samples (\d+) to (\d+), (?:avg|int) = (\d+)'];
%Replaced castExpr to be specific to NSW-IMOS PH NRT
%Note: also replace definitions below in 'case 9'
%BDM 24/01/2011
castExpr2    = 'Cast Time = (\w+ \d+ \d+ \d+:\d+:\d+)';
intervalExpr = 'interval = (.*): ([\d\.\+)$';
sbe38Expr    = 'SBE 38 = (yes|no), Gas Tension Device = (yes|no)';
optodeExpr   = 'OPTODE = (yes|no)';
voltCalExpr  = 'volt (\d): offset = (\S+), slope = (\S+)';
otherExpr    = '^\*\s*([^\s=]+)\s*=\s*([^\s=]+)\s*$';
firmExpr     = '<FirmwareVersion>(\S+)</FirmwareVersion>';
firmExpr2    = '^\*\s*FirmwareVersion:\s*(\S+)'; %SBE39plus
sensorId     = '<Sensor id=''(.*\S+.*)''>';
sensorType   = '<[tT]ype>(.*\S+.*)</[tT]ype>';
serialExpr   = '^\*\s*SerialNumber:\s*(\S+)'; %SBE39plus
serialExpr2  = '^\*\s*SEACAT PROFILER\s*V(\S+)\s*SN\s*(\S+)'; %SEACAT PROFILER

exprs = {...
    headerExpr   headerExpr2    headerExpr3    scanExpr     ...
    scanExpr2    memExpr      sampleExpr   ...
    sampleExpr2  profExpr       modeExpr     pressureExpr ...
    voltExpr     outputExpr   ...
    castExpr     castExpr2   intervalExpr ...
    sbe38Expr    optodeExpr   ...
    voltCalExpr  otherExpr ...
    firmExpr     sensorId   sensorType firmExpr2 serialExpr serialExpr2};



for k = 1:length(headerLines)
    filt = regexp(headerLines{k},'* System UTC =',"split" );
    if length(filt)>1
        
        date = regexp(cell2str(filt(2)),'\w*',"match");
    end
end

month = date{1};
day = str2num(date{2});
year = str2num(date{3});
hour = str2num(date{4});
minute = str2num(date{5});
second = str2num(date{6});


month_str  = {'January','February','March','April','May','June','July','August','September','October','November','December'};
month_st ={'Jan','Feb','Mar','Apr','May','Jun','Jul', 'Aug','Sept','Oct','Nov','Dec'};

for i =1:12
    if strcmp(month, month_str{i}) | strcmp(month,month_st{i})
        month_num = i;
        break
    end
end
header.castDate = datenum(year,month_num,day,hour,minute,second);


for k = 1:length(headerLines)
    
    % try each of the expressions
    for m = 1:length(exprs)
        
        % until one of them matches
        tkns = regexp(headerLines{k}, exprs{m}, 'tokens');
        if ~isempty(tkns)
            
            % yes, ugly, but easiest way to figure out which regex we're on
            switch m
                
                % header
                case 1
                    if ~isfield(header, 'instrument_model')
                        header.instrument_model = tkns{1}{1};
                    end
                    header.instrument_firmware  = tkns{1}{2};
                    header.instrument_serial_no = tkns{1}{3};
                    
                % header2
                case 2
                    if ~isfield(header, 'instrument_model')
                        header.instrument_model = tkns{1}{1};
                    end
                    header.instrument_serial_no = tkns{1}{2};
                    
                % header3
                case 3
                    header.instrument_model     = strrep(tkns{1}{1}, ' ', '');
                    
                % scan
                case 4
                    header.scanAvg = str2double(tkns{1}{1});
                    
                % scan2
                case 5
                    header.scanAvg = str2double(tkns{1}{1});
                    %%ADDED by Loz
                    header.castAvg = header.scanAvg;
                    
                % mem
                case 6
                    header.numSamples = str2double(tkns{1}{1});
                    header.freeMem    = str2double(tkns{1}{2});
                    header.numCasts   = str2double(tkns{1}{3});
                    
                % sample
                case 7
                    header.sampleInterval        = str2double(tkns{1}{1});
                    header.mesaurementsPerSample = str2double(tkns{1}{2});
                    
                % sample2
                case 8
                    header.castEnd = str2double(tkns{1}{1});
                
                % profile
                case 9
                    header.castNumber = str2double(tkns{1}{1});    
                    
                % mode
                case 10
                    header.mode         = tkns{1}{1};
                    header.minCondFreq  = str2double(tkns{1}{2});
                    header.pumpDelay    = str2double(tkns{1}{3});
                    
                % pressure
                case 11
                    header.pressureSensor = tkns{1}{1};
                    
                % volt
                case 12
                    for n = 1:length(tkns),
                        header.(['ExtVolt' tkns{n}{1}]) = tkns{n}{2};
                    end
                    
                % output
                case 13
                    header.outputFormat = tkns{1}{1};
                    
                % cast
                case 14
                    if ~isfield(header, 'castStart')
                        header.castNumber = str2double(tkns{1}{1});
                        header.castDate   = datenum(   tkns{1}{2}, 'dd mmm yyyy HH:MM:SS');
                        header.castStart  = str2double(tkns{1}{3});
                        header.castEnd    = str2double(tkns{1}{4});
                        header.castAvg    = str2double(tkns{1}{5});
                    else
                        % in timeSeries mode we only need the first occurence
                        % but in profile mode we require all cast dates
                        if strcmpi(mode, 'profile')
                            header.castNumber(end+1) = str2double(tkns{1}{1});
                            header.castDate(end+1)   = datenum(   tkns{1}{2}, 'dd mmm yyyy HH:MM:SS');
                            header.castStart(end+1)  = str2double(tkns{1}{3});
                            header.castEnd(end+1)    = str2double(tkns{1}{4});
                            header.castAvg(end+1)    = str2double(tkns{1}{5});
                        end
                    end
                    
                % cast2
                case 15                    
                    header.castDate   = datenum(tkns{1}{1}, 'mmm dd yyyy HH:MM:SS');
                    
                % interval
                case 16
                    header.resolution = tkns{1}{1};
                    header.interval   = str2double(tkns{1}{2});
                    
                % sbe38 / gas tension device
                case 17
                    header.sbe38 = tkns{1}{1};
                    header.gtd   = tkns{1}{2};
                    
                % optode
                case 18
                    header.optode = tkns{1}{1};
                    
                % volt calibration
                case 19
                    header.(['volt' tkns{1}{1} 'offset']) = str2double(tkns{1}{2});
                    header.(['volt' tkns{1}{1} 'slope'])  = str2double(tkns{1}{3});
                    
                % name = value
                case 20
                    header.(genvarname(tkns{1}{1})) = tkns{1}{2};
                    
                %firmware version
                case 21
                    header.instrument_firmware  = tkns{1}{1};
                
                %sensor id
                case 22
                    if ~isfield(header, 'sensorIds')
                        header.sensorIds = {};
                    end
                    header.sensorIds{end+1}  = tkns{1}{1};
                    
                %sensor type
                case 23
                    if ~isfield(header, 'sensorTypes')
                        header.sensorTypes = {};
                    end
                    header.sensorTypes{end+1}  = tkns{1}{1};

                %FirmwareVersion, SBE39plus cnv
                case 24
                    header.instrument_firmware  = tkns{1}{1};
                    
                % SerialNumber, SBE39plus cnv
                case 25
                    header.instrument_serial_no = tkns{1}{1};
                    
                % old SEACAT PROFILER serial number format
                % example "* SEACAT PROFILER V2.1a SN 597   10/15/11  10:02:56.721"
                case 26
                    % is tkns{1}{1} firmware version?
                    header.instrument_serial_no = tkns{1}{2};

            end
            break;
        end
    end
end
end

function header = parseProcessedHeader(headerLines)
%PARSEPROCESSEDHEADER Parses the data contained in the header added by SBE
% Data Processing. This includes the column layout of the data in the .cnv 
% file. 
%
% Inputs:
%   headerLines - Cell array of strings, the lines in the processed header 
%                 section.
%
% Outputs:
%   header      - struct containing information that was contained in the
%                 processed header section.
%

  header = struct;
  header.columns = {};
  
  nameExpr = 'name \d+ = (.+):';
  nvalExpr = 'nvalues = (\d+)';
  badExpr  = 'bad_flag = (.*)$';
  %BDM (18/02/2011) - added to get start time
  startExpr = 'start_time = (\w+ \d+ \d+ \d+:\d+:\d+)';
  volt0Expr = 'sensor \d+ = Extrnl Volt  0  (.+)';
  volt1Expr = 'sensor \d+ = Extrnl Volt  1  (.+)';
  volt2Expr = 'sensor \d+ = Extrnl Volt  2  (.+)';
  binExpr   = 'binavg_binsize = (\d+)';
  
  for k = 1:length(headerLines)
    
    % try name expr
    tkns = regexp(headerLines{k}, nameExpr, 'tokens');
    if ~isempty(tkns)
      header.columns{end+1} = tkns{1}{1};
      continue; 
    end
    
    % then try nvalues expr
    tkns = regexp(headerLines{k}, nvalExpr, 'tokens');
    if ~isempty(tkns)
      header.nValues = str2double(tkns{1}{1});
      continue;
    end
    
    % then try bad flag expr
    tkns = regexp(headerLines{k}, badExpr, 'tokens');
    if ~isempty(tkns)
      header.badFlag = str2double(tkns{1}{1});
      continue;
    end
    
    %BDM (18/02/2011) - added to get start time
    % then try startTime expr
    tkns = regexp(headerLines{k}, startExpr, 'tokens');
    if ~isempty(tkns)
      header.startTime = datenum(tkns{1}{1}, 'mmm dd yyyy HH:MM:SS');
      continue;
    end
    
    % then try volt exprs
    tkns = regexp(headerLines{k}, volt0Expr, 'tokens');
    if ~isempty(tkns)
      header.volt0Expr = tkns{1}{1};
      continue;
    end
    tkns = regexp(headerLines{k}, volt1Expr, 'tokens');
    if ~isempty(tkns)
      header.volt1Expr = tkns{1}{1};
      continue;
    end
    tkns = regexp(headerLines{k}, volt2Expr, 'tokens');
    if ~isempty(tkns)
      header.volt2Expr = tkns{1}{1};
      continue;
    end
    
    % then try bin expr
    tkns = regexp(headerLines{k}, binExpr, 'tokens');
    if ~isempty(tkns)
      header.binSize = str2double(tkns{1}{1});
      continue;
    end
  end
end


function time = genTimestamps(instHeader, procHeader,data)
%GENTIMESTAMPS Generates timestamps for the data. Horribly ugly. I shouldn't 
% have to have a function like this, but the .cnv files do not necessarily 
% provide timestamps for each sample.
%
  
  % time may have been present in the sample 
  % data - if so, we don't have to do any work
%   if isfield(data, 'TIME')
%       time = data.TIME;
%       return;
%   end
  
  % To generate timestamps for the CTD data, we need to know:
  %   - start time
  %   - sample interval
  %   - number of samples
  %
  % The SBE19 header information does not necessarily provide all, or any
  % of this information. .
  %
  start    = 0;
  interval = 0.25;
    
  % figure out number of samples by peeking at the 
  % number of values in the first column of 'data'
  f = fieldnames(data);
  nSamples = length(data.(f{1}));
  
  % try and find a start date - use castDate if present
  if isfield(instHeader, 'castDate')
    start = instHeader.castDate;
  end
  
  % if castStart is present then it means we have several cast records
  if isfield(instHeader, 'castStart')
      time = NaN(1, instHeader.castEnd(end));
      for i=1:length(instHeader.castNumber)
          for j=instHeader.castStart(i):instHeader.castEnd(i)
              time(j) = instHeader.castDate(i) + (j-instHeader.castStart(i))*instHeader.castAvg(i)/(3600*24);
          end
      end
      return;
  end
  
  % if scanAvg field is present, use it to determine the interval
  if isfield(instHeader, 'scanAvg')
    interval = (0.25 * instHeader.scanAvg) / 86400;
  end
  
  % if one of the columns is 'Scan Count', use the 
  % scan count number as the basis for the timestamps 
  if isfield(data, 'ScanCount')
    time = ((data.ScanCount - 1) ./ 345600) + cStart;
  % if scan count is not present, calculate the 
  % timestamps from start, end and interval
  else
    time = (start:interval:start + (nSamples - 1) * interval)';
  end
end
