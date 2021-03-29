function [data, comment] = readSBE9hex( dataLines, instHeader,procHeader )
%READSBE19HEX Parses the given data lines from a SBE19 .hex data file. 
%
% Currently, only raw hex (raw voltages and frequencies) output format is 
% supported.
%
% Inputs:
%   dataLines  - Cell array of strings, the lines from the .hex file which
%                contain data.
%   instHeader - Struct containing information contained in the .hex file
%                header.
%
% Outputs:
%   data       - Struct containing variable data.
%   comment    - Struct containing variable comment.
%
% Author:       Paul McCarthy <paul.mccarthy@csiro.au>
% Contributor:  Guillaume Galibert <guillaume.galibert@utas.edu.au>
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



  % boolean variables used to handle all of the optional entries
  
  % Currently, only raw hex output format is supported, 
  % so each pressure reading has an associated temperature 
  % compensation voltage. If raw eng is supported in the 
  % future, the logic defining the value of pressureVolt 
  % will need to be implemented.
  procHeader = procHeader{1};

  % preallocate space for the sample data
  data = struct;
  nLines = length(dataLines);
  preallocZeros = zeros(nLines, 1);
                   data.temperature  = preallocZeros;
                   data.conductivity = preallocZeros;
                   data.pressure = preallocZeros;
                   data.pst = preallocZeros;
  % read in the data
  for k = 1:length(dataLines)
    
    % l is an index into the current line
    l    = 1;
    line = dataLines{k};
    
                     data.temperature (k) = convert_freq(line(l:l+5)); l=l+6;
                     data.conductivity(k) = convert_freq(line(l:l+5)); l=l+6;
                     data.pressure    (k) = convert_freq(line(l:l+5)); l=l+6; 
                     data.pst         (k) = hex2dec(line(57:59));
                     
  end
  
  [data, comment] = convertData(data, procHeader,instHeader);
  
end
function freq = convert_freq(str)

freq = hex2dec(str(1:2))*256+hex2dec(str(3:4))+hex2dec(str(5:end))/256;

end

function [newData, comment] = convertData(data, procHeader,instHeader)
%CONVERTDATA Converts the data contained in the .hex file into IMOS
% compliant parameters.
%
  newData = struct;
  comment = struct;
  
  % temperature, pressure and conductivity will always be 
  % present, and conductivity conversion requires temperature 
  % and pressure, so we manually do all three
  newData.TEMP = convertTemperature(data.temperature, procHeader);
  comment.TEMP = '';
  
  newData.PRES_REL = convertPressure(data.pressure, data.pst, procHeader);
  comment.PRES_REL = '';
  
  newData.CNDC = convertConductivity(data.conductivity,newData.PRES_REL,newData.TEMP, procHeader);
  comment.CNDC = '';
  
  newData.TIME = find_time(instHeader);
  

  [newData.Latitude,newData.Longitude] = find_latlon(instHeader);
      
  
end

function time = find_time(intHeader)
    for i = 1:length(intHeader)
        res = regexp(intHeader{i},'[A-Z]?[a-z]* \d{2} \d{4}  \d{2}:\d{2}:\d{2}','match');
        if ~isempty(res)
            time = res;
        end
    end
end

function [lat1,lon1] = find_latlon(intHeader)
    for i = 1:length(intHeader)
       lat = regexp(intHeader{i},'Latitude = \d{2} \d{2}.\d{2}','match');

       if ~isempty(lat)
            lat = regexp(lat,'\d*','match');
            res = str2double(lat{1}(1))+str2double(lat{1}(2))./60+str2double(lat{1}(3))./3600;
            lat1 = res;
       end
       lon = regexp(intHeader{i},'Longitude = \d* \d{2}.\d{2}','match');
       if ~isempty(lon)
            lon = regexp(lon,'\d*','match');
            res = str2double(lon{1}(1))+str2double(lon{1}(2))./60+str2double(lon{1}(3))./3600;
            lon1 = res;
        end
    end

end



function temperatures = convertTemperature(temperature_freq, procHeader)
%CONVERTTEMPERATURE Converts temperature A/D counts to degrees celsius, via
% the convertion equation provided with SBE19 calibration sheets.
%
f = temperature_freq;
f0 = 1000;
R = f0./f;

pcut = [];
for i =1:length(procHeader)
    z = strfind(procHeader(i),'TemperatureSensor');
    if z{1}>0
        pcut = [pcut,i];
    end
end

if length(pcut) > 0, pcutTemperatureSensorHeader = procHeader(pcut(1):pcut(2)); end
sample  ='<([ABCDGHIJ])>(([-]?\d+.?\d*)[Ee]{1}[+-]{1}\d*)</[ABCDGHIJ]>';
param_temp = struct;
for i = 1:length([pcutTemperatureSensorHeader])
    % catch ABCDGHIJ ...
    res = regexp(pcutTemperatureSensorHeader{i},sample,'tokens');
    if ~isempty(res)
        param_temp = setfield(param_temp,res{1}{1},str2num(res{1}{2}));
    end

end
temperatures = 1 ./ ( param_temp.G + (param_temp.H .* log(R)) + (param_temp.I .* log(R).^2) + (param_temp.J * log(R).^3)) - 273.15;

end

function conductivity = convertConductivity(conductivity, pressure, temperature, procHeader)
%CONVERTCONDUCTIVITY Converts conductivity frequency to siemens per metre, 
% via the convertion equation provided with SBE19 calibration sheets.
% conductivity : conductivity freq
% temperature double
% pressure double
% procheader
f = conductivity;
f = f./1000; % KHz - Hz

pcut = [];
for i =1:length(procHeader)
    z = strfind(procHeader(i),'ConductivitySensor');
    if z{1}>0
        pcut = [pcut,i];
    end
end
if length(pcut) > 0, pcutconductivitySensorHeader = procHeader(pcut(1):pcut(2)); end
sample  ='<([GHIJ])>(([-]?\d+.?\d*)[Ee]{1}[+-]{1}\d*)</[GHIJ]>';
sample2  ='<(CPcor)>(([-]?\d+.?\d*)[Ee]{1}[+-]{1}\d*)</(CPcor)>';
sample3  ='<(CTcor)>(([-]?\d+.?\d*)[Ee]{1}[+-]{1}\d*)</(CTcor)>';
param_con = struct;
for i = 1:length([pcutconductivitySensorHeader])
    % catch ABCDGHIJ ...
    res = regexp(pcutconductivitySensorHeader{i},sample,'tokens');
    if ~isempty(res)
        param_con = setfield(param_con,res{1}{1},str2num(res{1}{2}));
    end
    
    res = regexp(pcutconductivitySensorHeader{i},sample2,'tokens');
    if ~isempty(res)
        param_con = setfield(param_con,res{1}{1},str2num(res{1}{2}));
    end
    res = regexp(pcutconductivitySensorHeader{i},sample3,'tokens');
    if ~isempty(res)
        param_con = setfield(param_con,res{1}{1},str2num(res{1}{2}));
    end

end


t = temperature;
p = pressure;
conductivity = (param_con.G + param_con.H.*(f.^2) + param_con.I.*(f.^3) + param_con.J.*(f.^4))./(10.*(1 + t.*param_con.CTcor + param_con.CPcor.*p));

end

function pressure = convertPressure(pressurefreq, pressuretemp, procHeader)
%CONVERTPRESSURE Converts pressure A/D counts to decibars, via the 
% convertion equation provided with SBE19 calibration sheets. Here, the
% constant value 14.7*0.689476 dbar for atmospheric pressure isn't 
% substracted like in the processed .cnv data.
%

pcut = [];
for i =1:length(procHeader)
    z = strfind(procHeader(i),'PressureSensor');
    if z{1}>0
        pcut = [pcut,i];
    end
end

if length(pcut) > 0, pcutPressureSensorHeader = procHeader(pcut(1):pcut(2)); end
sample  ='<([ABCDGHIJT][12345])>(([-]?\d+.?\d*)[Ee]{1}[+-]{1}\d*)</[ABCDGHIJT][12345]>';
sample1  ='<(Slope)>([-]?\d+.?\d*)</(Slope)>';
sample2  ='<(AD590M)>(([-]?\d+.?\d*)[Ee]{1}[+-]{1}\d*)</(AD590M)>';
sample3  ='<(AD590B)>(([-]?\d+.?\d*)[Ee]{1}[+-]{1}\d*)</(AD590B)>';
sample4  ='<(Offset)>([-]?\d+.?\d*)</(Offset)>';
param_pres = struct;
for i = 1:length([pcutPressureSensorHeader])
    % catch C1  C2 T1 T2 ...
    res = regexp(pcutPressureSensorHeader{i},sample,'tokens');
    if ~isempty(res)
%         param_pres.(res{1}(1))= res{1}{2};
        param_pres = setfield(param_pres,res{1}{1},str2num(res{1}{2}));
    end
    
    % catch 
    res = regexp(pcutPressureSensorHeader{i},sample1,'tokens');
    if ~isempty(res)
        param_pres = setfield(param_pres,res{1}{1},str2num(res{1}{2}));
    end
    res = regexp(pcutPressureSensorHeader{i},sample2,'tokens');
    if ~isempty(res)
        param_pres = setfield(param_pres,res{1}{1},str2num(res{1}{2}));
    end
    res = regexp(pcutPressureSensorHeader{i},sample3,'tokens');
    if ~isempty(res)
        param_pres = setfield(param_pres,res{1}{1},str2num(res{1}{2}));
    end
    res = regexp(pcutPressureSensorHeader{i},sample4,'tokens');
    if ~isempty(res)
        param_pres = setfield(param_pres,res{1}{1},str2num(res{1}{2}));
    end
end
pst = pressuretemp 
freq = pressurefreq


psi2dbar = 0.689476;
Td = param_pres.AD590M * pst + param_pres.AD590B;
c = param_pres.C1 + Td .* (param_pres.C2 + Td .* param_pres.C3);
d = param_pres.D1 + Td .* param_pres.D2;
t0 = param_pres.T1 + Td .* (param_pres.T2 + Td .* (param_pres.T3 + Td .* (param_pres.T4 + Td .* param_pres.T5)));
t0f = 1e-6 * t0 .* freq;
fact = 1 - (t0f .* t0f);
pres = psi2dbar .* (c .* fact.* (1 - d .* fact));
pressure = param_pres.Slope .* pres + param_pres.Offset-10.1353;
end




