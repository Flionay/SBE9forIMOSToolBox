function panel = tabbedPane( parent, tabs, tabNames, useButtons )
%TABBEDPANE Creates a tabbed pane containing the given tabs (uipanels).
%
% Creates a panel which contains a row of buttons along the top (labelled
% with the given tabNames); when the user pushes a button, the
% corresponding tab is displayed in the panel.
%
% Inputs:
%   parent      - Figure or uipanel to be used as the parent.
%   tabs        - Vector of uipanel handles.
%   tabNames    - Cell array of tab names, the same length as the tabs vector.
%   useButtons  - Optional. Boolean value. If true, each tab is displayed by 
%                 pushing a button. If false, each tab is displayed by 
%                 selecting from a drop-down menu instead.
%
% Outputs:
%   panel       - Handle to the tabbed pane uipanel.
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
narginchk(3,4);

if ~ishandle(parent),    error('parent must be a graphics handle');          end
if ~isvector(tabs) ||...
    isempty(tabs)  ||...
   ~any(ishandle(tabs)), error('tabs must be a vector of graphics handles'); end
if ~iscellstr(tabNames), error('tabNames must be a cell array of strings');  end
if length(tabs) ~= length(tabNames)
                         error('tabs and tabNames must be the same length'); end
if nargin == 3, useButtons = true;
elseif ~islogical(useButtons), error('useButtons must be logical');          end

  % it's up to the caller to set the position
  panel = uipanel(...
    'Parent',     parent,...
    'Units',      'normalized',...
    'BorderType', 'none');
  
  % set tab positions and make all tabs invisible
%   set(tabs,...
%     'Parent',   panel,...
%     'Units',    'normalized',...
%     'Position', [0.0, 0.0, 1.0, 0.95],...
%     'Visible',  'off');
  set(tabs,...
    'Parent',   panel,...
    'Units',    'normalized',...
    'Position', posUi2(panel, 100, 1, 6:100, 1, 0),...
    'Visible',  'off');
  
%   set(panel,   'Units', 'pixels');
%   set(tabs,    'Units', 'pixels');

  % create tab button row/popup menu along top
  if useButtons
    
    numTabs = length(tabNames);
    buttons = nan(numTabs, 1);
    for k = 1:numTabs
%       buttons(k) = uicontrol(...
%         'Parent',   panel,...
%         'Style',    'pushbutton',...
%         'String',   tabNames{k},...
%         'Units',    'normalized',...
%         'Position', [(k-1)/numTabs, 0.95, 1/numTabs, 0.05],...
%         'Callback', @tabCallback);
      buttons(k) = uicontrol(...
        'Parent',   panel,...
        'Style',    'pushbutton',...
        'String',   tabNames{k},...
        'Units',    'normalized',...
        'Position', posUi2(panel, 100, numTabs, 1:5, k, 0),...
        'Callback', @tabCallback);
    end
    
%     set(buttons, 'Units', 'pixels');
    
    tabCallback(buttons(1), []);
    
  %popup menu instead of buttons
  else
%     menu = uicontrol(...
%       'Parent',   panel,...
%       'Style',    'popupmenu',...
%       'String',   tabNames,...
%       'Value',    1,...
%       'Units',    'normalized',...
%       'Position', [0.0, 0.95, 1.0, 0.05],...
%       'Tag',      'exportPopUpMenu',...
%       'Callback', @tabCallback);
    menu = uicontrol(...
      'Parent',   panel,...
      'Style',    'popupmenu',...
      'String',   tabNames,...
      'Value',    1,...
      'Units',    'normalized',...
      'Position', posUi2(panel, 100, 1, 1:5, 1, 0),...
      'Tag',      'exportPopUpMenu',...
      'Callback', @tabCallback);
  
%     set(menu, 'Units', 'pixels');
    
    tabCallback(menu, []);
  end
  
  function tabCallback(source,ev)
  % Called when one of the tab buttons is clicked. Sets all but the selected 
  % tabs invisible.
  
    if useButtons, idx = find(buttons == source); 
    else           idx = get(source, 'Value');
    end
    
    set(tabs,      'Visible', 'off');
    set(tabs(idx), 'Visible', 'on');
  end
end