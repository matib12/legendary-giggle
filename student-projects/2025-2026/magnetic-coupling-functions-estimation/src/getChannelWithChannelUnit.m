function [channelData,channelUnit] = getChannelWithChannelUnit(varargin)
% EDITED BY LUNGHINI TO GET ALSO CHANNEL UNIT
% SEE GETCHANNEL FUNCTION FOR COMMON GETCHANNEL FUNCTION

%% input parsing

if nargin == 3 % default to raw.ffl
    source = 'raw';
    [channel, gstart, gstop_or_dur] = deal(varargin{:});
elseif nargin == 4
    [source, channel, gstart, gstop_or_dur] = deal(varargin{:});
else
    error('getChannel needs 3 or 4 inputs')
end

if ~strcmp(channel(3), ':') % default to V1 if prefix if missing
    channel = ['V1:', channel];
end

[gstart, ~, dur] = start_stop_dur(gstart, gstop_or_dur);

% protect against asking too much data from raw
if dur > 1200 && strcmp(source, 'raw')
    disp('*** You are probably asking too much data. Enter to continue, CTRL-C to cancel. ***')
    pause
end

% extend 'raw', 'trend', ... to usual files
source = expand_ffl(source);

%% get the data
[channelData,~,~,~,~,~,channelUnit] = frgetvectN(source, channel, gstart, dur, -1); 
channelData=channelData';



