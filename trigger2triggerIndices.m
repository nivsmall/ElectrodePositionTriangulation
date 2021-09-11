function [trigger_indices] = trigger2triggerIndices(triggers,ChannelDataTimeStamps)
%TRIGGER2TRIGGERINDICES Summary of this function goes here
trigger_indices=zeros(1, length(triggers));
cnt=1;
for trig = triggers(1,:)
   trigger_indices(cnt)=find(trig==ChannelDataTimeStamps);
   cnt=cnt+1;
end
end

