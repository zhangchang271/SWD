function [uxt,trace_index,trace_mask,distance] = offset_traces(seismo_v,gx,sx,dx,offset,side)
%OFFSET_TRACES Select receiver traces by physical source-receiver distance.
%
% gx and sx are model-grid coordinates.  dx converts their difference to
% the physical distance used by the Radon transform and the weighting step.
% trace_mask always refers to the original receiver order.  For the left
% side, uxt and trace_index are returned from near to far offset, matching
% the convention used by RTlADx and ADWDgrad_w.

distance_all = (gx(:).'-sx)*dx;

switch lower(side)
    case 'right'
        trace_mask = distance_all > 0 & distance_all <= offset;
        trace_index = find(trace_mask);
    case 'left'
        trace_mask = distance_all < 0 & -distance_all <= offset;
        trace_index = find(trace_mask);
        trace_index = fliplr(trace_index);
    otherwise
        error('offset_traces:InvalidSide', ...
            'side must be ''right'' or ''left''.');
end

distance = abs(distance_all(trace_index));
uxt = seismo_v(:,trace_index);
end
