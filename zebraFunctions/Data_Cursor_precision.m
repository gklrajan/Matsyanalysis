function outText = Data_Cursor_precision(obj, event)

pos = get(event, 'Position');

if length(pos) == 3
    outText = {sprintf('X: %10.0f', pos(1)), ...
        sprintf('Y: %.5f', pos(2)), ...
        sprintf('Z: %.5f', pos(3))};
    
else
    outText = {sprintf('X: %10.0f', pos(1)), ...
        sprintf('Y: %.5f', pos(2))};
end

end