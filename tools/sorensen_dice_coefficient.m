function dsc = sorensen_dice_coefficient(image_data, reconstruction)
    if ~isequal(size(image_data), size(reconstruction))
        error('sorensen_dice_coefficient:SizeMismatch', 'Both masks must have the same dimensions.');
    end
    tp = sum(sum((image_data == reconstruction) .* image_data));
    fp = sum(sum((image_data ~= reconstruction) .* reconstruction));
    fn = sum(sum((image_data ~= reconstruction) .* image_data));
    
    dsc = 2 * tp / (2 * tp + fp + fn);
end