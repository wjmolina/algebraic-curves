function image_data = load_noise_image(experiment_root, index)
    path = fullfile(experiment_root, "images", "temporary", "img_" + index + ".mat");
    if ~isfile(path)
        % Original external experiment files used this path and variable name.
        path = fullfile(experiment_root, "Images", "Temporary", "img" + index + ".mat");
    end
    input = load(path);
    if isfield(input, 'image_data')
        image_data = input.image_data;
    else
        image_data = input.('I');
    end
end
