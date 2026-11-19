function [qp_e, qp, pixelized_noisy] = vetterli_reconstruction(binary_image, x_a, x_b, y_a, y_b, x_n, y_n, target_snr_db)
    x_range = linspace(x_a, x_b, x_n);
    y_range = linspace(y_a, y_b, y_n);
    [pixelized, m_x, m_y] = pixelization(binary_image, x_range, y_range, 4);

    if target_snr_db ~= - 1
        pixelized_noisy = imnoise(pixelized, 'gaussian', 0, 10 ^ (- target_snr_db / 10) * mean(mean(pixelized .^ 2)));
    else
        pixelized_noisy = pixelized;
    end

    repo_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    coefficient_file = fullfile(repo_root, 'auxiliary', 'vetterli', 'moment_coefs', 'coefs_n_4_m_7_max_power_6_x_max_16.mat');
    load(coefficient_file, 'coef_g', 'coef_dg');
    [moments_dx_y, moments_x_dy, ~] = moment_generator(pixelized_noisy, m_x, m_y, coef_g, coef_dg);
    moment_matrix = matrix4nullspace(4, moments_dx_y, moments_x_dy);

    for row_ind = 1 : size(moment_matrix, 1)
        row_norm = norm(moment_matrix(row_ind, :));
        if row_norm > 0
            moment_matrix(row_ind, :) = moment_matrix(row_ind, :) / row_norm;
        end
    end

    x_power = 0;
    y_power = 0;
    x_power_vec = zeros(1, 15);
    y_power_vec = zeros(1, 15);
    term = 0;

    while x_power + y_power <= 4
        term = term + 1;
        x_power_vec(term) = x_power;
        y_power_vec(term) = y_power;

        if y_power == 0
            y_power = x_power + y_power + 1;
            x_power = 0;
        else
            y_power = y_power - 1;
            x_power = x_power + 1;
        end
    end

    poly_eval_mat = zeros(length(m_x) * length(m_y), 15);

    for m_x_ind = 1 : length(m_x)
        for m_y_ind = 1 : length(m_y)
            row_ind = m_y_ind + (m_x_ind - 1) * length(m_y);
            poly_eval_mat(row_ind, :) = (m_x(m_x_ind) .^ x_power_vec) .* (m_y(end - m_y_ind + 1) .^ y_power_vec);
        end
    end

    % The polynomial grid assigns the largest y coordinate to the first image row.
    b_eq = sign(pixelized_noisy(fliplr(m_y) == 0, m_x == 0) - 0.5);
    a_ineq = poly_eval_mat;
    b_ineq = ones(length(m_x) * length(m_y), 1) * inf;

    for m_x_ind = 1 : length(m_x)
        for m_y_ind = 1 : length(m_y)
            if pixelized_noisy(m_y_ind, m_x_ind) >= min(0.7, max(max(pixelized_noisy)))
                row_ind = m_y_ind + (m_x_ind - 1) * length(m_y);
                b_ineq(row_ind) = 0;
                a_ineq(row_ind, :) = - a_ineq(row_ind, :);
            elseif pixelized_noisy(m_y_ind, m_x_ind) <= 0.3
                row_ind = m_y_ind + (m_x_ind - 1) * length(m_y);
                b_ineq(row_ind) = 0;
            end
        end
    end

    constraints = a_ineq(isfinite(b_ineq), :);
    poly_coefs_qp = sign_constrained_fit(moment_matrix, poly_eval_mat, constraints, b_eq);
    binary_image_qp = ac_image(poly_coefs_qp, x_range, y_range);
    qp = binary_image_qp;
    pixelized_qp = pixelization(binary_image_qp , x_range , y_range , 4);
    max_iterations = 20;
    delta = 1e-4;
    current_err = 1;
    new_err = 0;
    iter = 0;

    while (((current_err - new_err) / current_err) > 0.001 ) && (iter < max_iterations)
        iter = iter + 1;
        gradient_matrix = zeros(length(m_x) * length(m_y), 15);

        for coef_ind = 1 : 15
            new_coef = poly_coefs_qp;
            new_coef(coef_ind) = new_coef(coef_ind) + delta;
            binary_image_new = ac_image(new_coef , x_range, y_range);
            pixelized_new = pixelization(binary_image_new, x_range, y_range, 4);
            aux = (pixelized_new - pixelized_qp) / delta;
            gradient_matrix(:, coef_ind) = aux(:);
        end

        diff = pixelized_qp - pixelized_noisy;
        diff = diff(:);
        current_err = norm(diff);
        change = - pinv(gradient_matrix) * diff;
        binary_image_qp = ac_image(poly_coefs_qp + change, x_range, y_range);
        pixelized_qp = pixelization(binary_image_qp, x_range, y_range, 4);
        new_err = sum(sum((pixelized_qp - pixelized_noisy) .^ 2)) .^ 0.5;

        while new_err > current_err * 1.001
            change = change / 2;
            binary_image_qp = ac_image(poly_coefs_qp + change, x_range, y_range);
            pixelized_qp = pixelization(binary_image_qp, x_range, y_range, 4);
            new_err = sum(sum((pixelized_qp - pixelized_noisy) .^ 2)) .^ 0.5;
        end

        poly_coefs_qp = poly_coefs_qp + change;
    end

    qp_e = binary_image_qp;
end
function coefficients = sign_constrained_fit(moments, grid, constraints, constant)
    % Eliminate the fixed constant and scale the remaining monomials on the sampling grid.
    [~, r] = qr(grid(:, 2:end), 0);
    design = moments(:, 2:end) / r;
    offset = constant * moments(:, 1);
    signs = constraints(:, 2:end) / r;
    limits = -constant * constraints(:, 1);
    lengths = sqrt(sum(signs .^ 2, 2));
    keep = lengths > 0;
    if any(limits(~keep) < 0)
        error('The sign at the origin contradicts the fixed coefficient.');
    end
    signs = signs(keep, :) ./ lengths(keep);
    limits = limits(keep) ./ lengths(keep);

    % An interior starting point simplifies the first steps of the active-set fit.
    n = size(design, 2);
    [start, ~] = linprog([zeros(n, 1); -1], [signs, ones(size(signs, 1), 1)], limits, ...
        [], [], [-inf(n, 1); 0], [inf(n, 1); 1]);
    if isempty(start) || any(~isfinite(start))
        z = zeros(n, 1);
    else
        z = start(1:n);
    end
    if any(signs * z - limits > 1e-7)
        % Some LP solvers return a finite but infeasible point. Project onto
        % violated half-spaces before using it as an active-set starting point.
        z = zeros(n, 1);
        for iteration = 1:20000
            [violation, index] = max(signs * z - limits);
            if violation <= 1e-10
                break
            end
            z = z - violation * signs(index, :)';
        end
        if any(signs * z - limits > 1e-7)
            error('The polynomial sign constraints have no verified feasible starting point.');
        end
    end
    active = false(size(limits));
    for iteration = 1:500
        rows = signs(active, :);
        if isempty(rows)
            directions = eye(n);
        else
            directions = null(rows);
        end
        residual = design * z + offset;
        if isempty(directions)
            step = zeros(n, 1);
        else
            step = -directions * (pinv(design * directions) * residual);
        end
        slopes = signs * step;
        candidates = find(~active & slopes > 1e-12 * max(1, norm(step)));
        fractions = (limits(candidates) - signs(candidates, :) * z) ./ slopes(candidates);
        [fraction, index] = min(fractions);
        if ~isempty(fraction) && fraction < 1 - 1e-10
            z = z + max(0, fraction) * step;
            active(candidates(index)) = true;
            continue
        end
        z = z + step;
        gradient = design' * (design * z + offset);
        multipliers = [];
        if ~isempty(rows)
            multipliers = -(rows' \ gradient);
        end
        [lowest, index] = min(multipliers);
        if isempty(lowest) || lowest >= -1e-10 * max(norm(gradient), eps)
            coefficients = [constant; r \ z];
            if any(~isfinite(coefficients)) || any(constraints * coefficients > 1e-6)
                error('Polynomial fit violates the sign constraints.');
            end
            return
        end
        indices = find(active);
        active(indices(index)) = false;
    end
    error('The sign-constrained polynomial fit did not converge.');
end
