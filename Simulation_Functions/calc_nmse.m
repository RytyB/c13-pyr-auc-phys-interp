function [nmse] = calc_nmse(data, fitted)
    % Input: Two matrices representing the pyruvate and lactate dynamic signal
    % Output: NMSE between the two matrices

    resid = data - fitted;
    resid = resid .^ 2;
    rss = sum(resid, 'all');

    nmse = rss / sum( data .^2, 'all' );
end

