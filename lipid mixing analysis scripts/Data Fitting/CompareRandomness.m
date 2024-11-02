function [rand_a, rand_b, nmin_a, nmin_b] = CompareRandomness(files_a, files_b, Options)
    % Compute bootstrap CI for randomness parameters and N_min
    % in two paired datasets.
    % Args:
    %    files_a, files_b: cell arrays of filenames for data
    %    options: options string
    % Rets:
    %   rand_a, rand_b: randomness parameters where rand_a(i, :) is the
    %       bootstrap distribution of the randomness parameter for the data
    %       in files_a(i).
    %   nmin_a, nmin_b: similar bootstrap distributions of n_min.
    
    Nboot = 100;
    
    % pre-allocate arrays
    rand_a(length(files_a), Nboot) = 0;
    rand_b(length(files_a), Nboot) = 0;
    nmin_a(length(files_a), Nboot) = 0;
    nmin_b(length(files_a), Nboot) = 0;
    
    for i = 1:length(files_a)
        dat = load(files_a{i});
        if isfield(Options, 'subdat')
            dat = dat.dat;
        end
        dist = Extract_Data(dat, 'Normal CDF-Improved Analysis', [], Options);
        res = bootstrp(Nboot, @calc_rand, dist);
        rand_a(i , :) = res(:, 1);
        nmin_a(i,:) = res(:, 2);
    end
    
    for i = 1:length(files_b)
        dat = load(files_b{i});
        dist = Extract_Data(dat, 'Normal CDF-Improved Analysis', [], Options);
        res = bootstrp(Nboot, @calc_rand, dist);
        rand_b(i, :) = res(:, 1);
        nmin_b(i, :) = res(:, 2);
    end
    
end

function res = calc_rand(waiting_times)
    % returns an array of [RandomnessParameter, Nmin]
    MeanWaitingTimeSquared = (mean(waiting_times))^2;
    MeanOfSquared = mean(waiting_times.^2);
    res(1) = (MeanOfSquared - MeanWaitingTimeSquared)/MeanWaitingTimeSquared;
    res(2) = 1 / res(1);
end
    