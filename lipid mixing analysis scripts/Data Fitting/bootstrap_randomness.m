function [pvals, ci] = bootstrap_randomness(time_lists)
% calculate bootstrapped randomness parameter intervals for several waiting time lists
n_boot = 100;
n_sam = length(time_lists);
figure; hold on
for sam_idx = 1:n_sam
  nmin_dist(sam_idx, :) = bootstrp(n_boot, @calc_nmin, time_lists{sam_idx});
  ci(sam_idx, 1) = prctile(nmin_dist(sam_idx, :), 5);
  ci(sam_idx, 2) = prctile(nmin_dist(sam_idx, :), 95);
  histogram(nmin_dist(sam_idx, :));
end
alphaval = 0.05 / ((n_sam * n_sam - n_sam) / 2)
for i=1:n_sam
  for j=1:n_sam
    [h, p] = ttest2(nmin_dist(i,:), nmin_dist(j,:), 'Alpha', alphaval, 'Vartype', 'unequal');
    pvals(i, j) = p;
  end
end
end

function nmin = calc_nmin(waiting_times)
    % returns an array of [RandomnessParameter, Nmin]
    MeanWaitingTimeSquared = (mean(waiting_times))^2;
    MeanOfSquared = mean(waiting_times.^2);
    rand_par  = (MeanOfSquared - MeanWaitingTimeSquared)/MeanWaitingTimeSquared;
    nmin = 1 / rand_par;
end
