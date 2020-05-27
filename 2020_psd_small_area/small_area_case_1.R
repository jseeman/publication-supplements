library(rmutil)
library(gtools)
library(hypergeo)
library(dplyr)
library(quadprog)
library(lpSolve)
library(modeest)
library(ggplot2)

ddl = function(v, eps) {
  # PMF of discrete laplace distribution
  p = exp(-eps / 2) 
  (1 - p) / (1 + p) * p**(abs(v))
}
rdl = function(n, eps) {
  # random samples from discrete laplace distribution
  p = 1 - exp(-eps/2)
  rgeom(n, p) - rgeom(n, p)
}

dgdl = function(v, r, eps) {
  # PMF of generalized discrete laplace distribution
  p = exp(-eps / 2) 
  Re(
    exp(lgamma(r + abs(v)) - lgamma(r) - lgamma(1 + abs(v))) * 
      (1 - p)**(2*r) * p**(abs(v)) * 
      hypergeo(r, r + abs(v), 1 + abs(v), p**2)
  )
}

dglcs = function(z, r, delta, eps) {
  # PMF of constrained generalized discrete laplace distribution
  dgdl(as.numeric(z[1]), r, eps) * 
    dgdl(as.numeric(z[2]), r, eps) * 
    dgdl(delta - as.numeric(z[1]) - as.numeric(z[2]), r, eps)
}

l2l1_opt = function(counts, total, total_deaths) {
  # post-process individual state deaths and total counts
  k = length(counts)
  frac_hist = solve.QP(
    diag(k), 
    counts,
    t(rbind(rep(1, k), 
            c(rep(1, k/2), rep(0, k/2)), 
            diag(k), 
            -1*diag(k),
            diag(c(rep(1, k/2), rep(0, k/2))),
            -1*diag(c(rep(1, k/2), rep(0, k/2))))),
    c(total, total_deaths, 
      rep(0, k), rep(-1*total, k),
      rep(0, k), rep(-1*total_deaths, k)),
    meq=2
  )$solution
  
  frac_hist_trunc = trunc(frac_hist)
  
  round_pert = lp(
    "min",
    frac_hist - frac_hist_trunc,
    matrix(rep(1, k), nrow=1), 
    c("="), 
    c(total - sum(frac_hist_trunc)),
    binary.vec=1:k
  )$solution
  
  frac_hist_trunc + round_pert
}


postprocess_state_counts = function(sdf) {
  # post-process privatized joint queries
  num_counties = dim(sdf)[1]
  
  sdf$ws_priv = sdf$wp_priv - sdf$wd_priv
  sdf$hs_priv = sdf$hp_priv - sdf$hd_priv
  sdf$bs_priv = sdf$bp_priv - sdf$bd_priv
  
  stacked_counts = cbind(sdf$wd_priv, sdf$hd_priv, sdf$bd_priv,
                         sdf$ws_priv, sdf$hs_priv, sdf$bs_priv)
  state_total = sum(sdf$tp)
  state_deaths = sum(sdf$td)
  
  pp_counts = l2l1_opt(stacked_counts, state_total, state_deaths)
  res = data.frame(matrix(pp_counts, ncol=6))
  names(res) <- c("wd_priv_pp", "hd_priv_pp", "bd_priv_pp", 
                  "ws_priv_pp", "hs_priv_pp", "bs_priv_pp")
  
  res$wp_priv_pp = res$wd_priv_pp + res$ws_priv_pp
  res$hp_priv_pp = res$hd_priv_pp + res$hs_priv_pp
  res$bp_priv_pp = res$bd_priv_pp + res$bs_priv_pp
  
  cbind(sdf, res[, c("wd_priv_pp", "wp_priv_pp",
                     "hd_priv_pp", "hp_priv_pp",
                     "bd_priv_pp", "bp_priv_pp")])
}

privatize_counts = function(fips_state_id) {
  ndf = mdf[mdf$state == fips_state_id, ]
  for (priv_col in c("wd", "wp", "hd", "hp", "bd", "bp")) {
    ndf[, paste(priv_col, "_priv", sep="")] = ndf[, priv_col] + rdl(dim(ndf)[1], .5)
  }
  postprocess_state_counts(ndf)
}

create_private_query = function(fips_state_id) {
  # resample privatization noise
  ndf = privatize_counts(fips_state_id)
  
  # aggregate to state level
  as.data.frame(
    ndf %>% group_by(state) %>% summarize(
      wd=sum(wd), wp=sum(wp),
      hd=sum(hd), hp=sum(hp),
      bd=sum(bd), bp=sum(bp),
      wd_priv=sum(wd_priv), wp_priv=sum(wp_priv),
      hd_priv=sum(hd_priv), hp_priv=sum(hp_priv),
      bd_priv=sum(bd_priv), bp_priv=sum(bp_priv),
      wd_priv_pp=sum(wd_priv_pp), wp_priv_pp=sum(wp_priv_pp),
      hd_priv_pp=sum(hd_priv_pp), hp_priv_pp=sum(hp_priv_pp),
      bd_priv_pp=sum(bd_priv_pp), bp_priv_pp=sum(bp_priv_pp),
      sd=sum(td), sp=sum(tp), nc=n()
    )
  )
}

posterior_sample_sum_query = function(n, n_tot, r, eps, dp_sample, prior) {
  samples = matrix(0, nrow=n, ncol=3)
  t = 1
  dlta = sum(dp_sample) - n_tot
  d0 = dglcs(c(dlta, 0), r, dlta, eps)
  while (t <= n) {
    theta_t = rdirichlet(1, prior)
    x_t = rmultinom(1, n_tot, theta_t)
    p_accept = dglcs(as.vector(dp_sample - x_t), r, dlta, eps) / d0
    if (!is.na(p_accept) & runif(1) <= p_accept) {
      samples[t, ] = theta_t
      t = t + 1
    }
  }
  samples
}

chi2_test_dist_state = function(mdfs,
                                samps_per_state=10,
                                chi2_samps=10000,
                                eps=.5,
                                empirical_prior_scale=.5) {
  
  pops = as.vector(mdfs[, c("wp_priv", "hp_priv", "bp_priv")])
  pops_tot = mdfs[, "sp"]
  pops_prior = pmax(1, pops) * empirical_prior_scale
  
  pops_samples = posterior_sample_sum_query(
    samps_per_state, 
    n_tot=pops_tot,
    r=mdfs[, "nc"],
    eps=eps,
    dp_sample=pops,
    prior=pops_prior
  )
  
  deaths = as.vector(mdfs[, c("wd_priv", "hd_priv", "bd_priv")])
  deaths_tot = mdfs[, "sd"]
  deaths_prior = pmax(1, deaths) * empirical_prior_scale
  deaths_samples = posterior_sample_sum_query(
    samps_per_state, 
    n_tot=deaths_tot,
    r=mdfs[, "nc"],
    eps=eps,
    dp_sample=deaths,
    prior=deaths_prior
  )
  
  t = 1
  chi2s = numeric(chi2_samps)
  while (t <= chi2_samps) {
    # sample population by race
    s_pop_theta = pops_samples[sample(1:samps_per_state, 1), ]
    s_pop = rmultinom(1, pops_tot, s_pop_theta)
    
    # sample deaths by race
    s_death_theta = deaths_samples[sample(1:samps_per_state, 1), ]
    s_death = rmultinom(1, deaths_tot, s_death_theta)
    
    # reject if constraint not satisfied
    if ( all(s_death <= s_pop)) {
      s_surv = s_pop - s_death
      s_death_exp = nmr * s_pop
      s_surv_exp = (1 - nmr) * s_pop
      
      chi2s[t] = (
        sum((s_surv - s_surv_exp)**2 / s_surv_exp) + 
          sum((s_death - s_death_exp)**2 / s_death_exp)
      )
      t = t + 1
    } 
  }
  chi2s
}

chi2_state_sample_stats = function(mdfs) {
  vs = chi2_test_dist_state(mdfs)
  
  true_exp = rbind(nmr*mdfs[, c("wp", "hp", "bp")],
                   (1-nmr)*mdfs[, c("wp", "hp", "bp")])
  
  true_obs = rbind(mdfs[, c("wd", "hd", "bd")],
                   -mdfs[, c("wd", "hd", "bd")] + 
                     mdfs[, c("wp", "hp", "bp")])
  
  true_chi2 = sum((true_obs - true_exp)**2 / true_exp)
  
  naive_exp = rbind(nmr*mdfs[, c("wp_priv", "hp_priv", "bp_priv")],
                    (1-nmr)*mdfs[, c("wp_priv", "hp_priv", "bp_priv")])
  
  naive_obs = rbind(mdfs[, c("wd_priv", "hd_priv", "bd_priv")],
                    -1 * mdfs[, c("wd_priv", "hd_priv", "bd_priv")] + 
                      mdfs[, c("wp_priv", "hp_priv", "bp_priv")])
  
  naive_chi2 = sum((naive_obs - naive_exp)**2 / naive_exp)
  
  pp_exp = rbind(nmr*mdfs[, c("wp_priv_pp", "hp_priv_pp", "bp_priv_pp")],
                 (1-nmr)*mdfs[, c("wp_priv_pp", "hp_priv_pp", "bp_priv_pp")])
  
  pp_obs = rbind(mdfs[, c("wd_priv_pp", "hd_priv_pp", "bd_priv_pp")],
                 -1 *mdfs[, c("wd_priv_pp", "hd_priv_pp", "bd_priv_pp")] + 
                   mdfs[, c("wp_priv_pp", "hp_priv_pp", "bp_priv_pp")])
  pp_chi2 = sum((pp_obs - pp_exp)**2 / pp_exp)
  
  c(
    mlv(vs, method="meanshift")[1], 
    true_chi2,
    naive_chi2,
    pp_chi2
  )
}


# simulation -----------------------------------------------------------------

# set RNG seed for replicability
set.seed(20200527)
n_sim = 100

# load original data
path_dir = "/Users/jeremyseeman/Desktop/small_area_draft/"
dh = read.csv(paste(path_dir, "deaths_hispanic.csv", sep=""))
dnhw = read.csv(paste(path_dir, "deaths_nhwhite.csv", sep=""))
dnhb = read.csv(paste(path_dir, "deaths_nhblack.csv", sep=""))

mdf = merge(
  merge(
    dnhw[, c("fips", "nhwhite_deaths", "nhwhite_population")],
    dh[, c("fips", "hispanic_deaths", "pop_cdc")],
    by="fips", all=TRUE
  ),
  dnhb[, c("fips", "nh_black_deaths", "nh_black_population")],
  by="fips", all=TRUE
)
names(mdf) = c("fips", "wd", "wp", "hd", "hp", "bd", "bp")
mdf[is.na(mdf)] = 0
mdf$state = floor(mdf$fips / 1000)
mdf$td = mdf$hd + mdf$bd + mdf$wd 
mdf$tp = mdf$hp + mdf$bp + mdf$wp

nmr = c(
  sum(mdf$wd) / sum(mdf$wp),
  sum(mdf$hd) / sum(mdf$hp),
  sum(mdf$bd) / sum(mdf$bp)
)

tres = data.frame(
  t(replicate(n_sim, chi2_state_sample_stats(create_private_query(2))))
)
names(tres) <- c("PosteriorMode", "TrueValue", "Naive", "PostProcessed")

chi2_true = tres[1, "TrueValue"]
chi2_summary = as.data.frame(
  stack(tres[, c("PosteriorMode", "Naive", "PostProcessed")]) %>%
    group_by(ind) %>% summarize(
      SampleMSE=mean((values - chi2_true)**2),
      SampleVar=var(values),
      SampleBias2=mean(values - chi2_true)**2
    )
)

chi2_summary[, c("SampleMSE", "SampleVar", "SampleBias2")] = (
  round(chi2_summary[, c("SampleMSE", "SampleVar", "SampleBias2")], 0)
)
print(xtable(chi2_summary), include.rownames=FALSE)

ggplot(stack(tres[, c("PosteriorMode", "Naive", "PostProcessed")])) + 
  geom_histogram(aes(x=values, y=..density..), bins=20) + 
  geom_vline(xintercept=tres[1, "TrueValue"], color="red") + 
  facet_grid(rows=vars(ind)) + 
  xlab("Replicated Chi-Square test statistics") + 
  ylab("Density") + 
  ggtitle("Test statistics for independence of (mortality x race) in Alaska")

