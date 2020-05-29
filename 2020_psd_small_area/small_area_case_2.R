library(rmutil)
library(gtools)
library(hypergeo)
library(dplyr)
library(quadprog)
library(lpSolve)
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

posterior_wd_sample = function(n, n_tot, eps, dp_sample, prior_a, prior_b) {
  samples = numeric(n)
  t = 1
  d0 = ddl(0, eps)
  while (t <= n) {
    theta_t = rbeta(1, prior_a, prior_b)
    x_t = rbinom(1, n_tot, theta_t)
    p_accept = ddl(dp_sample - x_t, eps) / d0
    if (runif(1) <= p_accept) {
      samples[t] = theta_t
      t = t + 1
    }
  }
  samples
}

postprocess_state = function(df, total) {
  
  counts = cbind(df$wd_priv, df$ws_priv)
  
  k = length(counts)
  frac_hist = solve.QP(
    diag(k), 
    counts,
    t(rbind(rep(1, k), 
            diag(k), 
            -1*diag(k))),
    c(total, rep(0, k), rep(-1*total, k)),
    meq=1
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
  
  pp_counts = frac_hist_trunc + round_pert
  pp_df = data.frame(matrix(pp_counts, ncol=2))
  names(pp_df) <- c("wd_priv_pp", "ws_priv_pp")
  cbind(df, pp_df)
}

create_county_query = function() {
  ndf = adf
  for (priv_col in c("wd", "ws")) {
    ndf[, paste(priv_col, "_priv", sep="")] = ndf[, priv_col] + rdl(dim(ndf)[1], .1)
  }
  postprocess_state(ndf, sum(adf$wp))
}

# simulation -----------------------------------------------------------------

# set RNG seed for replicability
set.seed(20200527)

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

# filter alaska only
adf = mdf[mdf$state == 2, c("fips", "wd", "wp")]
adf$ws = adf$wp - adf$wd
true_nmr = sum(mdf$wd) / sum(mdf$wp)


n_sim = 100

res = data.frame(
  sim_id=c(),
  fips=c(),
  method=c(),
  p1=c()
)

# for each replicate and each county
for (sim_id in 1:n_sim) {
  
  cq = create_county_query() 
  
  for (row_ix in 1:dim(cq)[1]) {
    # sample from posterior distribution
    post_samps = posterior_wd_sample(
      500, cq[row_ix, "wp"], .1, cq[row_ix, "wd_priv"], .1, 10
    )
    
    # posterior power
    post_p1 = mean(post_samps >= true_nmr)
    
    # naive posterior power
    naive_p1 = 1 - pbeta(true_nmr, 
                         .1 + max(cq[row_ix, "wd_priv"], 1), 
                         10 + max(cq[row_ix, "ws_priv"], -9))
    
    # post-processed power
    pp_p1 = 1 - pbeta(true_nmr, 
                      .1 + cq[row_ix, "wd_priv_pp"], 
                      10 + cq[row_ix, "ws_priv_pp"])
    
    res = rbind(
      res,
      data.frame(
        sim_id=rep(sim_id, 3),
        fips=rep(cq[row_ix, "fips"], 3),
        method=c("ConstrainedPosterior", "Naive", "PostProcessed"),
        p1=c(post_p1, naive_p1, pp_p1)
      )
    )
  }
}

true_rej = 1 - pbeta(true_nmr, .1 + adf[, "wd"], 10 + adf[, "ws"])
res = merge(res, data.frame(fips=adf$fips, true_p1=true_rej), by="fips")
fips_names = data.frame(
  fips=c(2020,2090,2100,2110,2122,2130,2150,2170,2180,2201,2220,2240,2261,2280),
  cname=c("Anchorage", "Fairbanks NS", "Haines", "Juneau", "Kenai", "Ketchikan", 
          "Kodiak", "Matanuska", "Nome", "Prince of Wales", "Sitka", 
          "SE Fairbanks", "VC", "WP")
)
res = merge(res, fips_names, by="fips")
res$method = factor(res$method, 
                    levels=c("Naive", "PostProcessed", "ConstrainedPosterior"))

# small county data 
small_counties = res[res$cname %in% c("Nome", "Haines", "Prince of Wales"), ]
small_county_summary = small_counties %>% group_by(cname, method) %>% summarize(
  sample_mse = mean((p1 - true_p1)**2),
  sample_var = var(p1),
  sample_bias2 = mean(p1 - true_p1)**2
)
names(small_county_summary) <- c("County", "Method", "MSE", "Variance", "Bias")
print(xtable(small_county_summary), include.rownames=FALSE)

ggplot(data=small_counties) + 
  geom_histogram(aes(x=p1, y=..density..), bins=20) + 
  geom_vline(aes(xintercept=true_p1), color="blue") + 
  facet_grid(vars(cname), vars(method)) + 
  xlab("Probability of rejecting H0: county mortality rate <= national average") + 
  ylab("Density") +
  ggtitle("Alaskan non-Hispanic White Mortality Rate Hypothesis Testing Posterior Power (Selected Small Counties)")

# all county data
all_county_summary = res %>% group_by(cname, method) %>% summarize(
  sample_mse = mean((p1 - true_p1)**2),
  sample_var = var(p1),
  sample_bias2 = mean(p1 - true_p1)**2
)
names(all_county_summary) <- c("County", "Method", "MSE", "Variance", "Bias")
print(xtable(all_county_summary), include.rownames=FALSE)

ggplot(data=res) + 
  geom_histogram(aes(x=p1, y=..density..), bins=20) + 
  geom_vline(aes(xintercept=true_p1), color="blue") + 
  facet_grid(vars(cname), vars(method)) + 
  xlab("Probability of rejecting H0: county mortality rate <= national average") + 
  ylab("Density") +
  theme(strip.text.y = element_text(size=6)) + 
  ggtitle("Alaskan non-Hispanic White Mortality Rate Hypothesis Testing Posterior Power")

