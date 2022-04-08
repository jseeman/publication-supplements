# 
# ARSIA 2022 Statistical Data Privacy Case Study Figures
# 
# User notes:
# * Data is available at https://data.world/data-society/game-of-thrones
#   in the file "character-predictions_pose.csv"
# * Conda environment specification available in requirements.yml
#

# ----------------------------------------------------------------------------------
# Libraries and setup 

library(tidyverse)
library(gridExtra)
library(boot)
library(MASS)
library(extraDistr)
library(svglite)

options(dplyr.summarise.inform = FALSE)

# Note: change this section to local directory containing 
setwd("/var/git/publication-supplements/2022_arsia/")

# ----------------------------------------------------------------------------------
# Data loading and pre-processing 

df = read.csv("character-predictions_pose.csv") %>%
  # standardize and shorten string variables
  # (helps with conjugated cultures, ex: Braavos vs. Braavosi)
  mutate(title=ifelse(
    title == "", "NONE", 
    str_replace_all(
      str_replace_all(
        str_to_upper(title), 
        "\\h+", ""
      ),
      regex("\\W+"), ""
    )
  ),
  culture=ifelse(
    culture == "", "NONE", 
    str_sub(
      str_replace_all(
        str_to_upper(culture), "\\h+", ""
      ),
      start=1,
      end=4
    )
  ),
  house=ifelse(
    house == "", "NONE", 
    str_sub(
      str_replace_all(
        str_to_upper(house),
        "HOUSE ",
        ""
      ), 
      start=1,
      end=6
    )
  ),
  gender=ifelse(male == 1, "M", "F"),
  nobility=ifelse(isNoble == 1, "NOBLE", "PEASANT"),
  alive=ifelse(isAlive == 1, "ALIVE", "DEAD")) %>%
  dplyr::select(name,
                title,
                culture,
                house,
                gender,
                nobility,
                alive)

# create reduced nominal representations
top_cultures = names(table(df$culture)[table(df$culture) >= 10])
top_titles = names(table(df$title)[table(df$title) >= 10])
top_houses = names(table(df$house)[table(df$house) >= 10])
df = df %>% mutate(
  culture_reduced = ifelse(
    culture %in% top_cultures, culture, "OTHER"
  ),
  title_reduced = ifelse(
    title %in% top_titles, title, "OTHER"
  ),
  house_reduced = ifelse(
    house %in% top_houses, house, "OTHER"
  )
)
n_char = dim(df)[1]


# ----------------------------------------------------------------------------------
# Figure 2: k-anonymity and t-closeness analysis

k_anon = function(s, quasis) {
  # 
  #' calculate k-anonymity of first s rows aggregated
  #' 
  #' @param s number of database rows to aggregate
  #' @param quasis quasi-identifiers (vector of strings)
  #' 
  df[1:s, ] %>% 
    group_by_at(vars(quasis)) %>% 
    summarise(nq=n()) %>%
    pull(nq) %>%
    min
}

k1 = sapply(1:n_char, k_anon, quasis=c("title_reduced"))
k2 = sapply(1:n_char, k_anon, quasis=c("title_reduced", "nobility"))
k3 = sapply(1:n_char, k_anon, quasis=c("culture_reduced"))
k4 = sapply(1:n_char, k_anon, quasis=c("culture_reduced", "nobility"))

sdc_risk_k = ggplot(data=data.frame(s=rep(1:n_char, 4),
                                    QuasiIdentifier=c(
                                      rep("TitleReduced", n_char),
                                      rep("TitleReduced+Nobility", n_char),
                                      rep("CultureReduced", n_char),
                                      rep("CultureReduced+Nobility", n_char)
                                    ),
                                    logk=log10(1 + c(k1, k2, k3, k4))) %>%
                      filter(s >= 20)) + 
  geom_line(aes(x=s, y=logk, color=QuasiIdentifier)) + 
  xlab("Database size") + 
  ylab("log10(1 + k)")

ggsave("count_sdc_risk_k_anon.svg", 
       plot=sdc_risk_k,
       device="svg", width=5, height=3, limitsize=FALSE)
ggsave("count_sdc_risk_k_anon.png", 
       plot=sdc_risk_k,
       device="png", width=5, height=3, limitsize=FALSE)

t_close = function(s, quasis, sens) {
  # 
  #' calculate L1 t-closeness of first s rows aggregated
  #' by comparing the distributions of sens within quasis
  #' 
  #' @param s number of database rows to aggregate
  #' @param quasis quasi-identifiers (vector of strings)
  #' @param sens sensitive attribute (string)
  #' 
  
  # first, construct global distribution of sens
  global_sens = df[1:s, ] %>% 
    group_by_at(vars(sens)) %>%
    summarise(global_sens_freq=n() / s)
  
  # then, construct the local distributions of sens
  quasi_sens = df[1:s, ] %>%
    group_by_at(vars(quasis, sens)) %>%
    summarise(quasi_sens=n())
  
  quasi_freq = df[1:s, ] %>%
    group_by_at(vars(quasis)) %>%
    summarise(quasi_tot=n())
  
  # rejoin the frames above, renormalize, and 
  # find the max L1 distance
  quasi_sens %>% 
    inner_join(quasi_freq, 
               by=quasis) %>%
    mutate(quasi_freq=quasi_sens / quasi_tot) %>%
    inner_join(global_sens, 
               by=sens) %>%
    mutate(sens_dist = abs(quasi_freq - global_sens_freq)) %>%
    group_by_at(vars(quasis)) %>% 
    summarise(sens_dist = sum(sens_dist)) %>%
    pull(sens_dist) %>%
    max
}

tclose_ga = sapply(1:250, t_close, quasis=c("gender"), sens=c("alive"))
tclose_na = sapply(1:250, t_close, quasis=c("nobility"), sens=c("alive"))
tclose_gna = sapply(1:250, t_close, quasis=c("gender", "nobility"), sens=c("alive"))

sdc_risk_tclose = ggplot(data=data.frame(s=rep(1:250, 3),
                                         QuasiIDSens=c(rep("Gender | Alive", 250),
                                                       rep("Nobility | Alive", 250),
                                                       rep("Gender + Nobility | Alive", 250)),
                                         tclose = c(tclose_ga, tclose_na, tclose_gna))) + 
  geom_line(aes(x=s, y=tclose, color=QuasiIDSens), alpha=.8) + 
  xlab("Database size") + 
  ylab("L1 t-closeness")
ggsave("count_sdc_risk_t_close.svg", 
       plot=sdc_risk_tclose,
       device="svg", width=5, height=3, limitsize=FALSE)
ggsave("count_sdc_risk_t_close.png",
       plot=sdc_risk_tclose,
       device="png", width=5, height=3, limitsize=FALSE)

# ----------------------------------------------------------------------------------
# Figure 3: k-anonymity and exact inference degradation

k_anon_conf_degrade = function(s) {
  #' 
  #' Construct k-anonymous Fishers exact test privacy-adjusted
  #' p-values for H0 as specified in the case study main text
  #' at different ks per database size
  #' 
  #' @param s database size
  #' 
  
  # first, construct the underlying contingency table and sample from the null
  conf_table = as.matrix(table(df$nobility[1:s], df$alive[1:s]))
  null_res = data.frame(
    an=rhyper(10000, 
              sum(conf_table[, 1]), 
              sum(conf_table[, 2]), 
              sum(conf_table[1, ])
    )
  )
  # next, for each sample, calculate their implied k-anonymities 
  # and the odds ratios (to compare to the observed ratio)
  obs_or = (conf_table[1, 1] / conf_table[1, 2]) / (conf_table[2, 1] / conf_table[2, 2])
  null_res$dn = sum(conf_table[1, ]) - null_res$an
  null_res$ap = sum(conf_table[, 1]) - null_res$an
  null_res$dp = sum(conf_table[, 2]) - null_res$dn
  null_res$k = pmin(null_res$an, 
                    null_res$dn, 
                    null_res$ap, 
                    null_res$dp)
  null_res$or = (null_res$an / null_res$dn) / (null_res$ap / null_res$dp)
  
  # calculate the privacy-adjusted p-value
  max_k = max(null_res$k)
  k_vs_fisher_pv = sapply(1:max_k, function(k) {
    mean(null_res[null_res$k >= k, "or"] <= obs_or)
  }
  )
  data.frame(
    DBsize=rep(paste0("n=", s), max_k),
    k=1:max_k,
    fpv = k_vs_fisher_pv,
    obs_pv = rep(k_vs_fisher_pv[1], max_k)
  )
}

sdc_util_exact = ggplot(data=rbind(
  k_anon_conf_degrade(200),
  k_anon_conf_degrade(300),
  k_anon_conf_degrade(400),
  k_anon_conf_degrade(500)
)) + geom_line(aes(k, fpv, color=DBsize)) + 
  geom_hline(aes(yintercept=obs_pv, color=DBsize), alpha=.5, linetype="dashed") + 
  xlab("k-anonymity") + 
  ylab("Privacy-adjusted p-value")

ggsave("count_sdc_util_exacttest.svg", 
       plot=sdc_util_exact, 
       device="svg", width=5, height=3, limitsize=FALSE)
ggsave("count_sdc_util_exacttest.png", 
       plot=sdc_util_exact, 
       device="png", width=5, height=3, limitsize=FALSE)
sdc_util_exact

# ----------------------------------------------------------------------------------
# Figure 5: DP count inference

# observed statistics from the data
x_obs = sum(df$alive == "DEAD")
test_xs = 350:650
ltxs = length(test_xs)

ddl = function(v, eps) {
  #' PMF of discrete laplace distribution
  #' 
  #' @param v evaluation points
  #' @param eps PLB
  #' 
  p = exp(-eps) 
  (1 - p) / (1 + p) * p**(abs(v))
}


ddlap_cond_xv = function(v, eps, xv) {
  #'
  #' discrete laplace conditional density
  #' 
  #' @param v evaluation points
  #' @param eps PLB
  #' @param xv confidential value
  #' 
  p = exp(-eps)
  (1 - p) / (1 + p) * p**(abs(v - xv))
}

expmech_cond_xv = function(v, eps, xv) {
  #'
  #' count exponential mechanism density
  #' 
  #' @param v evaluation points
  #' @param eps PLB
  #' @param xv confidential value
  #' 
  const = sum(ddlap_cond_xv(0:n_char, eps, xv))
  ifelse(v >= 0 & v <= n_char, 
         ddlap_cond_xv(v, eps, xv) / const, 
         0)
}

test_thetas = seq(.001, .999, .001)
ltt = length(test_thetas)

likelihood_expmech = function(eps, xv) {
  #' 
  #' Integrated likelihood for count exponential mechanism
  #' evaluated at `test_thetas`
  #' 
  #' @param eps PLB
  #' @param xv
  #'  
  vs = sapply(
    test_thetas, 
    function(t) { 
      sum(
        dbinom(0:n_char, n_char, t) * 
          expmech_cond_xv(0:n_char, eps, xv)
      )
    }
  ) 
  vs
}

exp_cond = ggplot(data=data.frame(x=rep(test_xs, 4),
                                  Epsilon=c(rep("Eps=.01", ltxs),
                                            rep("Eps=.03", ltxs),
                                            rep("Eps=.10", ltxs),
                                            rep("Eps=.30", ltxs)),
                                  y=c(expmech_cond_xv(test_xs, .01, x_obs),
                                      expmech_cond_xv(test_xs, .03, x_obs),
                                      expmech_cond_xv(test_xs, .10, x_obs),
                                      expmech_cond_xv(test_xs, .30, x_obs)))) + 
  geom_line(aes(x=x, y=y, color=Epsilon), alpha=.7) + 
  geom_vline(xintercept=x_obs, alpha=.5, linetype="dashed") + 
  xlab("S(D)") + 
  ylab("P(S(D) | D)")

lm10 = likelihood_expmech(.1, x_obs)
lm03 = likelihood_expmech(.03, x_obs)
lm01 = likelihood_expmech(.01, x_obs)

lm10rescale = lm10/sum(lm10)
lm03rescale = lm03/sum(lm03)
lm01rescale = lm01/sum(lm01)

nonpriv_mean = x_obs / n_char
nonpriv_ci = sqrt(nonpriv_mean * (1 - nonpriv_mean) / n_char) * 2

ci_widths = data.frame(
  lb = c(
    nonpriv_mean - nonpriv_ci,
    test_thetas[which(cumsum(lm10rescale) > .025)[1]],
    test_thetas[which(cumsum(lm03rescale) > .025)[1]],
    test_thetas[which(cumsum(lm01rescale) > .025)[1]]
  ),
  ub = c(
    nonpriv_mean + nonpriv_ci,
    test_thetas[which(cumsum(lm10rescale) > .975)[1]],
    test_thetas[which(cumsum(lm03rescale) > .975)[1]],
    test_thetas[which(cumsum(lm01rescale) > .975)[1]]
  ),
  Epsilon=c("Eps=Infty",
            "Eps=.10",
            "Eps=.03",
            "Eps=.01")
)

exp_lik = ggplot(data.frame(ths=rep(test_thetas, 4),
                            lk=c(likelihood_expmech(.01, x_obs),
                                 likelihood_expmech(.03, x_obs),
                                 likelihood_expmech(.10, x_obs),
                                 dbinom(x_obs, n_char, test_thetas)),
                            Epsilon=c(rep("Eps=.01", ltt),
                                      rep("Eps=.03", ltt),
                                      rep("Eps=.10", ltt),
                                      rep("Eps=Infty", ltt)))) + 
  geom_line(aes(ths, lk, color=Epsilon)) + 
  geom_vline(data=ci_widths, 
             aes(xintercept=lb, color=Epsilon),
             linetype="dashed", alpha=.7) + 
  geom_vline(data=ci_widths, 
             aes(xintercept=ub, color=Epsilon),
             linetype="dashed", alpha=.7) + 
  xlab("Theta") + 
  ylab("L(Theta; S(D))") + 
  xlim(.05, .45)

posterior_n_i_sample = function(
  #' 
  #' Posterior samples from theta | S(D) 
  #' with prior theta ~ Beta(a, b)
  #' 
  #' @param n number of posterior samples to draw
  #' @param n_tot total number of units
  #' @param eps PLB
  #' @param dp_sample observed DP sample S(D)
  #' @param prior_a prior alpha parameter
  #' @param prior_b prior beta parameter
  n, 
  eps,
  n_tot=n_char, 
  dp_sample=x_obs, 
  prior_a=.5, 
  prior_b=.5
) {
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

post10 = posterior_n_i_sample(1e4, .1)
post03 = posterior_n_i_sample(1e4, .03)
post01 = posterior_n_i_sample(1e4, .01)
nonpriv_post = rbeta(1e4, x_obs + .5, n_char - x_obs + .5)

post_cis = data.frame(
  rbind(
    quantile(post01, c(.025, .975)),
    quantile(post03, c(.025, .975)),
    quantile(post10, c(.025, .975)),
    qbeta(c(.025, .975), x_obs + .5, n_char - x_obs + .5)
  )
)
post_cis$Epsilon = c("Eps=.01", "Eps=.03", "Eps=.10", "Eps=Infty")
names(post_cis) = c("upper", "lower", "Epsilon")

exp_posts = ggplot(data=data.frame(Epsilon=rep(c("Eps=.01", "Eps=.03", "Eps=.10", "Eps=Infty"), 
                                               each=1e4),
                                   samp=c(post01, post03, post10, nonpriv_post))) +
  geom_density(aes(x=samp, color=Epsilon)) +
  geom_vline(data=post_cis,
             mapping=aes(xintercept=upper, 
                         color=Epsilon),
             linetype="dashed", alpha=.7) + 
  geom_vline(data=post_cis,
             mapping=aes(xintercept=lower, 
                         color=Epsilon),
             linetype="dashed", alpha=.7) + 
  xlim(0, .5) + 
  xlab("Theta") + 
  ylab("p(Theta | S(D))")

exp_util_plots = grid.arrange(exp_cond, exp_lik, exp_posts, ncol=1)
ggsave("count_fp_util_lik.svg", 
       plot=exp_util_plots,
       device="svg", width=4, height=5, limitsize=FALSE)
ggsave("count_fp_util_lik.png", 
       plot=exp_util_plots,
       device="png", width=4, height=5, limitsize=FALSE)

# ----------------------------------------------------------------------------------
# Figure 5: DP posterior disclosure risks

expmech_cond_xv_vec = function(n, eps, xv) {
  #'
  #' vector-valued version of conditional exponential count distribution
  #' 
  #' @param n total number of counts
  #' @param eps PLB
  #' @param xv confidential value
  #' 
  const = sum(ddlap_cond_xv(0:n, eps, xv))
  ddlap_cond_xv(0:n, eps, xv) / const
}

expmech_post_risk = function(n, eps, xp, y, prior) {
  #'
  #' posterior disclosure risk as outlined in case study
  #' 
  #' @param n total number of counts
  #' @param eps PLB
  #' @param xp confidential value
  #' @param y DP value
  #' @param prior prior probability that sens=1 for unobserved attribute
  #'
  if (xp == 0) {
    xp0 = expmech_cond_xv_vec(n, eps, xp)
    xp1 = expmech_cond_xv_vec(n, eps, xp + 1)
    post_ratio = prior / (1 - prior) * xp1[y + 1] / xp0[y + 1]
    1 - (post_ratio / (1 + post_ratio))
  } else {
    xp0 = expmech_cond_xv_vec(n, eps, xp - 1)
    xp1 = expmech_cond_xv_vec(n, eps, xp)
    post_ratio = prior / (1 - prior) * xp1[y + 1] / xp0[y + 1]
    post_ratio / (1 + post_ratio)
  }
}

set.seed(123)
gnc_agg = df %>% 
  group_by(gender, nobility, culture_reduced) %>%
  summarise(n=n(), 
            xp=sum(alive == "DEAD"))

global_mean = mean(df$alive == "DEAD")

create_post_risk_df = function(conf_counts, grouping_name) {
  #' construct eps-DP contingency table and estimate posterior
  #' disclosure risks for each cell
  #' 
  #' @param conf_counts data.frame (aggregated counts)
  #' @param grouping_name str (for labeling)
  #' 
  conf_counts$d = 1
  sim_df = conf_counts %>%
    inner_join(expand.grid(sim_ix=1:2000, 
                           eps=c(.1, .3, 1, 2), 
                           prior=c(global_mean, .5),
                           d=1),
               by="d")
  sim_df$y = mapply(
    function(n, xp, eps) {
      sample(0:n, 1, prob=expmech_cond_xv_vec(n, eps, xp))
    },
    n=sim_df$n,
    xp=sim_df$xp,
    eps=sim_df$eps
  )
  
  sim_df$post_risk = mapply(
    expmech_post_risk,
    n=sim_df$n,
    xp=sim_df$xp,
    eps=sim_df$eps,
    y=sim_df$y,
    prior=sim_df$prior
  )
  sim_df$grouping_name = grouping_name
  sim_df
}

# run experiment and plot results
gnc_post = create_post_risk_df(gnc_agg, "Gender+Nobility+CultureReduced")
gnc_post$prior_lab = ifelse(gnc_post$prior == .5, "UniformedPrior", "ConfMeanPrior")
gnc_post$eps_lab = paste0("Eps=", sprintf("%03.1f", gnc_post$eps))

fp_risk_post = ggplot(data=gnc_post %>% 
                        # filters below for ARSIA graphics guidelines
                        filter(culture_reduced %in% top_cultures[1:5]) %>%
                        filter(eps_lab != "Eps=2.0")) + 
  geom_boxplot(aes(x=culture_reduced, y=post_risk, fill=culture_reduced), outlier.alpha=0) + 
  geom_hline(aes(yintercept=prior), linetype="dashed", alpha=.5) + 
  facet_grid(prior_lab ~ eps_lab) + 
  theme(axis.text.x=element_blank()) +
  guides(fill=guide_legend(title="Culture")) + 
  xlab("") + 
  ylab("Posterior Disclosure Risk")

ggsave("count_fp_risk_post.svg", 
       plot=fp_risk_post,
       device="svg", width=6, height=4, limitsize=FALSE)
ggsave("count_fp_risk_post.png", 
       plot=fp_risk_post,
       device="png", width=6, height=4, limitsize=FALSE)
fp_risk_post
