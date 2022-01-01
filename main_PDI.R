library(msm); library(readr)

source("utility_fun.R")

t0 = 1
n = 1500
n.train = 1000
n.test = 500

mat.q = rbind(c(-1, 1, 1, 0),
              c(0, -1, 0, 1),
              c(0, 0, -1, 1),
              c(0, 0, 0, 0))

simdata = read_rds("simdata")
simdata.train = read_rds("simdata_train")

#MSM for weights and parameter estimation 
fitd <- msm(state ~ obstime, subject=id, data=simdata.train, qmatrix=mat.q, covariates = ~x1 + x2,
            control=list(fnscale=100000, reltol=1e-4, maxit=50000), center = F)

param.est = fitd$estimates.t
lam01_h = param.est[1]; lam02_h = param.est[2]; lam13_h = param.est[3];lam23_h = param.est[4]
b01_h = c(param.est[5], param.est[9]);  b02_h = c(param.est[6], param.est[10])
b13_h = c(param.est[7], param.est[11]); b23_h = c(param.est[8], param.est[12])

# expanding into the pseudo-dataset
exp.df = pseudo.df.gen(lamda = c(lam01_h, lam02_h, lam13_h, lam23_h), b01 = b01_h,
                       b02 = b02_h, b13 = b13_h, b23 = b23_h, t0 = t0, simdata = simdata)

exp.df.train = exp.df[exp.df$id %in% unique(simdata.train$id),]
exp.df.test = exp.df[!exp.df$id %in% unique(simdata.train$id),]
exp.df.test$id = rep(1:(length(unique(exp.df.test$id))), each = 3)
exp.df.test$id.fake = c(1:nrow(exp.df.test))

prob.mat = NULL
for(i in 1:nrow(exp.df.test)){
  indv.i = exp.df.test[i,]
  est.p = p.mat(lamda = c(lam01_h, lam02_h, lam13_h, lam23_h), b01 = b01_h,
                b02 = b02_h, b13 = b13_h, b23 = b23_h, x = c(indv.i$x1, indv.i$x2), t = t0)
  
  prob.mat = rbind(prob.mat, c(est.p[1,1], (est.p[1,2] + est.p[1,3]), est.p[1,4]))
}

line1.prob = prob.mat[seq(1, nrow(prob.mat), 3),] # every 3rd row 

out1 = lapply(exp.df.test$id[exp.df.test$label.fake==1], PDI1.fun, exp.df = exp.df.test, prob.mat = prob.mat,
              line1.prob = line1.prob, perc_miss = 999)
out2 = lapply(exp.df.test$id[exp.df.test$label.fake==2], PDI2.fun, exp.df = exp.df.test, prob.mat = prob.mat,
                line1.prob = line1.prob, perc_miss = 999)
out3 = lapply(exp.df.test$id[exp.df.test$label.fake==3], PDI3.fun, exp.df = exp.df.test, prob.mat = prob.mat,
                line1.prob = line1.prob, perc_miss = 999)

PDI.1 = apply(do.call(rbind, out1), 2, sum)/(sum(exp.df.test$wi[exp.df.test$label.fake == 1])*sum(exp.df.test$wi[exp.df.test$label.fake == 2])*sum(exp.df.test$wi[exp.df.test$label.fake == 3]))
PDI.2 = apply(do.call(rbind, out2), 2, sum)/(sum(exp.df.test$wi[exp.df.test$label.fake == 1])*sum(exp.df.test$wi[exp.df.test$label.fake == 2])*sum(exp.df.test$wi[exp.df.test$label.fake == 3]))
PDI.3 = apply(do.call(rbind, out3), 2, sum)/(sum(exp.df.test$wi[exp.df.test$label.fake == 1])*sum(exp.df.test$wi[exp.df.test$label.fake == 2])*sum(exp.df.test$wi[exp.df.test$label.fake == 3]))


# print out all PDIs
print(c(PDI.1, PDI.2, PDI.3))

# overall PDI
print(mean(c(PDI.1, PDI.2, PDI.3)))


