p.mat<- function(lamda,t, b01, b02, b13, b23, x){
  x = matrix(x, ncol = 1)
  b01 = matrix(b01, nrow = 1); b02 = matrix(b02, nrow = 1)
  b13 = matrix(b13, nrow = 1); b23 = matrix(b23, nrow = 1)
  
  lam01 = lamda[1]*exp(b01 %*% x); lam02 = lamda[2]*exp(b02 %*% x)
  lam13 = lamda[3]*exp(b13 %*% x); lam23 = lamda[4]*exp(b23 %*% x)
  
  Q.both <-   matrix(c(-( lam01 + lam02 ), lam01, lam02, 0,
                       0,-lam13, 0, lam13,
                       0,0,-lam23, lam23,
                       0,0,0,0),nrow=4,byrow=T)
  p.both <- expm(Q.both*t)
  return(p.both)
}

pseudo.df.gen = function(lamda, b01, b02, b13, b23, t0, simdata){
  # expand dataset 
  exp.df = NULL
  id.fake = c(0, 0, 0)
  for(i in 1:n){
    #print(i)
    indv.i = simdata[simdata$id==i,][1,]
    wi = rep(0, 3)
    pi = rep(0, 3)
    if(indv.i$miss.i == 0){
      wi[indv.i$label] = 1
      id.fake = max(id.fake) + c(1, 2, 3)
      interm = rbind(indv.i, indv.i, indv.i)
      interm$id.fake = id.fake
      interm$wi = wi
      interm$label.fake = c(1, 2, 3)
      
      exp.df = rbind(exp.df, interm)
    }else{
      tL = indv.i$timeL; tU = indv.i$timeU
      sL = indv.i$stageL; sU = indv.i$stageU
      
      prob_tLt0 = p.mat(lamda = lamda, b01 = b01, b02 = b02,
                        b13 = b13, b23 = b23, x = c(indv.i$x1, indv.i$x2), t = (t0 - tL))
      
      if(is.infinite(tU)){
        prob_t0tU =  prob_tLtU = matrix(0, nrow = 4, ncol = 4)
        prob_t0tU[,4] = prob_tLtU[,4] = 1
      }else{
        prob_t0tU = p.mat(lamda = lamda, b01 = b01, b02 = b02,
                          b13 = b13, b23 = b23, x = c(indv.i$x1, indv.i$x2), t = (tU - t0))
        prob_tLtU = p.mat(lamda = lamda, b01 = b01, b02 = b02,
                          b13 = b13, b23 = b23, x = c(indv.i$x1, indv.i$x2), t = (tU - tL))
      }
      
      if(sL == 2 || sL == 3){ 
        py1 = 0; wi[1] = py1
        py2 = prob_tLt0[sL, sL]* prob_t0tU[sL, sU]/prob_tLtU[sL, sU]; wi[2] = py2
        py3 = prob_tLt0[sL, sU]* prob_t0tU[sU, sU]/prob_tLtU[sL, sU]; wi[3] = py3
        
      }
      if(sL == 1){ 
        if(sU == 4){
          py1 = prob_tLt0[sL, 1] * prob_t0tU[1, sU]/prob_tLtU[sL, sU]; wi[1] = py1
          p12 = prob_tLt0[sL, 2] * prob_t0tU[2, sU]/prob_tLtU[sL, sU]
          p13 = prob_tLt0[sL, 3] * prob_t0tU[3, sU]/prob_tLtU[sL, sU]
          py2 = p12 + p13; wi[2] = py2
          py3 = prob_tLt0[sL, 4] * prob_t0tU[4, sU]/prob_tLtU[sL, sU]; wi[3] = py3
        }
        if(sU == 2 || sU == 3){
          py1 = prob_tLt0[sL, 1] * prob_t0tU[1, sU]/prob_tLtU[sL, sU]; wi[1] = py1
          py2 = prob_tLt0[sL, sU] * prob_t0tU[sU, sU]/prob_tLtU[sL, sU]; wi[2] = py2
          py3 = 0; wi[3] = py3
        }
      }
      
      id.fake = max(id.fake) + c(1, 2, 3)
      interm = rbind(indv.i, indv.i, indv.i)
      interm$id.fake = id.fake
      interm$wi = wi
      interm$label.fake = c(1, 2, 3)
      exp.df = rbind(exp.df, interm)
    }
  }
  
  return(exp.df)
}

PDI1.fun = function(i, exp.df, perc_miss, prob.mat, line1.prob){
  
  pp = prob.mat
  
  id.rm = i
  df.c = exp.df[-which(exp.df$id==id.rm),]
  
  weight.i = exp.df[which(exp.df$id==id.rm),]$wi[1]

  idx = i
  idx.pp2 = df.c$id.fake[df.c$label.fake==2] 
  idx.pp3 = df.c$id.fake[df.c$label.fake==3] 
  
  PDI.1 = sum((line1.prob[idx, 1] > pp[idx.pp2, 1])*sqrt(weight.i)*df.c$wi[df.c$label.fake==2]) * sum((line1.prob[idx, 1] > pp[idx.pp3, 1])*sqrt(weight.i)*df.c$wi[df.c$label.fake==3])
  
  return(PDI.1)
}

PDI2.fun = function(i, exp.df, perc_miss, prob.mat, line1.prob){
  pp = prob.mat
  id.rm = i
  df.c = exp.df[-which(exp.df$id==id.rm),]
  
  weight.i = exp.df[which(exp.df$id==id.rm),]$wi[2]

  idx = i
  idx.pp1 = df.c$id.fake[df.c$label.fake==1] 
  idx.pp3 = df.c$id.fake[df.c$label.fake==3] 
    
  PDI.2 = sum((line1.prob[idx, 2] > pp[idx.pp1, 2])*sqrt(weight.i)*df.c$wi[df.c$label.fake==1]) * sum((line1.prob[idx, 2] > pp[idx.pp3, 2])*sqrt(weight.i)*df.c$wi[df.c$label.fake==3])
  
  return(PDI.2)
}

PDI3.fun = function(i, exp.df, perc_miss, prob.mat, line1.prob){
  
  pp = prob.mat
  id.rm = i
  df.c = exp.df[-which(exp.df$id==id.rm),]
  
  weight.i = exp.df[which(exp.df$id==id.rm),]$wi[3]

  idx = i
  idx.pp2 = df.c$id.fake[df.c$label.fake==2] 
  idx.pp1 = df.c$id.fake[df.c$label.fake==1] 

  PDI.3 = sum((line1.prob[idx, 3] > pp[idx.pp1, 3])*sqrt(weight.i)*df.c$wi[df.c$label.fake==1]) * sum((line1.prob[idx, 3] > pp[idx.pp2, 3])*sqrt(weight.i)*df.c$wi[df.c$label.fake==2])
  
  return(PDI.3)
}