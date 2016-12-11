
## code from Ryan # species abundance changes

load("/Users/nkmstar/Dropbox/secplot_recovery/abundFitBT_sec_23Nov2016.rdata")
load("/Users/KangMin/Dropbox/secplot_recovery/abundFitBT_sec_23Nov2016.rdata")
library(nlme)
dbh_threshold = "1cm"

filenames=("/Users/nkmstar/Dropbox/secplot_recovery/abundFitBT_sec_23Nov2016.rdata")
filenames=("/Users/KangMin/Dropbox/secplot_recovery/abundFitBT_sec_23Nov2016.rdata")

CTFSabundmodel_new = read.table("/Users/nkmstar/Dropbox/secplot_recovery/CTFS_new_fits_BT.txt", header=T)
CTFSabundmodel_new = read.table("/Users/KangMin/Dropbox/secplot_recovery/CTFS_new_fits_BT.txt", header=T)

CTFSdemog_list = lapply(filenames, function(filename) {
  cat(filename,'\n')
  temp = load(filename)
  eval(parse(text=paste('x =',temp[1])))
  demog_all = list()
  if ( dbh_threshold == "1cm" ) y = x[[1]] else y = x[[2]]
  # exclude the last element, which should be the overall change from first to last census
  for ( i in 1:max(1,length(y)) )
  {
    #    cat('census interval ',i,'\n',sep='')
    #    print(table(substr(fromjulian(y[[i]]$demog$date1),7,10)))
    #    print(table(substr(fromjulian(y[[i]]$demog$date2),7,10)))
    # fix for Mudumalai
    if ( 'lowermean.1' %in% names(y[[i]]$demog) )
    {
      y[[i]]$demog$lowermeanR = y[[i]]$demog$lowermean.1
      y[[i]]$demog$uppermeanR = y[[i]]$demog$uppermean.1
    }
    S = dim(y[[i]]$demog)[1]    
    site_name = substr(filename,49,nchar(filename)-15) ###### change this accordingly #################
    demog = data.frame(N1=y[[i]]$demog$N1,N2=y[[i]]$demog$N2,S=y[[i]]$demog$S,
                       time=y[[i]]$demog$time)
    demog_all[[i]] = demog
    names(demog_all)[i] = paste(site_name,substr(names(y)[i],7,9),sep='')
  }
  return(demog_all)
})

CTFSdemog = unlist(CTFSdemog_list,recursive=F)

## rasymexp function
phi_1 = 10
phi_2 = 7
Gamma = -0.07
phi = phi_1*phi_2/(phi_1+phi_2)

rasymexp <- function( n, Gamma, phi_1, phi_2 )
{
  return (qasymexp(runif(n), Gamma, phi_1, phi_2 ))
}

dasymexp <- function(x, Gamma, phi_1, phi_2)
{
  phi = phi_1*phi_2/(phi_1+phi_2)
  is_lo = (x < Gamma)
  rval = numeric(length(x))
  rval_lo = phi*exp(-phi_2*(Gamma-x))
  rval_hi = phi*exp(-phi_1*(x-Gamma))
  rval[is_lo] = rval_lo[is_lo]
  rval[!is_lo] = rval_hi[!is_lo]
  rval 
}

qasymexp <- function(y, Gamma, phi_1, phi_2 )
{
  phi = phi_1*phi_2/(phi_1+phi_2)
  is_lo = (y < 1-phi/phi_1)
  rval = numeric(length(y))
  rval_lo = Gamma+(1/phi_2)*log(phi_2*y/phi)
  rval_hi = Gamma-(1/phi_1)*log(phi_1*(1-y)/phi)
  rval[is_lo] = rval_lo[is_lo]
  rval[!is_lo] = rval_hi[!is_lo]
  rval
}

# Rick's dasymexp (should give same result as dasymexp)
dasymexp_Rick=function(x,center,rate1,rate2,log=FALSE)
{
  logy=numeric()
  right=x>=center
  left=x<center
  k=rate1*rate2/(rate1+rate2)
  
  logy[right]=log(k)-rate1*(x[right]-center)
  logy[left]=log(k)-rate2*(center-x[left])
  
  if(log) return(logy)
  return(exp(logy))
}



# code for plotting most recent census does not draw the environmental variance line & CI
# refer to printed codes for more details

# code for longest census interval

rho_model_type = 'asymexp'

# these are the parameters of the asymmetric power law distribution
c = CTFSabundmodel_new$hyperR
r1 = 1/CTFSabundmodel_new$hyperSDlow
r2 = 1/CTFSabundmodel_new$hyperSDup
CTFSabundmodel_new$c = c
CTFSabundmodel_new$r1 = r1
CTFSabundmodel_new$r2 = r2

if ( rho_model_type == 'asympower' )
{
  E_R = c+(r2-r1)*(3+r1+r2)/((2+r1)*(2+r2)*(2+r1+r2))
  E_R_sq = ((1+r1)*(1+r2)/(2+r1+r2))*((2+c*(3+r2)*(-2+c*(2+r2)))/((1+r2)*(2+r2)*(3+r2))+(2+c*(3+r1)*(2+c*(2+r1)))/((1+r1)*(2+r1)*(3+r1)))
  CTFSabundmodel_new$env_var = E_R_sq - E_R^2
} else if ( rho_model_type == 'asymexp' )
{
  E_R = c+1/r2-1/r1
  var_R = 1/r1^2 + 1/r2^2
  CTFSabundmodel_new$env_var = var_R
}

#verify_rhos(CTFSdemog,CTFSabundmodel_new,c,r1,r2,E_R,var_R+E_R^2)

E_mu = exp(CTFSabundmodel_new$hyperMu+CTFSabundmodel_new$hyperSD^2/2)
var_mu = (exp(CTFSabundmodel_new$hyperSD^2)-1)*exp(2*CTFSabundmodel_new$hyperMu+CTFSabundmodel_new$hyperSD^2)
CTFSabundmodel_new$lognorm_sigma = CTFSabundmodel_new$hyperSD
CTFSabundmodel_new$lognorm_mu = CTFSabundmodel_new$hyperMu
CTFSabundmodel_new$dem_var = E_R + 2*E_mu# - var_mu - E_mu^2
CTFSabundmodel_new$E_mu = E_mu

#verify_mus(CTFSdemog,CTFSabundmodel_new,CTFSabundmodel_new$lognorm_mu,CTFSabundmodel_new$lognorm_sigma,E_mu,var_mu+E_mu^2)

vars_from_mu_rho <- function(dT,mu_mu,sigma_mu,c,r1,r2,rho_model_type,N1_range)
{
  n_rep = 100000
  mu = rlnorm(n_rep,mu_mu,sigma_mu)
  
  # note that r2 and r1 are switched deliberately in the function calls below
  if ( rho_model_type == 'asympower' )
    rho = rasympower(n_rep,c,r2,r1)
  else if ( rho_model_type == 'asymexp' )
    rho = rasymexp(n_rep,c,r2,r1)
  else
    stopifnot(F)
  
  delta = 1-exp(-dT*mu)
  lambda = exp(dT*rho)
  
  E_delta = mean(delta); var_delta = var(delta)
  E_lambda = mean(lambda); var_lambda = var(lambda)
  
  env_var = var_lambda
  dem_var = E_lambda + 2*E_delta - var_delta - E_delta^2 - 1
  
  E_log_dN_sq = numeric(length(N1_range))
  q_hi_log_dN_sq = numeric(length(N1_range))
  q_lo_log_dN_sq = numeric(length(N1_range))
  beta = pmax(0,lambda+delta-1)
  for ( i in 1:length(N1_range) )
  {
    N1 = N1_range[i]
    N2 = rbinom(n_rep,N1,1-delta) + rpois(n_rep,N1*beta) 
    dN = N2-N1
    log_dN_sq = log(dN[dN!=0]^2) 
    E_log_dN_sq[i] = mean(log_dN_sq)
    #sd_log_dN_sq = sd(log_dN_sq)
    #q_lo_log_dN_sq[i] = E_log_dN_sq[i]-1.96*sd_log_dN_sq
    #q_hi_log_dN_sq[i] = E_log_dN_sq[i]+1.96*sd_log_dN_sq
    q_lo_log_dN_sq[i] = quantile(log_dN_sq,0.025)
    q_hi_log_dN_sq[i] = quantile(log_dN_sq,0.975)
  }
  
  return(list(env_var=env_var,dem_var=dem_var,
              E_log_dN_sq=E_log_dN_sq,q_lo_log_dN_sq=q_lo_log_dN_sq,q_hi_log_dN_sq=q_hi_log_dN_sq))
}


# graph plotting 

do_graph_censuses = c('BT_sec_1.2')
do_graph_titles = c('BT sec 2004-2008')
# do_graph_censuses = c('BT_sec_2.3')
# do_graph_titles = c('BT sec 2008-2012')
do_graph_sites = c('Bukit Timah')
filename_base = 'figure_single_interval'

do_graph_censuses = c('BT_sec_1.3')
do_graph_titles = c('BT sec 2004-2012')
do_graph_sites = c('Bukit Timah')
filename_base = 'figure_longest_interval'

n_graph = length(do_graph_censuses)

if ( filename_base == 'figure_single_censuses' | filename_base == 'figure_longest_interval' )
{
  site_lat_names_index_0 = c('Bukit Timah')
  site_lats_index_0 = c(1.20)
  
  iii = order(site_lats_index_0)
  site_lat_names_index = site_lat_names_index_0[iii]
  site_lats_index = site_lats_index_0[iii]
  jjj = match(site_lat_names_index,do_graph_sites)
  jjj = jjj[!is.na(jjj)]
} else
{
  jjj = 1:n_graph
}

# maximum abundance in any census at any plot
max_abund = max(sapply(CTFSdemog,function(x) {max(x$N1)}))
max_dN_sq = max(sapply(CTFSdemog,function(x) {max((x$N2-x$N1)^2)}))

graphics.off()
dev.new(width=12,height=10)
#x11(width=12,height=10,type='cairo')
n_graph_per_row = 1

par(mfrow=c(ceiling(n_graph/n_graph_per_row),n_graph_per_row),mar=c(3,2,2,2),oma=c(5,7,0,0))

all_exponents = numeric(n_graph)

all_x1 = c()

do_synch_env = F

for ( i in 1:n_graph )
{
  data_name = do_graph_censuses[jjj[i]]
  i2 = match(data_name,names(CTFSdemog))
  demog = CTFSdemog[[i2]]
  
  temp = CTFSabundmodel_new[data_name,]
  dT = mean(demog$time)
  N1s_for_plot = round(exp(seq(log(1),log(max(demog$N1)),length.out=100)))
  vvv = vars_from_mu_rho(dT,temp$lognorm_mu,temp$lognorm_sigma,temp$c,temp$r1,temp$r2,
                         rho_model_type,N1s_for_plot)
  env_var = vvv$env_var
  dem_var = vvv$dem_var
  E_log_dN_sq = vvv$E_log_dN_sq
  q_lo_log_dN_sq = vvv$q_lo_log_dN_sq
  q_hi_log_dN_sq = vvv$q_hi_log_dN_sq
  
  S = dim(demog)[1]
  J1 = sum(demog$N1)
  J2 = sum(demog$N2)
  Z_tot = sum(demog$S)
  zeta = Z_tot/J1
  beta = (J2-Z_tot)/J1
  
  E_N2 = demog$N1*(zeta+beta)
  
  #demog$N2 = rbinom(S,demog$N1,zeta) + rpois(S,demog$N1*beta)
  
  J1 = sum(demog$N1)
  J2 = sum(demog$N2)
  
  x1 = demog$N1#/J1
  x2 = demog$N2#/J2
  y = x2-x1
  
  # neutral model
  n_point = 20
  x1s = exp(seq(log(min(x1)),log(max(x1)),length.out=n_point))
  #N1s = min(demog$N1):max(demog$N1)
  n_point = length(x1s)
  y_mean = numeric(n_point)
  y_mean_log_scale = numeric(n_point)
  y_lo = numeric(n_point)
  y_hi = numeric(n_point)
  for ( j in 1:n_point )
  {
    n_rep = 1000
    if ( j == 1 | j == n_point )
      x1s_j = rep(x1s[j],n_rep) else
        x1s_j = runif(n_rep,x1s[j-1],x1s[j+1])
      N1s_j = round(x1s_j)
      if ( do_synch_env )
      {
        # allow the community size to change
        beta2 = beta
        zeta2 = zeta
      }
      else
      {
        # force the community size to be constant
        beta2 = (beta + (1-zeta))/2
        zeta2 = 1-beta2
      }
      N2_pred = rbinom(n_rep,N1s_j,zeta2) + rpois(n_rep,N1s_j*beta2)
      #y_pred = (N2_pred-(beta+zeta)*N1s_j)^2
      y_pred = (N2_pred-N1s_j)^2
      y_mean[j] = mean(y_pred[y_pred>0])
      y_mean_log_scale[j] = exp(mean(log(y_pred[y_pred>0])))
      y_lo[j] = max(1e-10,quantile(y_pred[y_pred>0],0.025))
      y_hi[j] = max(1e-10,quantile(y_pred[y_pred>0],0.975))
  }
  
  all_x1 = c(all_x1,x1)
  
  # draw data points
  par(mar=c(0,0,0,0))
  x1y <- cbind(x1,y)
  neg.gro <- x1y[sign(y)==-1,]
  zero.gro <- x1y[sign(y)==0,]
  pos.gro <- x1y[sign(y)==1,]
  plot(neg.gro[,1],(neg.gro[,2])^2,log='xy',xlab='',ylab='', col="black", pch=25, bg="black", cex.axis=1.5, xlim=c(1,max_abund), ylim=c(1,max_dN_sq), yaxt='n', xaxt='n')
  points(pos.gro[,1],(pos.gro[,2])^2,log='xy',xlab='',ylab='', col="black", pch=24, bg="black", cex.axis=1.5, xlim=c(1,max_abund), ylim=c(1,max_dN_sq), yaxt='n', xaxt='n')
  
  # x axes
  if ( i>n_graph-n_graph_per_row )
    axis(1,at=c(10^0,10^1,10^2,10^3,10^4),
         labels=c(1,10,100,1000,10000),cex.axis=1.5)
  else
    axis(1,at=c(10^0,10^1,10^2,10^3,10^4),labels=c('','','','',''))
  
  # y axes
  if ( n_graph_per_row==1 )
    axis(2,at=c(10^0,10^2,10^4,10^6,10^8),
         labels=c(1,expression(10^2),expression(10^4),expression(10^6),expression(10^8)),cex.axis=1.5)
  else
    axis(2,at=c(10^0,10^2,10^4,10^6,10^8),labels=c('','','','',''))
  #text(sqrt(max_abund),max_dN_sq^0.95,do_graph_titles[jjj[i]],cex=2)
  
  # draw the fitted model
  #lines(xs,env_var*(xs^2)+dem_var*(xs),col='black',lwd=2)
  lines(N1s_for_plot,exp(E_log_dN_sq),col='black',lwd=2)
  #lines(xs,temp$env_var*(dT^2)*(xs^2)+temp$dem_var*dT*(xs),col='black',lwd=2,lty=2)
  polygon(c(N1s_for_plot,rev(N1s_for_plot),N1s_for_plot[1]),
          c(exp(q_lo_log_dN_sq),rev(exp(q_hi_log_dN_sq)),exp(q_lo_log_dN_sq)[1]),
          # col=rgb(0,0,0,0.15),border=rgb(0,0,0,0.15))
          col=rgb(220,220,220,100, maxColorValue=255),border=rgb(220,220,220,100, maxColorValue=255))
  
  # draw neutral prediction
  xs = exp(seq(log(min(x1)),log(max(x1)),length.out=100))
  #lines(xs,(zeta*(1-zeta)+beta)*xs,col='blue',lwd=2,lty=2)
  #lines(x1s,y_mean,col='blue',lwd=2,lty=2)
  lines(x1s,y_mean_log_scale,col='black',lwd=2,lty=2)
  polygon(c(x1s,rev(x1s),x1s[1]),c(y_lo,rev(y_hi),y_lo[1]),
          # col=rgb(220,220,220,100, maxColorValue=255),border=rgb(220,220,220,100, maxColorValue=255))
          col=rgb(0,0,0,0.15),border=rgb(0,0,0,0.15))
  
  x_min = min(x1)
  y_min = (zeta*(1-zeta)+beta)*x_min
  
  # draw line of best fit
  y_ok = y[y>0]
  x1_ok = x1[y>0]
  mygls = gls(log(y_ok)~log(x1_ok))
  #lines(xs,exp(mygls$coef[1]+mygls$coef[2]*log(xs)),col='red',lwd=2,lty=1)
  #text(max(x1)/100,max(y)/2,paste('exponent = ',round(mygls$coef[2],2),sep=''),col='red')
  
  # draw environmental variance prediction
  mygnls = gnls(log(y_ok)~log_a+2*log(x1_ok),
                start=list(log_a=-10))#as.numeric(mygls$coef[1])))
  # lines(xs,exp(mygnls$coef[1])*(xs^2),col=rgb(0.1,0.9,0.1),lwd=2,lty=2)
  # lines(xs,exp(mygnls$coef[1])*(xs^2),col="black",lwd=2,lty=3)
  
  if ( i == 1 )
  {
    if ( do_synch_env )
    {
      legend_text = c('data','synch env.','env.','fitted')
    }
    else
    {
      # legend_text = c('data','neutral','env.','fitted')
      legend_text = c('positive change','negative change','neutral','fitted')
    }
    
    legend(min(x1),max_dN_sq^0.99,legend=legend_text,cex=1.5,
           col=c('black',col='black'), bty="n",
           lty=c(0,0,2,1),
           lwd=c(NA,NA,2,2),
           pch=c(24,25,NA,NA), 
           text.width=1)
    
    mtext(expression(paste('Initial abundance ',(italic(N[1])))),1,outer=T,cex=1.5,padj=1.6)
    mtext(expression(paste('Squared abundance change ',(italic((N[2]-N[1])^2)))),2,outer=T,cex=1.5,padj=-1.5)
  }
  
  all_exponents[i] = mygls$coef[2]
  
}








# larger trees 
dbh_threshold = "10cm"
CTFSabundmodel_new = read.table("/Users/nkmstar/Dropbox/secplot_recovery/CTFS_new_fits_BT.10.txt", header=T)
CTFSabundmodel_new = read.table("/Users/KangMin/Dropbox/secplot_recovery/CTFS_new_fits_BT.10.txt", header=T)







