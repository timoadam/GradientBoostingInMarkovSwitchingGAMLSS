## Load required packages
library(MSwM) # contains the data
library(RColorBrewer)

## Load data
data(energy)
x = matrix(Oil[25:1784], ncol = 1)
y = Price[25:1784]

## Fit linear model
mod1 = FitMarkovSwitchingGAMLSS(x = x, y = y, m.stop = c(100, 200), type = "MSGLMSS", conv.tol = 1e-02, max.iter = 50)

## Fit non-linear model
mod2 = fit.msgamboostlss(x = x, y = y, m.stop = c(1600, 200), type = "MSGAMLSS", conv.tol = 1e-02, max.iter = 50)

## Get colors
pal = brewer.pal(6, "RdBu")
col1 = pal[6]
col2 = pal[1]

## fig2:3

scaleFUN = function(x) {
  sprintf("%.1f", x)
} 

## fig2:6
mod = mod1
states = apply(mod$state.probs, 2, which.max)
cols = rep(col1, 1760)
cols[states==2] = col2
minx = min(energy$Oil[25:1784])
maxx = max(energy$Oil[25:1784])
miny = min(energy$Oil[25:1784])
maxy = max(energy$Oil[25:1784])
breaks = seq(minx, maxx, length = 9)[-9]
sf = 4.25 # scaling factor
lw = 0.25 # lwd
lw2 = 0.5
muv1 = as.vector(predict(mod$mod[[1]], newdata = data.frame(x1 = breaks))$mu)
muv2 = as.vector(predict(mod$mod[[2]], newdata = data.frame(x1 = breaks))$mu)
sigmav1 = as.vector(exp(predict(mod$mod[[1]], newdata = data.frame(x1 = breaks))$sigma))
sigmav2 = as.vector(exp(predict(mod$mod[[2]], newdata = data.frame(x1 = breaks))$sigma))
m1 = as.vector(predict(mod$mod[[1]], newdata = data.frame(x1 = seq(minx, maxx, length = 1760)))$mu)
m2 = as.vector(predict(mod$mod[[2]], newdata = data.frame(x1 = seq(minx, maxx, length = 1760)))$mu)
s1 = as.vector(exp(predict(mod$mod[[1]], newdata = data.frame(x1 = seq(minx, maxx, length = 1760)))$sigma))
s2 = as.vector(exp(predict(mod$mod[[2]], newdata = data.frame(x1 = seq(minx, maxx, length = 1760)))$sigma))

d1 = data.frame(x=Oil[25:1784],y=Price[25:1784],u=seq(miny,maxy,length=1760),v1=dnorm(seq(miny,maxy,length=1760),mean=muv1[1],sd=sigmav1[1])*sf+breaks[1],v2=dnorm(seq(miny,maxy,length=1760),mean=muv1[2],sd=sigmav1[2])*sf+breaks[2],v3=dnorm(seq(miny,maxy,length=1760),mean=muv1[3],sd=sigmav1[3])*sf+breaks[3],v4=dnorm(seq(miny,maxy,length=1760),mean=muv1[4],sd=sigmav1[4])*sf+breaks[4],v5=dnorm(seq(miny,maxy,length=1760),mean=muv1[5],sd=sigmav1[5])*sf+breaks[5],v6=dnorm(seq(miny,maxy,length=1760),mean=muv1[6],sd=sigmav1[6])*sf+breaks[6],v7=dnorm(seq(miny,maxy,length=1760),mean=muv1[7],sd=sigmav1[7])*sf+breaks[7],v8=dnorm(seq(miny,maxy,length=1760),mean=muv1[8],sd=sigmav1[8])*sf+breaks[8],
w1 = dnorm(seq(miny,maxy,length=1760),mean=muv2[1],sd=sigmav2[1])*sf+breaks[1],w2=dnorm(seq(miny,maxy,length=1760),mean=muv2[2],sd=sigmav2[2])*sf+breaks[2],w3=dnorm(seq(miny,maxy,length=1760),mean=muv2[3],sd=sigmav2[3])*sf+breaks[3],w4=dnorm(seq(miny,maxy,length=1760),mean=muv2[4],sd=sigmav2[4])*sf+breaks[4],w5=dnorm(seq(miny,maxy,length=1760),mean=muv2[5],sd=sigmav2[5])*sf+breaks[5],w6=dnorm(seq(miny,maxy,length=1760),mean=muv2[6],sd=sigmav2[6])*sf+breaks[6],w7=dnorm(seq(miny,maxy,length=1760),mean=muv2[7],sd=sigmav2[7])*sf+breaks[7],w8=dnorm(seq(miny,maxy,length=1760),mean=muv2[8],sd=sigmav2[8])*sf+breaks[8])
d3 = data.frame(x=seq(minx,maxx,length=1760),y1=as.vector(predict(mod$mod[[1]],newdata=data.frame(x1=seq(minx,maxx,length=1760)))$mu),y2=as.vector(predict(mod$mod[[2]],newdata=data.frame(x1=seq(minx,maxx,length=1760)))$mu))

d4 = data.frame(x=seq(minx,maxx,length=1760),y1=qnorm(p=0.95,mean=m1,sd=s1),y2=qnorm(p=0.85,mean=m1,sd=s1),y3=qnorm(p=0.75,mean=m1,sd=s1),y4=qnorm(p=0.05,mean=m1,sd=s1),y5=qnorm(p=0.15,mean=m1,sd=s1),y6=qnorm(p=0.25,mean=m1,sd=s1),y7=qnorm(p=0.95,mean=m2,sd=s2),y8=qnorm(p=0.85,mean=m2,sd=s2),y9=qnorm(p=0.75,mean=m2,sd=s2),y10=qnorm(p=0.05,mean=m2,sd=s2),y11=qnorm(p=0.15,mean=m2,sd=s2),y12=qnorm(p=0.25,mean=m2,sd=s2))
p1 = ggplot() +
  geom_point(data = d1, aes(x = x, y = y), col = cols, shape = 21, size = 1.5, alpha = 1 / 5, fill = cols) +
  geom_polygon(data = d1, aes(y = u, x = v1), fill = col1, alpha = 1 / 5, size = 1 / 2) +
  geom_polygon(data = d1, aes(y = u, x = v2), fill = col1, alpha = 1 / 5, size = 1 / 2) +
  geom_polygon(data = d1, aes(y = u, x = v3), fill = col1, alpha = 1 / 5, size = 1 / 2) +
  geom_polygon(data = d1, aes(y = u, x = v4), fill = col1, alpha = 1 / 5, size = 1 / 2) +
  geom_polygon(data = d1, aes(y = u, x = v5), fill = col1, alpha = 1 / 5, size = 1 / 2) +
  geom_polygon(data = d1, aes(y = u, x = v6), fill = col1, alpha = 1 / 5, size = 1 / 2) +
  geom_polygon(data = d1, aes(y = u, x = v7), fill = col1, alpha = 1 / 5, size = 1 / 2) +
  geom_polygon(data = d1, aes(y = u, x = v8), fill = col1, alpha = 1 / 5, size = 1 / 2) +
  geom_path(data=d1, aes(y=u, x=v1), col=col1, size=1/2) +
  geom_path(data=d1, aes(y=u, x=v2), col=col1, size=1/2) +
  geom_path(data=d1, aes(y=u, x=v3), col=col1, size=1/2) +
  geom_path(data=d1, aes(y=u, x=v4), col=col1, size=1/2) +
  geom_path(data=d1, aes(y=u, x=v5), col=col1, size=1/2) +
  geom_path(data=d1, aes(y=u, x=v6), col=col1, size=1/2) +
  geom_path(data=d1, aes(y=u, x=v7), col=col1, size=1/2) +
  geom_path(data=d1, aes(y=u, x=v8), col=col1, size=1/2) +
  geom_polygon(data=d1,aes(y=u, x=w1),fill=col2,alpha=1/5,size=1/2) +
  geom_polygon(data=d1,aes(y=u, x=w2),fill=col2,alpha=1/5,size=1/2) +
  geom_polygon(data=d1,aes(y=u, x=w3),fill=col2,alpha=1/5,size=1/2) +
  geom_polygon(data=d1,aes(y=u, x=w4),fill=col2,alpha=1/5,size=1/2) +
  geom_polygon(data=d1,aes(y=u, x=w5),fill=col2,alpha=1/5,size=1/2) +
  geom_polygon(data=d1,aes(y=u, x=w6),fill=col2,alpha=1/5,size=1/2) +
  geom_polygon(data=d1,aes(y=u, x=w7),fill=col2,alpha=1/5,size=1/2) +
  geom_polygon(data=d1,aes(y=u, x=w8),fill=col2,alpha=1/5,size=1/2) +
  geom_path(data=d1, aes(y=u, x=w1), col=col2, size=1/2) +
  geom_path(data=d1, aes(y=u, x=w2), col=col2, size=1/2) +
  geom_path(data=d1, aes(y=u, x=w3), col=col2, size=1/2) +
  geom_path(data=d1, aes(y=u, x=w4), col=col2, size=1/2) +
  geom_path(data=d1, aes(y=u, x=w5), col=col2, size=1/2) +
  geom_path(data=d1, aes(y=u, x=w6), col=col2, size=1/2) +
  geom_path(data=d1, aes(y=u, x=w7), col=col2, size=1/2) +
  geom_path(data=d1, aes(y=u, x=w8), col=col2, size=1/2) +
  geom_line(data = d3, aes(x = x, y = y1), size = 1/2) + 
  geom_line(data = d3, aes(x = x, y = y2), size = 1/2) + 
  geom_line(data = d4, aes(x = x, y = y1), size = 1/2, linetype = "dashed") +
  geom_line(data = d4, aes(x = x,y = y4), size = 1/2, linetype = "dashed") +
  geom_line(data = d4, aes(x = x, y = y7), size = 1/2, linetype = "dashed") +
  geom_line(data = d4, aes(x = x, y = y10), size = 1/2, linetype = "dashed") +
  ggtitle("Linear model") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(name = expression(paste(OilPrice[t]," (EUR/barrel)", sep = "")), labels = scaleFUN1) +
  scale_y_continuous(name = expression(paste(EnergyPrice[t]," (c/kWh)", sep = "")), limits = c(miny, maxy))

d2 = data.frame(x=seq(1,1760,by=1),y=d1$y)
p2 = ggplot()+
  geom_path(data=d2,aes(x=x,y=y),col=cols,size=1/2,alpha=1)+
  ggtitle("Decoded time series") +
  theme(plot.title=element_text(hjust=0.5)) +
  scale_x_continuous(name=expression(t), breaks=c(1,500,1000,1500), labels=c("Feb. 4, 2002", "Jan. 2, 2004", "Dec. 4, 2005", "Nov. 4, 2007")) +
  scale_y_continuous(name=expression(paste(EnergyPrice[t]," (c/kWh)",sep="")),limits=c(miny,maxy))

## Plot Figure 6
grid.arrange(p1, p2, nrow=1)

mod = mod2
  states = apply(mod$state.probs,2,which.max)
  cols = rep(col1, 1760)
  cols[states==2] = col2
  minx = min(Oil[25:1784])
  maxx = max(Oil[25:1784])
  miny = min(Price[25:1784])
  maxy = max(Price[25:1784])
  breaks = seq(minx, maxx, length=9)[-9]
  sf = 4.25 # scaling factor
  lw = 0.25 # lwd
  lw2 = 0.5
  muv1 = as.vector(predict(mod$mod[[1]],newdata=data.frame(x1=breaks))$mu)
  muv2 = as.vector(predict(mod$mod[[2]],newdata=data.frame(x1=breaks))$mu)
  sigmav1 = as.vector(exp(predict(mod$mod[[1]],newdata=data.frame(x1=breaks))$sigma))
  sigmav2 = as.vector(exp(predict(mod$mod[[2]],newdata=data.frame(x1=breaks))$sigma))
  m1 = as.vector(predict(mod$mod[[1]],newdata=data.frame(x1=seq(minx,maxx,length=1760)))$mu)
  m2 = as.vector(predict(mod$mod[[2]],newdata=data.frame(x1=seq(minx,maxx,length=1760)))$mu)
  s1 = as.vector(exp(predict(mod$mod[[1]],newdata=data.frame(x1=seq(minx,maxx,length=1760)))$sigma))
  s2 = as.vector(exp(predict(mod$mod[[2]],newdata=data.frame(x1=seq(minx,maxx,length=1760)))$sigma))
  
  d1 = data.frame(x=Oil[25:1784],y=Price[25:1784],u=seq(miny,maxy,length=1760),v1=dnorm(seq(miny,maxy,length=1760),mean=muv1[1],sd=sigmav1[1])*sf+breaks[1],v2=dnorm(seq(miny,maxy,length=1760),mean=muv1[2],sd=sigmav1[2])*sf+breaks[2],v3=dnorm(seq(miny,maxy,length=1760),mean=muv1[3],sd=sigmav1[3])*sf+breaks[3],v4=dnorm(seq(miny,maxy,length=1760),mean=muv1[4],sd=sigmav1[4])*sf+breaks[4],v5=dnorm(seq(miny,maxy,length=1760),mean=muv1[5],sd=sigmav1[5])*sf+breaks[5],v6=dnorm(seq(miny,maxy,length=1760),mean=muv1[6],sd=sigmav1[6])*sf+breaks[6],v7=dnorm(seq(miny,maxy,length=1760),mean=muv1[7],sd=sigmav1[7])*sf+breaks[7],v8=dnorm(seq(miny,maxy,length=1760),mean=muv1[8],sd=sigmav1[8])*sf+breaks[8],
                  w1=dnorm(seq(miny,maxy,length=1760),mean=muv2[1],sd=sigmav2[1])*sf+breaks[1],w2=dnorm(seq(miny,maxy,length=1760),mean=muv2[2],sd=sigmav2[2])*sf+breaks[2],w3=dnorm(seq(miny,maxy,length=1760),mean=muv2[3],sd=sigmav2[3])*sf+breaks[3],w4=dnorm(seq(miny,maxy,length=1760),mean=muv2[4],sd=sigmav2[4])*sf+breaks[4],w5=dnorm(seq(miny,maxy,length=1760),mean=muv2[5],sd=sigmav2[5])*sf+breaks[5],w6=dnorm(seq(miny,maxy,length=1760),mean=muv2[6],sd=sigmav2[6])*sf+breaks[6],w7=dnorm(seq(miny,maxy,length=1760),mean=muv2[7],sd=sigmav2[7])*sf+breaks[7],w8=dnorm(seq(miny,maxy,length=1760),mean=muv2[8],sd=sigmav2[8])*sf+breaks[8])
  d3 = data.frame(x=seq(minx,maxx,length=1760),y1=as.vector(predict(mod$mod[[1]],newdata=data.frame(x1=seq(minx,maxx,length=1760)))$mu),y2=as.vector(predict(mod$mod[[2]],newdata=data.frame(x1=seq(minx,maxx,length=1760)))$mu))
  
  d4 = data.frame(x=seq(minx,maxx,length=1760),y1=qnorm(p=0.95,mean=m1,sd=s1),y2=qnorm(p=0.85,mean=m1,sd=s1),y3=qnorm(p=0.75,mean=m1,sd=s1),y4=qnorm(p=0.05,mean=m1,sd=s1),y5=qnorm(p=0.15,mean=m1,sd=s1),y6=qnorm(p=0.25,mean=m1,sd=s1),y7=qnorm(p=0.95,mean=m2,sd=s2),y8=qnorm(p=0.85,mean=m2,sd=s2),y9=qnorm(p=0.75,mean=m2,sd=s2),y10=qnorm(p=0.05,mean=m2,sd=s2),y11=qnorm(p=0.15,mean=m2,sd=s2),y12=qnorm(p=0.25,mean=m2,sd=s2))
  p1 = ggplot()+geom_point(data=d1,aes(x=x,y=y),col=cols,shape=21,size=1.5,alpha=1/5,fill=cols)+
    geom_polygon(data=d1,aes(y=u, x=v1),fill=col1,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=v2),fill=col1,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=v3),fill=col1,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=v4),fill=col1,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=v5),fill=col1,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=v6),fill=col1,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=v7),fill=col1,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=v8),fill=col1,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=w1),fill=col2,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=w2),fill=col2,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=w3),fill=col2,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=w4),fill=col2,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=w5),fill=col2,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=w6),fill=col2,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=w7),fill=col2,alpha=1/5,size=1/2)+
    geom_polygon(data=d1,aes(y=u, x=w8),fill=col2,alpha=1/5,size=1/2)+
    
    geom_path(data=d1, aes(y=u, x=v1), col=col1, size=1/2)+
    geom_path(data=d1, aes(y=u, x=v2), col=col1, size=1/2)+
    geom_path(data=d1, aes(y=u, x=v3), col=col1, size=1/2)+
    geom_path(data=d1, aes(y=u, x=v4), col=col1, size=1/2)+
    geom_path(data=d1, aes(y=u, x=v5), col=col1, size=1/2)+
    geom_path(data=d1, aes(y=u, x=v6), col=col1, size=1/2)+
    geom_path(data=d1, aes(y=u, x=v7), col=col1, size=1/2)+
    geom_path(data=d1, aes(y=u, x=v8), col=col1, size=1/2)+
    geom_path(data=d1, aes(y=u, x=w1), col=col2, size=1/2)+
    geom_path(data=d1, aes(y=u, x=w2), col=col2, size=1/2)+
    geom_path(data=d1, aes(y=u, x=w3), col=col2, size=1/2)+
    geom_path(data=d1, aes(y=u, x=w4), col=col2, size=1/2)+
    geom_path(data=d1, aes(y=u, x=w5), col=col2, size=1/2)+
    geom_path(data=d1, aes(y=u, x=w6), col=col2, size=1/2)+
    geom_path(data=d1, aes(y=u, x=w7), col=col2, size=1/2)+
    geom_path(data=d1, aes(y=u, x=w8), col=col2, size=1/2)+
    
    geom_line(data=d3,aes(x=x,y=y1),size=1/2)+geom_line(data=d3,aes(x=x,y=y2),size=1/2)+
    geom_line(data=d4,aes(x=x,y=y1),size=1/2,linetype="dashed")+
    geom_line(data=d4,aes(x=x,y=y4),size=1/2,linetype="dashed")+
    geom_line(data=d4,aes(x=x,y=y7),size=1/2,linetype="dashed")+
    geom_line(data=d4,aes(x=x,y=y10),size=1/2,linetype="dashed")+
    ggtitle("Non-linear model")+
    theme(plot.title=element_text(hjust=0.5))+
    scale_x_continuous(name=expression(paste(OilPrice[t]," (EUR/barrel)",sep="")),labels=scaleFUN1)+scale_y_continuous(name=expression(paste(EnergyPrice[t], " (c/kWh)",sep="")),limits=c(miny,maxy))
  
  d2 = data.frame(x=seq(1,1760,by=1),y=d1$y)
p2 = ggplot()+
  geom_path(data = d2,aes(x = x, y = y), col = cols, size = 1/2, alpha = 1)+
  ggtitle("Decoded time series") + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(name = expression(t), breaks=c(1, 500, 1000, 1500), labels=c("Feb. 4, 2002", "Jan. 2, 2004", "Dec. 4, 2005", "Nov. 4, 2007"))+
  scale_y_continuous(name = expression(paste(EnergyPrice[t], " (c/kWh)",sep="")), limits = c(miny, maxy))
  
## Plot Figure 7
p3 = grid.arrange(p1, p2, nrow = 1)


## fig2:8
## linear model
ps1 = pseudoResiduals(mod1)
q1 = qqnorm(ps1$Res)
a1 = acf(ps1$Res, na.action = na.pass, lag = 40)[[1]][,,1]
ci1 = rep(qnorm((1 + 0.95) / 2)/sqrt(acf(ps1$Res, na.action = na.pass, lag = 40)$n.used), 41)

## nonlinear model
ps2 = pseudoResiduals(mod2)
q2 = qqnorm(ps2$Res)
a2 = acf(ps2$Res, na.action=na.pass, lag = 40)[[1]][,,1]
ci2 = rep(qnorm((1 + 0.95) / 2) / sqrt(acf(ps2$Res, na.action = na.pass, lag = 40)$n.used), 41)

states = apply(mod1$state.probs, 2, which.max)
cols = rep(col1, 1760)
cols[states==2] = col2

## qq plot linear model
df1 = data.frame(x = q1$x, y = q1$y)
p1 = ggplot(data = df1) +
  geom_point(aes(x = x, y = y), col = cols, shape = 21, size = 1.5, alpha = 1 / 2, fill = cols) +
  geom_line(aes(x = seq(min(df1$x, na.rm = TRUE), max(df1$x, na.rm = TRUE), length = 1760), y = seq(min(df1$x,na.rm = TRUE), max(df1$x, na.rm = TRUE), length = 1760)), size = 1 / 2, linetype = "dashed") +
  xlab("Theoretical quantiles") +
  ylab("Sample quantiles") +
  scale_x_continuous(breaks = seq(-4, 4, by = 2), labels=c("-4.0", "-2.0", "0.0", "2.0", "4.0"), limits = c(-4, 4)) + 
  scale_y_continuous(breaks = seq(-4, 4, by = 2), labels=c("-4.0", "-2.0", "0.0", "2.0", "4.0"), limits = c(-4, 4)) +
  ggtitle("Linear model")+
  theme(plot.title = element_text(hjust = 0.5))

states = apply(mod2$state.probs,2,which.max)
cols = rep(col1,1760)
cols[states==2] = col2

## qq plot nonlinear model
df2 = data.frame(x=q2$x, y=q2$y)
p2 = ggplot(data=df2)+geom_point(aes(x=x,y=y),col=cols,shape=21,size=1.5,alpha=1/2,fill=cols)+geom_line(aes(x=seq(min(df2$x,na.rm=TRUE),max(df2$x,na.rm=TRUE),length=1760),y=seq(min(df2$x,na.rm=TRUE),max(df2$x,na.rm=TRUE),length=1760)),size=1/2,linetype="dashed")+xlab("Theoretical quantiles")+ylab("Sample quantiles")+scale_x_continuous(breaks=seq(-4,4,by=2), labels=c("-4.0","-2.0","0.0","2.0","4.0"), limits=c(-4,4))+scale_y_continuous(breaks=seq(-4,4,by=2),labels=c("-4.0","-2.0","0.0","2.0","4.0"), limits=c(-4,4))+
  ggtitle("Non-linear model")+
  theme(plot.title=element_text(hjust=0.5))

## ACF linear model
df3 = data.frame(x=0:40,y=a1,z=ci1)
p3 = ggplot(data=df3)+geom_segment(aes(x=x,y=rep(0,41),xend=x,yend=y),size=0.5)+geom_line(aes(x=x,y=z),color="black",linetype="dashed",size=0.5)+geom_line(aes(x=x,y=-z),color="black",linetype="dashed",size=0.5)+xlab("Lag")+ylab("Sample ACF")+scale_x_continuous(labels=scaleFUN0)+scale_y_continuous(labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), breaks=seq(0, 1, by=0.2))+
  ggtitle("Linear model")+
  theme(plot.title=element_text(hjust=0.5))

## ACF nonlinear model
df4 = data.frame(x=0:40,y=a2,z=ci2)
p4 = ggplot(data=df4)+geom_segment(aes(x=x,y=rep(0,41),xend=x,yend=y),size=0.5)+geom_line(aes(x=x,y=z),color="black",linetype="dashed",size=0.5)+geom_line(aes(x=x,y=-z),color="black",linetype="dashed",size=0.5)+xlab("Lag")+ylab("Sample ACF")+scale_x_continuous(labels=scaleFUN0)+scale_y_continuous(labels=c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), breaks=seq(0, 1, by=0.2))+
  ggtitle("Non-linear model")+
  theme(plot.title=element_text(hjust=0.5))

p5 = grid.arrange(p1, p2, p3, p4, nrow=2)