### symmetric unimodal
quaddecons <- function(y,h,support=NULL,d=NULL,plot=TRUE,lam=NULL,mode0=TRUE){
  n = length(y)
  me=0
  if(mode0==FALSE){me=median(y);y=y-me}
  if(is.null(support)){
    support=c(0,0)
    dis=quantile(y,0.99)-quantile(y,0.01)
    support[1]=quantile(y,0.01)-dis/3
    support[2]=quantile(y,0.99)+dis/3
    if(-support[1]>support[2]){support[2]=-support[1]}else{support[1]=-support[2]}
  }
  if(!is.null(d)) {
    l = h/d
    if(d>h|round(h/d)!=h/d) d=NULL
  }
  if(is.null(d)){
    l = round((22*n^(1/9))*h/diff(support))
    d = h/l 
  }
  if(l<2) {
    l = 2
    d = h/l
  }
  support[1] = floor(support[1]/d)*d
  support[2] = ceiling(support[2]/d)*d
  tk = seq(support[1],support[2],d)
  m = length(tk)-3
  yp = seq(support[1]-h,support[2]+h,length.out = 2000); np = length(yp)
  delta  = matrix(0,nrow=np,ncol=m)
  for(i in 1:m){
    delta[yp>tk[i]&yp<=tk[i+1],i] = 2*(yp[yp>tk[i]&yp<=tk[i+1]]-tk[i])^2/d^2/3
    delta[yp>=tk[i+1]&yp<=tk[i+2],i] = 1-4*(yp[yp>=tk[i+1]&yp<=tk[i+2]]-tk[i+1]-d/2)^2/d^2/3
    delta[yp>tk[i+2]&yp<=tk[i+3],i] = 2*(yp[yp>tk[i+2]&yp<=tk[i+3]]-tk[i+3])^2/d^2/3
  }
  gamma = matrix(0,nrow=np,ncol=m)
  for(j in 1:m){
    gamma[yp>=tk[j]-h&yp<tk[j+1]-h,j] = (yp[yp>=tk[j]-h&yp<tk[j+1]-h]+h-tk[j])^3/9/h/d^2
    gamma[yp>=tk[j+1]-h&yp<tk[j+2]-h,j] = -2*(yp[yp>=tk[j+1]-h&yp<tk[j+2]-h]+h-tk[j+1]-d/2)^3/9/h/d^2+(yp[yp>=tk[j+1]-h&yp<tk[j+2]-h]+h-tk[j+1])/2/h+d/h/12
    gamma[yp>=tk[j+2]-h&yp<tk[j+3]-h,j] = (yp[yp>=tk[j+2]-h&yp<tk[j+3]-h]+h-tk[j+3])^3/9/h/d^2+2*d/3/h
    gamma[yp>=tk[j+3]-h&yp<=tk[j]+h,j] = 2*d/3/h
    gamma[yp>tk[j]+h&yp<=tk[j+1]+h,j] = -(yp[yp>tk[j]+h&yp<=tk[j+1]+h]-h-tk[j])^3/9/h/d^2+2*d/3/h
    gamma[yp>tk[j+1]+h&yp<=tk[j+2]+h,j] = 2*(yp[yp>tk[j+1]+h&yp<=tk[j+2]+h]-h-tk[j+1]-d/2)^3/9/h/d^2+(-yp[yp>tk[j+1]+h&yp<=tk[j+2]+h]+tk[j+2]+h)/2/h+d/12/h
    gamma[yp>tk[j+2]+h&yp<=tk[j+3]+h,j] = (-yp[yp>tk[j+2]+h&yp<=tk[j+3]+h]+tk[j+3]+h)^3/9/h/d^2
  }
  if(is.null(lam)){
    gammay = matrix(0,nrow=n,ncol=m)
    for(j in 1:m){
      gammay[y>=tk[j]-h&y<tk[j+1]-h,j] = (y[y>=tk[j]-h&y<tk[j+1]-h]+h-tk[j])^3/9/h/d^2
      gammay[y>=tk[j+1]-h&y<tk[j+2]-h,j] = -2*(y[y>=tk[j+1]-h&y<tk[j+2]-h]+h-tk[j+1]-d/2)^3/9/h/d^2+(y[y>=tk[j+1]-h&y<tk[j+2]-h]+h-tk[j+1])/2/h+d/h/12
      gammay[y>=tk[j+2]-h&y<tk[j+3]-h,j] = (y[y>=tk[j+2]-h&y<tk[j+3]-h]+h-tk[j+3])^3/9/h/d^2+2*d/3/h
      gammay[y>=tk[j+3]-h&y<=tk[j]+h,j] = 2*d/3/h
      gammay[y>tk[j]+h&y<=tk[j+1]+h,j] = -(y[y>tk[j]+h&y<=tk[j+1]+h]-h-tk[j])^3/9/h/d^2+2*d/3/h
      gammay[y>tk[j+1]+h&y<=tk[j+2]+h,j] = 2*(y[y>tk[j+1]+h&y<=tk[j+2]+h]-h-tk[j+1]-d/2)^3/9/h/d^2+(-y[y>tk[j+1]+h&y<=tk[j+2]+h]+tk[j+2]+h)/2/h+d/12/h
      gammay[y>tk[j+2]+h&y<=tk[j+3]+h,j] = (-y[y>tk[j+2]+h&y<=tk[j+3]+h]+tk[j+3]+h)^3/9/h/d^2
    }
  }
  hmat = matrix(0,nrow=m,ncol=m)
  for(i in 1:m){
    for(j in 1:i)
      hmat[i,j] = hmat[j,i]=sum(gamma[,i]*gamma[,j])*(yp[2]-yp[1])
  }
  avec = 1:m*0+4*d/3
  nz=floor(length(tk)/2)+1
  smat=matrix(0,nrow=m/2,ncol=m)
  for(i in 1:(m/2)){
    smat[i,i]=-1;smat[i,m-i+1]=1
  }
  wmat = matrix(0,nrow=m,ncol=m/2-1)
  ww=cbind(t(smat),avec)
  pw=-ww%*%solve(t(ww)%*%ww)%*%t(ww)
  for(i in 1:m){pw[i,i]=pw[i,i]+1}
  for(i in 1:(m/2-1)){
    xx=rnorm(m)
    wmat[,i]=pw%*%xx
  }
  wmat = matrix(0,nrow=m,ncol=m/2-1)
  for(i in 1:(m/2-1)){
    wmat[i,i]=wmat[m-i+1,i]=1
  }
  wmat[m/2,]=wmat[m/2+1,]=-1
  b0 = 1:m*0+3/m/d/4
  
  cvec = makec(y,tk,d,h)
  D <- matrix(0,m+1,m)
  for(i in 1:m){
    D[i,i]=3
    D[i+1,i]=-3
  }
  for (i in 1:(m-1)) D[i,i+1]=-1
  for(i in 1:(m-1)) D[i+2,i]=1
  D = 4*D/3*sqrt(d)
  amat = matrix(0,m/2,m)
  amat[1,1] = 1
  for(i in 2:(m/2)){
    amat[i,i-1] = -1
    amat[i,i] = 1
  }
  amat=rbind(c(rep(0,m/2-1),1,-1,rep(0,m/2-1)),amat)
  
  c0vec=t(wmat)%*%(cvec-hmat%*%b0)
  q0mat=t(wmat)%*%(hmat)%*%wmat
  nv <- ceiling(n/10)
  fold <- rep(1:10,nv)[1:n]
  ##if(nv<n/10) {fold <- c(fold,1:(n-10*nv))}
  nf <- 1:10
  for(i in 1:10){nf[i] <- sum(fold==i)}
  fold <- sample(fold,n)
  risk=NULL
  id=NULL
  if(is.null(lam)){
    penvals <- 2^(1:20)/2^20*n^(-5/9)*100
    # seq(1,diff(range(y)),length.out = 20)* n^(-11/18)
    risk <- sapply(penvals,function(penval){
      c0vec=t(wmat)%*%(cvec-hmat%*%b0-penval*t(D)%*%D%*%b0)
      q0mat=t(wmat)%*%(hmat+penval*t(D)%*%D)%*%wmat
      ans <- quadprog::solve.QP(q0mat,c0vec,t(amat%*%wmat),-amat%*%b0,meq=1)
      alphahat=ans$solution
      bhat=wmat%*%alphahat+b0
      tmpr=sum(sapply(1:10, function(j){xcv <- y[fold!=j];nm <- n-nf[j];cvecf <- makec(xcv,tk,d,h);
      solf <- quadprog::solve.QP(q0mat, t(wmat)%*%(cvecf-hmat%*%b0-penval*t(D)%*%D%*%b0), t(amat%*%wmat),-amat%*%b0,meq=1)
      alphahatf <- solf$solution;bhatf <- wmat%*%alphahatf + b0;(gammay%*%bhatf)[fold==j]}))
      t(bhat)%*%hmat%*%bhat-2/n*tmpr
    })
    rm=1:20*0
    if(risk[1]<risk[2]){rm[1]=1}
    for(i in 2:19){
      if(risk[i-1]>=risk[i]&risk[i+1]>=risk[i]){rm[i]=1}
    }
    if(risk[20]<risk[19]){rm[20]=1}
    #id <- max(which(rm==1))
    riskrange = diff(range(risk))
    id=1
    for(i in sort(which(rm==1),decreasing = TRUE)){
      if( (risk[i]-min(risk)) < riskrange/2) {id=i;break}
    }
    lam <- penvals[id]
    #risk <- risk[id]
  } 
  c0vec=t(wmat)%*%(cvec-hmat%*%b0-lam*t(D)%*%D%*%b0)
  q0mat=t(wmat)%*%(hmat+lam*t(D)%*%D)%*%wmat  
  ans=solve.QP(q0mat,c0vec,t(amat%*%wmat),-amat%*%b0,meq=1)
  
  crit=ans$value
  alphahat=ans$solution
  bhat=wmat%*%alphahat+b0
  #risk <- t(bhat)%*%hmat%*%bhat-2/n*sum(t(gamma)%*%bhatf)
  #alphatil=solve(q0mat)%*%c0vec
  #btil=wmat%*%alphatil+b0
  
  
  if(plot){
    par(mfrow=c(1,2))
    hist(y+me,freq=F,br=40,main="y density")
    lines(yp+me,gamma%*%bhat,col=4,lwd=3)
    plot(support+me,c(0,max(delta%*%bhat)*1.1),type="n",xlab="x",ylab = "Density",main = "x density")
    lines(yp+me,delta%*%bhat,col=4,lwd=3)
    rug(tk+me)
  }
  ans=new.env()
  ans$yp = yp+me
  ans$knots = tk+me
  ans$fhat = delta%*%bhat
  ans$ghat = gamma%*%bhat
  ans$bhat=bhat
  ans$support = support+me
  ans$d = d
  ans$amat=amat
  #  ans$mode=ifelse(mode0,0,yp[which.max(delta%*%bhat)])
  ans$lam=lam
  ans$risk=risk
  ans$id=id
  ans$crit=crit
  ans$me=me
  ans
  
}




####make c vector
makec=function(y,tk,d,h){
  nc=length(y)
  nm=length(tk)-3
  cvec=1:nm
  if(h==0){
    for(j in 1:nm){
      gt=rep(0,nc)
      gt[y>tk[j]&y<=tk[j+1]] = 2*(y[y>tk[j]&y<=tk[j+1]]-tk[j])^2/d^2/3
      gt[y>=tk[j+1]&y<=tk[j+2]] = 1-4*(y[y>=tk[j+1]&y<=tk[j+2]]-tk[j+1]-d/2)^2/d^2/3
      gt[y>tk[j+2]&y<=tk[j+3]] = 2*(y[y>tk[j+2]&y<=tk[j+3]]-tk[j+3])^2/d^2/3
      cvec[j]=sum(gt)/nc
    }
  }else{
    for(j in 1:nm){
      gt=rep(0,nc)
      gt[y>tk[j]-h&y<tk[j+1]-h]=(y[y>tk[j]-h&y<tk[j+1]-h]+h-tk[j])^3/9/h/d^2
      gt[y>=tk[j+1]-h&y<tk[j+2]-h]=-2*(y[y>=tk[j+1]-h&y<tk[j+2]-h]+h-tk[j+1]-d/2)^3/9/h/d^2+(y[y>=tk[j+1]-h&y<tk[j+2]-h]+h-tk[j+1])/2/h+d/h/12
      gt[y>=tk[j+2]-h&y<tk[j+3]-h]=(y[y>=tk[j+2]-h&y<tk[j+3]-h]+h-tk[j+3])^3/9/d/h^2+2*d/3/h
      gt[y>=tk[j+3]-h&y<=tk[j]+h]=2*d/3/h
      gt[y>tk[j]+h&y<=tk[j+1]+h]= -(y[y>tk[j]+h&y<=tk[j+1]+h]-h-tk[j])^3/9/d/h^2+2*d/3/h
      gt[y>tk[j+1]+h&y<=tk[j+2]+h]=2*(y[y>tk[j+1]+h&y<=tk[j+2]+h]-h-tk[j+1]-d/2)^3/9/h/d^2+(-y[y>tk[j+1]+h&y<=tk[j+2]+h]+tk[j+2]+h)/2/h+d/12/h
      gt[y>tk[j+2]+h&y<tk[j+3]+h]=(-y[y>tk[j+2]+h&y<tk[j+3]+h]+tk[j+3]+h)^3/9/h/d^2
      cvec[j]=sum(gt)/nc
    }
  }
  cvec
}
