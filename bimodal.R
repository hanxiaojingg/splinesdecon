
############## bimodal ###############
## y:       observed data with measurement error ~ U[-h,h]
## h:       size of error density
## support: the support of the true density to be estimated
## d :      length between knots
## plot :   whether to plot the estimated densities
## lam:     the penalty parameter
quaddeconbm <- function(y,h,support=NULL,d=NULL,plot=TRUE,lam=NULL){
  n = length(y)
  if(is.null(support)){
    support=c(0,0)
    dis=quantile(y,0.99)-quantile(y,0.01)
    support[1]=quantile(y,0.01)-dis/3
    support[2]=quantile(y,0.99)+dis/3
  }
  if(!is.null(d)) {
    l = h/d
    if(d>h|round(h/d)!=h/d) d=NULL
  }
  if(is.null(d)){
    l = round((22*n^(1/9))*h/diff(support))
    d = h/l 
  }
  if(l<1) {
    l = 1
    d = h/l
  }
  
  ran=diff(support)
  support[1]=support[1]-(ceiling(ran/d)*d-ran)/2
  support[2]=support[2]+(ceiling(ran/d)*d-ran)/2
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
  wmat = matrix(0,nrow=m,ncol=m-1)
  for(i in 1:(m-1)){
    wmat[i,i] = -1;wmat[i+1,i]=1
  }
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
  amatl=list()
  iter=0
  if(tk[1]<(min(y)+h)){s1=max(which.min(tk<(min(y)+h)),2)}else{s1=2}
  if(tk[m+3]>(max(y)-h)){s2=min(which.max(tk>(max(y)-h)),m)}else{s2=m}
  for(i in s1:(s2-6)){
    for(j in (i+3):(s2-3)){
      for(k in (j+3):(s2)){
        tmpm=matrix(0,m+2,m)
        tmpm[1,1]=1
        tmpm[m+1,m]=1
        for(id in 2:(i-1)) {tmpm[id,id-1]=-1;tmpm[id,id]=1}
        for(id in i:(j-1))  {tmpm[id,id-1]=1;tmpm[id,id]=-1}
        for(id in j:(k-1)) {tmpm[id,id-1]=-1;tmpm[id,id]=1}
        for(id  in k:m) {tmpm[id,id-1]=1;tmpm[id,id]=-1}
        #for(id in (m+2):(2*m-1)) {tmpm[id,id-m]=1}
        tmpm[m+2,j-1]=1
        iter=iter+1
        amatl[[iter]]=tmpm
      }
    }
  }
  c0vec=t(wmat)%*%(cvec-hmat%*%b0-0.1*n^(-5/9)*t(D)%*%D%*%b0)
  q0mat=t(wmat)%*%(hmat+0.1*n^(-5/9)*t(D)%*%D)%*%wmat
  #q0mat=t(wmat)%*%(hmat)%*%wmat
  #for(i in 1:iter){
  #  ans <- quadprog::solve.QP(q0mat,c0vec,t(amatl[[i]]$mat%*%wmat),-amatl[[i]]$mat%*%b0)
  #  crit[i] <- ans$value
  #}
  crit<-lapply(amatl, function(x){ans <- quadprog::solve.QP(q0mat,c0vec,t(x%*%wmat),-x%*%b0);ans$value})
  ind <- which.min(crit)
  amat <- amatl[[ind]]
  nv <- ceiling(n/10)
  fold <- rep(1:10,nv)[1:n]
  ##if(nv<n/10) {fold <- c(fold,1:(n-10*nv))}
  nf <- 1:10
  for(i in 1:10){nf[i] <- sum(fold==i)}
  fold <- sample(fold,n)
  risk=NULL
  id=NULL
  if(is.null(lam)){
    penvals <- 2^(1:20)/2^20*n^(-5/9)
    # seq(1,diff(range(y)),length.out = 20)* n^(-11/18)
    risk <- 1:20*0
    for(i in 1:20){
      c0vec=t(wmat)%*%(cvec-hmat%*%b0-penvals[i]*t(D)%*%D%*%b0)
      q0mat=t(wmat)%*%(hmat+penvals[i]*t(D)%*%D)%*%wmat
      ans <- quadprog::solve.QP(q0mat,c0vec,t(amat%*%wmat),-amat%*%b0)
      alphahat=ans$solution
      bhat=wmat%*%alphahat+b0
      risk[i]=sum(sapply(1:10, function(j){xcv <- y[fold!=j];nm <- n-nf[j];cvecf <- makec(xcv,tk,d,h);
      solf <- quadprog::solve.QP(q0mat, t(wmat)%*%(cvecf-hmat%*%b0-penvals[i]*t(D)%*%D%*%b0), t(amat%*%wmat),-amat%*%b0)
      alphahatf <- solf$solution;bhatf <- wmat%*%alphahatf + b0;sum((gammay%*%bhatf)[fold==j])}))
      risk[i] <- t(bhat)%*%hmat%*%bhat-2/n*risk[i]
    }
    rm=1:20*0
    if(risk[1]<risk[2]){rm[1]=1}
    for(i in 2:19){
      if(risk[i-1]>=risk[i]&risk[i+1]>=risk[i]){rm[i]=1}
    }
    if(risk[20]<risk[19]){rm[20]=1}
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
  ans=solve.QP(q0mat,c0vec,t(amat%*%wmat),-amat%*%b0)
  alphahat=ans$solution
  bhat=wmat%*%alphahat+b0
  crit=ans$value
  #risk <- t(bhat)%*%hmat%*%bhat-2/n*sum(t(gamma)%*%bhatf)
  
  #alphatil=solve(q0mat)%*%c0vec
  #btil=wmat%*%alphatil+b0
  
  
  if(plot){
    par(mfrow=c(1,2))
    hist(y,freq=F,br=40,main="y density")
    lines(yp,gamma%*%bhat,col=4,lwd=3)
    plot(support,c(0,max(delta%*%bhat)+0.05),type="n",xlab="x",ylab = "Density",main = "x density")
    lines(yp,delta%*%bhat,col=4,lwd=3)
    rug(tk)
  }
  ans=new.env()
  ans$yp = yp
  ans$knots = tk
  ans$fhat = delta%*%bhat
  ans$ghat = gamma%*%bhat
  ans$bhat=bhat
  ans$support = support
  ans$d = d
  ans$amat=amat
  ans$lam=lam
  ans$risk=risk
  ans$id=id
  ans$crit=crit
  ans
}