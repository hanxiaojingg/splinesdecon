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
