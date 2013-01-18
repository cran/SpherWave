c     Decomposition

      subroutine cov(aa,temp,x,pp0,pp,site1,site2,
     +               eta,eta1,eta2,netlab,JJ,p)
c      dll_export cov
      integer JJ,p
      integer netlab(JJ)
      real*8 pp0,eta1,eta2
      real*8 aa(JJ,JJ),temp(JJ,JJ),x(JJ,JJ),pp(JJ,JJ),
     +       site1(JJ),site2(JJ),eta(p) 
      aa(1,1)=0.0
      x(1,1)=0.0
      temp(1,1)=0.0
      do 10 i=1,JJ
         do 11 j=i,JJ
            x(i,j)=sin(site1(i))*sin(site1(j))
            x(i,j)=x(i,j)+cos(site1(i))*cos(site1(j))*
     +              cos(site2(i)-site2(j))
            temp(i,j)=4*3.141593*1.414214
            eta1=eta(netlab(i))
            pp0=0.0
            do 12 k=1,p
               eta2=eta(k)
               pp0=(1+eta1*eta2)*((1-eta1)*(1-eta2))**2
               pp0=pp0/(1+eta1)/(1+eta2)/(1-eta1*eta2)**2
               if(netlab(j).eq.k) then
                  pp(i,j)=pp0*(1-eta1*eta2)**3/(1-2*eta1*eta2*
     +                    x(i,j)+(eta1*eta2)**2)**1.5 
               endif
 12         continue   
            temp(i,j)=pp(i,j)*4*3.141593*1.414214
            aa(i,j)=temp(i,j)
            aa(j,i)=aa(i,j)
 11      continue
 10   continue
      return
      end

c     BETA1
      subroutine ridge(gg,temp,ssite1,ssite2,site1,site2,
     +                 x,snet,seta,JJ,KK,p,lam)
c      dll_export ridge
      integer JJ,p,snet(JJ)
      real*8 gg(KK+JJ,JJ),temp(JJ,KK),ssite1(JJ),ssite2(JJ),
     +       site1(KK),site2(KK),x(JJ,KK),seta(p),lam
      gg(1,1)=0.0
      do 10 j=1,JJ
         x(1,1)=0.0
         temp(1,1)=0.0
         do 11 k=1,KK
            x(j,k)=sin(ssite1(j))*sin(site1(k))
            x(j,k)=x(j,k)+cos(ssite1(j))*cos(site1(k))
     +               *cos(ssite2(j)-site2(k))
             do 12 i=1,p
                if(snet(j).eq.i) then
                temp(j,k)=(1-seta(i))**3/(1-2*seta(i)*x(j,k)+seta(i)
     +                    **2)**1.5 
             endif   
 12          continue
         gg(k,j)=temp(j,k)
 11      continue
 10   continue
      do 20 i=KK+1,KK+JJ
         do 21 j=1,KK
            if(i-KK.eq.j) then
               gg(i,j)=lam
            endif
 21      continue
 20   continue 
      return
      end

c    BETA1
      subroutine ls(gg,temp,ssite1,ssite2,site1,site2,
     +              x,snet,seta,JJ,KK,p)
c      dll_export ls
      integer JJ,p,snet(JJ)
      real*8 gg(KK,JJ),temp(JJ,KK),ssite1(JJ),
     +       ssite2(JJ),site1(KK),site2(KK),x(JJ,KK),seta(p)
      gg(1,1)=0.0
      do 10 j=1,JJ
         x(1,1)=0.0
         temp(1,1)=0.0
         do 11 k=1,KK
            x(j,k)=sin(ssite1(j))*sin(site1(k))
            x(j,k)=x(j,k)+cos(ssite1(j))*cos(site1(k))
     +               *cos(ssite2(j)-site2(k))
             do 12 i=1,p
                if(snet(j).eq.i) then
                temp(j,k)=(1-seta(i))**3/(1-2*seta(i)*x(j,k)+seta(i)
     +                    **2)**1.5 
             endif   
 12          continue
         gg(k,j)=temp(j,k)
 11      continue
 10   continue
      return
      end

c     BETA1

      subroutine sbf(aa,bb,stemp,temp,spp,pp,ppp,point1,point2,site1,
     +               site2,x,coef,beta,netlab,para,eta,m,n,JJ,p,p0)
c      dll_export sbf
      integer m,n,JJ,p,p0
      integer netlab(JJ),para(JJ)
      real*8 aa(m,n),bb(m,n),stemp,temp(JJ),spp,pp(JJ),ppp(JJ),
     +       point1(n),point2(m),site1(JJ),site2(JJ),x(JJ),
     +       coef(JJ),eta(p),beta
      aa(1,1)=0.0
      bb(1,1)=0.0
      do 10 i=1,m
         do 20 j=1,n
         x(1)=0.0
         temp(1)=0.0
         pp(1)=0.0
         ppp(1)=0.0
         stemp=0.0
         spp=0.0
         do 11 k=1,JJ
            x(k)=sin(point1(j))*sin(site1(k))
            x(k)=x(k)+cos(point1(j))*cos(site1(k))
     +               *cos(point2(i)-site2(k))
             do 12 l=p0+1,p
                if(netlab(k).eq.l) then
                beta=coef(k)
                pp(k)=(1-eta(l))**3/(1-2*eta(l)*x(k)+eta(l)**2)**1.5 
                temp(k)=beta*pp(k)
                ppp(k)=para(k)*pp(k)
             endif   
 12          continue
             spp=spp+ppp(k)
             stemp=stemp+temp(k)
 11      continue
         aa(i,j)=stemp
         bb(i,j)=spp
 20      continue
 10   continue
      return
      end

c     NETWORK
      subroutine disboun(xm,center,dot,dist,min,netlab,theta,phi,
     +                   terri,spot,m,n,l)
      integer m,n,l
      integer spot,netlab(m)
      real*8 dist(m),dot(m),theta,phi,terri,min,xm(m,2),center(n,2)
      min=0.0
      theta=0.0
      phi=0.0
      do 10 i=1,n
         min=3.141593
         spot=0
         dist(1)=0.0
         dot(1)=0.0
         theta=(90-center(i,1))*3.141593/180
         phi=center(i,2)*3.141593/180
         do 20 j=1,m
            dot(j)=cos(theta)*cos((90-xm(j,1))*3.141593/180)+sin(theta)*
     +                       sin((90-xm(j,1))*3.141593/180)*
     +                       cos(phi-(-xm(j,2)*3.141593/180))
            dist(j)=acos(dot(j))
            if(dist(j).le.terri) then
               if(dist(j).le.min) then
                  min=dist(j)
                  spot=j
               else
                  min=min
               endif
            else
               go to 20
            endif   
 20      continue
         netlab(spot)=l
 10   continue
      return
      end
      
c     NETWORK
      subroutine sedisboun(xm,center,dot,dist,min,netlab,theta,phi,
     +                   terri,spot,m,n,l)
      integer m,n,l
      integer spot,netlab(m)
      real*8 dist(m),dot(m),theta,phi,terri,min,xm(m,2),center(n,2)
      min=0.0
      theta=0.0
      phi=0.0
      do 10 i=1,n
         min=3.141593
         spot=0
         dist(1)=0.0
         dot(1)=0.0
         theta=(90-center(i,1))*3.141593/180
         phi=center(i,2)*3.141593/180
         do 20 j=1,m
         if(netlab(j).eq.0) then
            dot(j)=cos(theta)*cos((90-xm(j,1))*3.141593/180)+sin(theta)*
     +                       sin((90-xm(j,1))*3.141593/180)*
     +                       cos(phi-(-xm(j,2)*3.141593/180))
            dist(j)=acos(dot(j))
c            if(netlab(j).eq.0) then
            if(dist(j).le.terri) then
               if(dist(j).le.min) then
                  min=dist(j)
                  spot=j
               else
                  min=min
               endif
             else
               go to 20
             endif
         endif   
 20      continue
         netlab(spot)=l
 10   continue
      return
      end
      
c     NETWORK
      subroutine seselc(xm,ym,dot,dist,netlab,theta,phi,
     +                   bb,m,n,l)
      integer m,n,l
      integer netlab(m)
      real*8 dist(m),dot(m),theta,phi,bb,xm(m,2),ym(n,2)
      theta=0.0
      phi=0.0
      do 10 i=1,n
         dist(1)=0.0
         dot(1)=0.0
         theta=(90-ym(i,1))*3.141593/180
         phi=-ym(i,2)*3.141593/180
         do 20 j=1,m
         if(netlab(j).eq.l) then   
            dot(j)=cos(theta)*cos((90-xm(j,1))*3.141593/180)+sin(theta)*
     +                       sin((90-xm(j,1))*3.141593/180)*
     +                       cos(phi-(-xm(j,2)*3.141593/180))
            dist(j)=acos(dot(j))
            if(dist(j).le.bb) then
               netlab(j)=0
            else
               go to 20
            endif   
         endif   
 20      continue
 10   continue
      return
      end
