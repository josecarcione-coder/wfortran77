C-----------------------------------------------------------------------        
C       PRIME FACTOR FFT ROUTINE BASED ON TEMPERTON ALGORITHM                
C       J. COMPUT. PHYS. VOL 58, 1985, PAGE 283.                                             
C                                                                               
C       VARIABLES:                                                              
C                 N -       TRANSFORM LENGTH                                    
C                 A -       REAL PART                    
C                 B -       IMAGINARY PART                    
C                 IFAX -    ARRAY CONTAINING FACTORIZATIONS OF N        
C                           ACCORDING TO N=n1*n2*....*nk                        
C                 NFAC -    NUMBER OF FACTORS IN N          
C                 ISKIP -   STRIDE OF FFT (e.g IF=2 SKIP EVERY SECOND           
C                           SAMPLE)                                             
C                 ISIGN -   FFT SIGN                                            
C-----------------------------------------------------------------------        
C      File: p11.f                                                            
C      Date: Nov 1985                                                        
C      Written by Jose' M. Carcione.                                               
       subroutine pfft(a,b,n,ifax,nfac,iskip,isign)                             
       integer ifax(1)                                                          
       dimension a(1),b(1),c(20),d(20),e(20),f(20)                              
       if(n.eq.1) return                                                        
       data sin60/0.8660254038/,sin72/0.9510565163/,sin36/0.5877852523/         
       data sq54/0.5590169944/                                                  
       data c71/0.6234898018/,c72/-0.2225209342/,c73/-0.9009688680/             
       data s71/0.7818314825/,s72/0.9749279121/,s73/0.4338837388/              
       data c91/0.7660444431/,c92/0.1736481775/,c94/-0.939692621/               
       data s91/0.6427876097/,s92/0.9848077530/,s94/0.3420201430/               
       data c111/0.8412535328/,c112/0.4154150129/,c113/-0.1423148385/           
       data c114/-0.6548607342/,c115/-0.9594929737/                             
       data s111/0.5406408175/,s112/0.9096319954/,s113/0.9898214418/            
       data s114/0.7557495742/,s115/0.2817325565/                               
       sr2=1./sqrt(2.)                                                         
C                                                                               
       ns=n*iskip                                                               
        do 1000 k=1,nfac                                                        
C                                                                               
        ifac=ifax(k)                                                            
        m=n/ifac                                                                
       do 100 j=1,ifac                                                          
        mu=j                                                                    
       mm=j*m                                                                   
       if(mod(mm,ifac).eq.1) go to 110                                          
100     continue                                                                
110     mm=mm*iskip                                                             
C                                                                               
C       PERFORM THE INVERSE TRANSFORM                                           
        if(isign.eq.-1) mu=ifac-mu                                              
C       -----------------------------                                           
C       MU IS THE REQUIRED ROTATION FOR THE DFT MODULE OF ORDER IFAC            
C       NOW COMPUTE THE ADDRESSES IA,IB ETC. AND SELECT THE CODING              
C       FOR THE CURRENT FACTOR                                                  
C                                                                               
       ia=1                                                                     
        iaa=1                                                                   
       ib=ia+mm                                                                 
        ibb=ib+ib-1                                                             
        if(ifac.eq.2) go to 200                                                 
       ic=ib+mm                                                                 
        if(ic.gt.ns)ic=ic-ns                                                    
        icc=ic+ic-1                                                             
       if(ifac.eq.3) go to 300                                                  
       id=ic+mm                                                                 
       if(id.gt.ns)id=id-ns                                                     
        idd=id+id-1                                                             
       if(ifac.eq.4) go to 400                                                  
       ie=id+mm                                                                 
       if(ie.gt.ns)ie=ie-ns                                                     
        iee=ie+ie-1                                                             
       if(ifac.eq.5) go to 500                                                  
       if=ie+mm                                                                 
        if(if.gt.ns)if=if-ns                                                    
       iff=if+if-1                                                              
        if(ifac.eq.6) go to 600                                                 
        ig=if+mm                                                                
        if(ig.gt.ns)ig=ig-ns                                                    
       igg=ig+ig-1                                                              
       if(ifac.eq.7) go to 700                                                  
       ih=ig+mm                                                                 
       if(ih.gt.ns)ih=ih-ns                                                     
       ihh=ih+ih-1                                                              
       if(ifac.eq.8) go to 800                                                  
       ii=ih+mm                                                                 
       if(ii.gt.ns)ii=ii-ns                                                     
       iii=ii+ii-1                                                              
       if(ifac.eq.9) go to 900                                                  
       ij=ii+mm                                                                 
       if(ij.gt.ns)ij=ij-ns                                                     
       ijj=ij+ij-1                                                              
       if(ifac.eq.10) stop 'no module for factor 10'                            
       ik=ij+mm                                                                 
       if(ik.gt.ns)ik=ik-ns                                                     
       ikk=ik+ik-1                                                              
       if(ifac.eq.11) go to 1100                                                
       stop 'no modules for factors greater than 11'                            
C                                                                               
C       CODING FOR FACTOR 2                                                     
C       --------------------                                                    
200     continue                                                                
       do 210 l=1,m                                                             
        aa=a(iaa)                                                               
        ba=b(iaa)                                                               
        ab=a(ibb)                                                               
        bb=b(ibb)                                                               
       a(iaa)=aa+ab                                                             
       b(iaa)=ba+bb                                                             
       a(ibb)=aa-ab                                                             
       b(ibb)=ba-bb                                                             
        ix=ib+iskip                                                             
        ib=ia+iskip                                                             
       ia=ix                                                                    
       iaa=ia+ia-1                                                              
       ibb=ib+ib-1                                                              
210     continue                                                                
        go to 1000                                                              
C                                                                               
C       CODING FOR FACTOR 3                                                     
C       --------------------                                                    
300     continue                                                                
       z3=sin60                                                                 
       if(mu.eq.2)z3=-z3                                                        
       do 310 l=1,m                                                             
        aa=a(iaa)                                                               
        ba=b(iaa)                                                               
        ab=a(ibb)                                                               
       bb=b(ibb)                                                                
        ac=a(icc)                                                               
        bc=b(icc)                                                               
       t1r=ab+ac                                                                
       t1i=bb+bc                                                                
       t2r=aa-0.5*t1r                                                           
       t2i=ba-0.5*t1i                                                           
       t3r=z3*(ab-ac)                                                           
       t3i=z3*(bb-bc)                                                           
       a(iaa)=aa+t1r                                                            
       b(iaa)=ba+t1i                                                            
       a(ibb)=t2r-t3i                                                           
       b(ibb)=t2i+t3r                                                           
        a(icc)=t2r+t3i                                                          
       b(icc)=t2i-t3r                                                           
       ix=ic+iskip                                                              
       ic=ib+iskip                                                              
       ib=ia+iskip                                                              
       ia=ix                                                                    
       iaa=ia+ia-1                                                              
       ibb=ib+ib-1                                                              
       icc=ic+ic-1                                                              
310     continue                                                                
       go to 1000                                                               
C                                                                               
C       CODING FOR FACTOR 4                                                     
C       --------------------                                                    
400     continue                                                                
       z4=1.0                                                                   
       if(mu.eq.3)z4=-z4                                                        
       do 410 l=1,m                                                             
       aa=a(iaa)                                                                
       ba=b(iaa)                                                                
       ab=a(ibb)                                                                
       bb=b(ibb)                                                                
       ac=a(icc)                                                                
       bc=b(icc)                                                                
       ad=a(idd)                                                                
       bd=b(idd)                                                                
       t1r=aa+ac                                                                
       t1i=ba+bc                                                                
       t2r=ab+ad                                                                
       t2i=bb+bd                                                                
       t3r=aa-ac                                                                
       t3i=ba-bc                                                                
       t4r=z4*(ab-ad)                                                           
       t4i=z4*(bb-bd)                                                           
       a(iaa)=t1r+t2r                                                           
       b(iaa)=t1i+t2i                                                           
       a(ibb)=t3r-t4i                                                           
       b(ibb)=t3i+t4r                                                           
       a(icc)=t1r-t2r                                                           
       b(icc)=t1i-t2i                                                           
       a(idd)=t3r+t4i                                                           
       b(idd)=t3i-t4r                                                           
       ix=id+iskip                                                              
       id=ic+iskip                                                              
       ic=ib+iskip                                                              
       ib=ia+iskip                                                              
       ia=ix                                                                    
       iaa=ia+ia-1                                                              
       ibb=ib+ib-1                                                              
       icc=ic+ic-1                                                              
       idd=id+id-1                                                              
410     continue                                                                
        go to 1000                                                              
C                                                                               
C       CODING FOR FACTOR 5                                                     
C       --------------------                                                    
500     continue                                                                
       go to(501,502,503,504) mu                                                
501     c1=sq54                                                                 
       c2=sin72                                                                 
       c3=sin36                                                                 
        go to 505                                                               
502     c1=-sq54                                                                
       c2=sin36                                                                 
       c3=-sin72                                                                
       go to 505                                                                
503     c1=-sq54                                                                
       c2=-sin36                                                                
       c3=sin72                                                                 
       go to 505                                                                
504     c1=sq54                                                                 
       c2=-sin72                                                                
       c3=-sin36                                                                
505     do 510 l=1,m                                                            
       aa=a(iaa)                                                                
       ba=b(iaa)                                                                
       ab=a(ibb)                                                                
       bb=b(ibb)                                                                
       ac=a(icc)                                                                
       bc=b(icc)                                                                
       ad=a(idd)                                                                
       bd=b(idd)                                                                
       ae=a(iee)                                                                
       be=b(iee)                                                                
        t1r=ab+ae                                                               
        t1i=bb+be                                                               
       t2r=ac+ad                                                                
       t2i=bc+bd                                                                
       t3r=ab-ae                                                                
       t3i=bb-be                                                                
       t4r=ac-ad                                                                
       t4i=bc-bd                                                                
       t5r=t1r+t2r                                                              
       t5i=t1i+t2i                                                              
       t6r=c1*(t1r-t2r)                                                         
       t6i=c1*(t1i-t2i)                                                         
       t7r=aa-0.25*t5r                                                          
       t7i=ba-0.25*t5i                                                          
       t8r=t7r+t6r                                                              
       t8i=t7i+t6i                                                              
       t9r=t7r-t6r                                                              
       t9i=t7i-t6i                                                              
       t10r=c2*t3r+c3*t4r                                                       
       t10i=c2*t3i+c3*t4i                                                       
       t11r=c3*t3r-c2*t4r                                                       
       t11i=c3*t3i-c2*t4i                                                       
       a(iaa)=aa+t5r                                                            
       b(iaa)=ba+t5i                                                            
       a(ibb)=t8r-t10i                                                          
       b(ibb)=t8i+t10r                                                          
       a(icc)=t9r-t11i                                                          
       b(icc)=t9i+t11r                                                          
       a(idd)=t9r+t11i                                                          
       b(idd)=t9i-t11r                                                          
       a(iee)=t8r+t10i                                                          
       b(iee)=t8i-t10r                                                          
       ix=ie+iskip                                                              
       ie=id+iskip                                                              
       id=ic+iskip                                                              
       ic=ib+iskip                                                              
       ib=ia+iskip                                                              
       ia=ix                                                                    
       iaa=ia+ia-1                                                              
       ibb=ib+ib-1                                                              
       icc=ic+ic-1                                                              
       idd=id+id-1                                                              
       iee=ie+ie-1                                                              
510     continue                                                                
        go to 1000                                                              
C                                                                               
C       CODING FOR FACTOR 6                                                     
C       --------------------                                                    
600     continue                                                                
        z3=sin60                                                                
        if(mu.eq.5)z3=-z3                                                       
        do 610 l=1,m                                                            
       aa=a(iaa)                                                                
       ba=b(iaa)                                                                
       ab=a(ibb)                                                                
       bb=b(ibb)                                                                
       ac=a(icc)                                                                
       bc=b(icc)                                                                
       ad=a(idd)                                                                
       bd=b(idd)                                                                
       ae=a(iee)                                                                
       be=b(iee)                                                                
       af=a(iff)                                                                
       bf=b(iff)                                                                
        t1r=aa+ad                                                               
        t1i=ba+bd                                                               
        t2r=ac+af                                                               
        t2i=bc+bf                                                               
        t3r=ae+ab                                                               
        t3i=be+bb                                                               
       t4r=aa-ad                                                                
       t4i=ba-bd                                                                
       t5r=ac-af                                                                
       t5i=bc-bf                                                                
       t6r=ae-ab                                                                
       t6i=be-bb                                                                
       t7r=t2r+t3r                                                              
       t7i=t2i+t3i                                                              
        t8r=t2r-t3r                                                             
       t8i=t2i-t3i                                                              
       t9r=t5r+t6r                                                              
       t9i=t5i+t6i                                                              
       t10r=t5r-t6r                                                             
       t10i=t5i-t6i                                                             
       t11r=t1r-0.5*t7r                                                         
       t11i=t1i-0.5*t7i                                                         
       t12r=t4r-0.5*t9r                                                         
        t12i=t4i-0.5*t9i                                                        
       t13r=z3*t8r                                                              
       t13i=z3*t8i                                                              
       t14r=z3*t10r                                                             
       t14i=z3*t10i                                                             
       a(iaa)=t1r+t7r                                                           
       b(iaa)=t1i+t7i                                                           
       a(ibb)=t12r-t14i                                                         
       b(ibb)=t12i+t14r                                                         
       a(icc)=t11r+t13i                                                         
       b(icc)=t11i-t13r                                                         
       a(idd)=t4r+t9r                                                           
       b(idd)=t4i+t9i                                                           
       a(iee)=t11r-t13i                                                         
       b(iee)=t11i+t13r                                                         
       a(iff)=t12r+t14i                                                         
       b(iff)=t12i-t14r                                                         
       ix=if+iskip                                                              
       if=ie+iskip                                                              
       ie=id+iskip                                                              
       id=ic+iskip                                                              
       ic=ib+iskip                                                              
       ib=ia+iskip                                                              
       ia=ix                                                                    
       iaa=ia+ia-1                                                              
       ibb=ib+ib-1                                                              
       icc=ic+ic-1                                                              
       idd=id+id-1                                                              
       iee=ie+ie-1                                                              
       iff=if+if-1                                                              
610     continue                                                                
        go to 1000                                                              
C                                                                               
C       CODING FOR FACTOR 7                                                     
C       --------------------                                                    
700     continue                                                                
        go to(701,702,703,704,705,706) mu                                       
701     c1=c71                                                                  
       c2=c72                                                                   
       c3=c73                                                                   
       c4=c73                                                                   
       c6=c71                                                                   
       c9=c72                                                                   
       s1=s71                                                                   
       s2=s72                                                                   
       s3=s73                                                                   
       s4=-s73                                                                  
       s6=-s71                                                                  
       s9=s72                                                                   
       go to 707                                                                
702     c1=c72                                                                  
       c2=c73                                                                   
       c3=c71                                                                   
       c4=c71                                                                   
       c6=c72                                                                   
       c9=c73                                                                   
       s1=s72                                                                   
       s2=-s73                                                                  
        s3=-s71                                                                 
       s4=s71                                                                   
       s6=-s72                                                                  
       s9=-s73                                                                  
       go to 707                                                                
703     c1=c73                                                                  
       c2=c71                                                                   
       c3=c72                                                                   
       c4=c72                                                                   
       c6=c73                                                                   
       c9=c71                                                                   
       s1=s73                                                                   
       s2=-s71                                                                  
       s3=s72                                                                   
       s4=-s72                                                                  
       s6=-s73                                                                  
       s9=-s71                                                                  
       go to 707                                                                
704     c1=c73                                                                  
       c2=c71                                                                   
       c3=c72                                                                   
       c4=c72                                                                   
       c6=c73                                                                   
       c9=c71                                                                   
       s1=-s73                                                                  
       s2=s71                                                                   
       s3=-s72                                                                  
       s4=s72                                                                   
       s6=s73                                                                   
       s9=s71                                                                   
        go to 707                                                               
705       c1=c72                                                                
       c2=c73                                                                   
       c3=c71                                                                   
       c4=c71                                                                   
       c6=c72                                                                   
       c9=c73                                                                   
       s1=-s72                                                                  
       s2=s73                                                                   
       s3=s71                                                                   
       s4=-s71                                                                  
       s6=s72                                                                   
       s9=s73                                                                   
       go to 707                                                                
706     c1=c71                                                                  
       c2=c72                                                                   
       c3=c73                                                                   
       c4=c73                                                                   
       c6=c71                                                                   
       c9=c72                                                                   
       s1=-s71                                                                  
       s2=-s72                                                                  
       s3=-s73                                                                  
       s4=s73                                                                   
       s6=s71                                                                   
       s9=-s72                                                                  
707     do 710 l=1,m                                                            
       aa=a(iaa)                                                                
       ba=b(iaa)                                                                
       ab=a(ibb)                                                                
       bb=b(ibb)                                                                
       ac=a(icc)                                                                
       bc=b(icc)                                                                
       ad=a(idd)                                                                
       bd=b(idd)                                                                
       ae=a(iee)                                                                
       be=b(iee)                                                                
       af=a(iff)                                                                
       bf=b(iff)                                                                
       ag=a(igg)                                                                
       bg=b(igg)                                                                
        t1r=ab+ag                                                               
       t1i=bb+bg                                                                
       t2r=ac+af                                                                
       t2i=bc+bf                                                                
       t3r=ad+ae                                                                
       t3i=bd+be                                                                
       t4r=ab-ag                                                                
       t4i=bb-bg                                                                
       t5r=ac-af                                                                
       t5i=bc-bf                                                                
       t6r=ad-ae                                                                
       t6i=bd-be                                                                
       r1r=aa+t1r*c1+t2r*c2+t3r*c3                                              
       r1i=ba+t1i*c1+t2i*c2+t3i*c3                                              
       r2r=t4r*s1+t5r*s2+t6r*s3                                                 
       r2i=t4i*s1+t5i*s2+t6i*s3                                                 
       r3r=aa+t1r*c2+t2r*c4+t3r*c6                                              
       r3i=ba+t1i*c2+t2i*c4+t3i*c6                                              
       r4r=t4r*s2+t5r*s4+t6r*s6                                                 
       r4i=t4i*s2+t5i*s4+t6i*s6                                                 
       r5r=aa+t1r*c3+t2r*c6+t3r*c9                                              
       r5i=ba+t1i*c3+t2i*c6+t3i*c9                                              
       r6r=t4r*s3+t5r*s6+t6r*s9                                                 
       r6i=t4i*s3+t5i*s6+t6i*s9                                                 
       a(iaa)=aa+t1r+t2r+t3r                                                    
       b(iaa)=ba+t1i+t2i+t3i                                                    
       a(ibb)=r1r-r2i                                                           
       b(ibb)=r1i+r2r                                                           
       a(icc)=r3r-r4i                                                           
       b(icc)=r3i+r4r                                                           
       a(idd)=r5r-r6i                                                           
       b(idd)=r5i+r6r                                                           
       a(iee)=r5r+r6i                                                           
       b(iee)=r5i-r6r                                                           
       a(iff)=r3r+r4i                                                           
       b(iff)=r3i-r4r                                                           
       a(igg)=r1r+r2i                                                           
       b(igg)=r1i-r2r                                                           
       ix=ig+iskip                                                              
       ig=if+iskip                                                              
       if=ie+iskip                                                              
       ie=id+iskip                                                              
       id=ic+iskip                                                              
       ic=ib+iskip                                                              
       ib=ia+iskip                                                              
       ia=ix                                                                    
       iaa=ia+ia-1                                                              
       ibb=ib+ib-1                                                              
       icc=ic+ic-1                                                              
       idd=id+id-1                                                              
       iee=ie+ie-1                                                              
       iff=if+if-1                                                              
       igg=ig+ig-1                                                              
710     continue                                                                
        go to 1000                                                              
C                                                                               
C       CODING FOR FACTOR 8                                                     
C       --------------------                                                    
800     continue                                                                
        z4=1.                                                                   
       if(mu.eq.3.or.mu.eq.7)z4=-z4                                             
       mu=(mu+1)/2                                                              
       go to (801,802,803,804)mu                                                
801     d1=sr2                                                                  
       d2=sr2                                                                   
       d4=1                                                                     
       d5=-sr2                                                                  
       d6=sr2                                                                   
       go to 805                                                                
802     d1=-sr2                                                                 
       d2=sr2                                                                   
       d4=-1                                                                    
       d5=sr2                                                                   
       d6=sr2                                                                   
       go to 805                                                                
803     d1=-sr2                                                                 
       d2=-sr2                                                                  
       d4=1                                                                     
       d5=sr2                                                                   
       d6=-sr2                                                                  
       go to 805                                                                
804     d1=sr2                                                                  
       d2=-sr2                                                                  
       d4=-1                                                                    
       d5=-sr2                                                                  
       d6=-sr2                                                                  
805       do 810 l=1,m                                                          
       e(1)=a(iaa)                                                              
       f(1)=b(iaa)                                                              
       e(3)=a(ibb)                                                              
       f(3)=b(ibb)                                                              
       e(5)=a(icc)                                                              
       f(5)=b(icc)                                                              
       e(7)=a(idd)                                                              
       f(7)=b(idd)                                                              
       e(9)=a(iee)                                                              
       f(9)=b(iee)                                                              
       e(11)=a(iff)                                                             
       f(11)=b(iff)                                                             
       e(13)=a(igg)                                                             
       f(13)=b(igg)                                                             
       e(15)=a(ihh)                                                             
       f(15)=b(ihh)                                                             
       aa=e(1)                                                                  
       ba=f(1)                                                                  
       ab=e(5)                                                                  
       bb=f(5)                                                                  
       ac=e(9)                                                                  
       bc=f(9)                                                                  
       ad=e(13)                                                                 
       bd=f(13)                                                                 
       t1r=aa+ac                                                                
       t1i=ba+bc                                                                
       t2r=ab+ad                                                                
       t2i=bb+bd                                                                
       t3r=aa-ac                                                                
       t3i=ba-bc                                                                
       t4r=z4*(ab-ad)                                                           
       t4i=z4*(bb-bd)                                                           
       c(1)=t1r+t2r                                                             
       d(1)=t1i+t2i                                                             
       c(3)=t3r-t4i                                                             
       d(3)=t3i+t4r                                                             
       c(5)=t1r-t2r                                                             
        d(5)=t1i-t2i                                                            
       c(7)=t3r+t4i                                                             
       d(7)=t3i-t4r                                                             
       aa=e(3)                                                                  
       ba=f(3)                                                                  
       ab=e(7)                                                                  
       bb=f(7)                                                                  
       ac=e(11)                                                                 
       bc=f(11)                                                                 
       ad=e(15)                                                                 
       bd=f(15)                                                                 
       t1r=aa+ac                                                                
       t1i=ba+bc                                                                
       t2r=ab+ad                                                                
       t2i=bb+bd                                                                
       t3r=aa-ac                                                                
       t3i=ba-bc                                                                
       t4r=z4*(ab-ad)                                                           
       t4i=z4*(bb-bd)                                                           
       c(9)=t1r+t2r                                                             
        d(9)=t1i+t2i                                                            
       x1r=t3r-t4i                                                              
       x1i=t3i+t4r                                                              
       x2r=t1r-t2r                                                              
       x2i=t1i-t2i                                                              
       x3r=t3r+t4i                                                              
       x3i=t3i-t4r                                                              
       c(11)=d1*x1r-d2*x1i                                                      
       d(11)=d2*x1r+d1*x1i                                                      
       c(13)=-d4*x2i                                                            
       d(13)=d4*x2r                                                             
       c(15)=d5*x3r-d6*x3i                                                      
        d(15)=d6*x3r+d5*x3i                                                     
       i=1                                                                      
       j=1                                                                      
       do 820 jt=1,4                                                             
       aa=c(i)                                                                  
       ba=d(i)                                                                  
        ikb=i+8                                                                 
       ab=c(ikb)                                                                
       bb=d(ikb)                                                                
       e(j)=aa+ab                                                               
       f(j)=ba+bb                                                               
       jjb=j+8                                                                  
       e(jjb)=aa-ab                                                             
       f(jjb)=ba-bb                                                             
       i=i+2                                                                    
       j=j+2                                                                    
820     continue                                                                
       a(iaa)=e(1)                                                              
       b(iaa)=f(1)                                                              
       a(ibb)=e(3)                                                              
       b(ibb)=f(3)                                                              
       a(icc)=e(5)                                                              
       b(icc)=f(5)                                                              
       a(idd)=e(7)                                                              
       b(idd)=f(7)                                                              
       a(iee)=e(9)                                                              
       b(iee)=f(9)                                                              
       a(iff)=e(11)                                                             
       b(iff)=f(11)                                                             
       a(igg)=e(13)                                                             
       b(igg)=f(13)                                                             
       a(ihh)=e(15)                                                             
       b(ihh)=f(15)                                                             
       ix=ih+iskip                                                              
       ih=ig+iskip                                                              
       ig=if+iskip                                                              
       if=ie+iskip                                                              
       ie=id+iskip                                                              
       id=ic+iskip                                                              
       ic=ib+iskip                                                              
       ib=ia+iskip                                                              
       ia=ix                                                                    
       iaa=ia+ia-1                                                              
       ibb=ib+ib-1                                                              
       icc=ic+ic-1                                                              
       idd=id+id-1                                                              
       iee=ie+ie-1                                                              
       iff=if+if-1                                                              
       igg=ig+ig-1                                                              
       ihh=ih+ih-1                                                              
810     continue                                                                
        go to 1000                                                              
C                                                                               
C       CODING FOR FACTOR 9                                                     
C       --------------------                                                    
900     continue                                                                
       z3=sin60                                                                 
       if(mu.eq.2.or.mu.eq.5.or.mu.eq.8)z3=-z3                                  
       go to (901,902,903,904,905,906,907,908)mu                                
901     d1=c91                                                                  
       d2=s91                                                                   
       d3=c92                                                                   
       d4=s92                                                                   
       d7=c94                                                                   
       d8=s94                                                                   
       go to 911                                                                
902     d1=c92                                                                  
       d2=s92                                                                   
       d3=c94                                                                   
       d4=s94                                                                   
       d7=c91                                                                   
       d8=-s91                                                                  
       go to 911                                                                
903     continue                                                                
904     d1=c94                                                                  
       d2=s94                                                                   
       d3=c91                                                                   
       d4=-s91                                                                  
       d7=c92                                                                   
       d8=-s92                                                                  
       go to 911                                                                
905     d1=c94                                                                  
       d2=-s94                                                                  
       d3=c91                                                                   
       d4=s91                                                                   
       d7=c92                                                                   
       d8=s92                                                                   
       go to 911                                                                
906     continue                                                                
907     d1=c92                                                                  
       d2=-s92                                                                  
       d3=c94                                                                   
       d4=-s94                                                                  
       d7=c91                                                                   
       d8=s91                                                                   
       go to 911                                                                
908     d1=c91                                                                  
       d2=-s91                                                                  
       d3=c92                                                                   
       d4=-s92                                                                  
       d7=c94                                                                   
       d8=-s94                                                                  
911       do 910 l=1,m                                                          
       e(1)=a(iaa)                                                              
       f(1)=b(iaa)                                                              
       e(3)=a(ibb)                                                              
       f(3)=b(ibb)                                                              
       e(5)=a(icc)                                                              
       f(5)=b(icc)                                                              
       e(7)=a(idd)                                                              
       f(7)=b(idd)                                                              
       e(9)=a(iee)                                                              
       f(9)=b(iee)                                                              
       e(11)=a(iff)                                                             
       f(11)=b(iff)                                                             
       e(13)=a(igg)                                                             
       f(13)=b(igg)                                                             
       e(15)=a(ihh)                                                             
       f(15)=b(ihh)                                                             
       e(17)=a(iii)                                                             
       f(17)=b(iii)                                                             
       aa=e(1)                                                                  
       ba=f(1)                                                                  
       ab=e(7)                                                                  
       bb=f(7)                                                                  
       ac=e(13)                                                                 
       bc=f(13)                                                                 
       t1r=ab+ac                                                                
       t1i=bb+bc                                                                
       t2r=aa-0.5*t1r                                                           
       t2i=ba-0.5*t1i                                                           
       t3r=z3*(ab-ac)                                                           
       t3i=z3*(bb-bc)                                                           
       c(1)=aa+t1r                                                              
       d(1)=ba+t1i                                                              
       c(3)=t2r-t3i                                                             
       d(3)=t2i+t3r                                                             
       c(5)=t2r+t3i                                                             
       d(5)=t2i-t3r                                                             
       i=3                                                                      
       j=7                                                                      
       do 912 jt=1,2                                                             
       aa=e(i)                                                                  
       ba=f(i)                                                                  
       ikb=i+6                                                                  
       ab=e(ikb)                                                                
       bb=f(ikb)                                                                
       ikc=i+12                                                                 
       ac=e(ikc)                                                                
       bc=f(ikc)                                                                
       t1r=ab+ac                                                                
       t1i=bb+bc                                                                
       t2r=aa-0.5*t1r                                                           
       t2i=ba-0.5*t1i                                                           
       t3r=z3*(ab-ac)                                                           
       t3i=z3*(bb-bc)                                                           
       c(j)=aa+t1r                                                              
       d(j)=ba+t1i                                                              
       x1r=t2r-t3i                                                              
       x1i=t2i+t3r                                                              
       x2r=t2r+t3i                                                              
       x2i=t2i-t3r                                                              
       jjb=j+2                                                                  
       if(jt.eq.2) go to 913                                                     
       c(jjb)=d1*x1r-d2*x1i                                                     
       d(jjb)=d2*x1r+d1*x1i                                                     
       jjc=j+4                                                                  
       c(jjc)=d3*x2r-d4*x2i                                                     
       d(jjc)=d4*x2r+d3*x2i                                                     
       go to 914                                                                
913       c(jjb)=d3*x1r-d4*x1i                                                  
       d(jjb)=d4*x1r+d3*x1i                                                     
       jjc=j+4                                                                  
       c(jjc)=d7*x2r-d8*x2i                                                     
       d(jjc)=d8*x2r+d7*x2i                                                     
914     i=i+2                                                                   
       j=j+6                                                                    
912     continue                                                                
       i=1                                                                      
       j=1                                                                      
       do 915 jt=1,3                                                             
       aa=c(i)                                                                  
       ba=d(i)                                                                  
       ikb=i+6                                                                  
       ab=c(ikb)                                                                
       bb=d(ikb)                                                                
       ikc=i+12                                                                 
       ac=c(ikc)                                                                
       bc=d(ikc)                                                                
       t1r=ab+ac                                                                
       t1i=bb+bc                                                                
       t2r=aa-0.5*t1r                                                           
       t2i=ba-0.5*t1i                                                           
       t3r=z3*(ab-ac)                                                           
       t3i=z3*(bb-bc)                                                           
       e(j)=aa+t1r                                                              
       f(j)=ba+t1i                                                              
       jjb=j+6                                                                  
       e(jjb)=t2r-t3i                                                           
       f(jjb)=t2i+t3r                                                           
       jjc=j+12                                                                 
       e(jjc)=t2r+t3i                                                           
       f(jjc)=t2i-t3r                                                           
       i=i+2                                                                    
       j=j+2                                                                    
915     continue                                                                
       a(iaa)=e(1)                                                              
       b(iaa)=f(1)                                                              
       a(ibb)=e(3)                                                              
       b(ibb)=f(3)                                                              
       a(icc)=e(5)                                                              
       b(icc)=f(5)                                                              
       a(idd)=e(7)                                                              
       b(idd)=f(7)                                                              
       a(iee)=e(9)                                                              
       b(iee)=f(9)                                                              
       a(iff)=e(11)                                                             
       b(iff)=f(11)                                                             
       a(igg)=e(13)                                                             
       b(igg)=f(13)                                                             
       a(ihh)=e(15)                                                             
       b(ihh)=f(15)                                                             
       a(iii)=e(17)                                                             
       b(iii)=f(17)                                                             
       ix=ii+iskip                                                              
       ii=ih+iskip                                                              
       ih=ig+iskip                                                              
       ig=if+iskip                                                              
       if=ie+iskip                                                              
       ie=id+iskip                                                              
       id=ic+iskip                                                              
       ic=ib+iskip                                                              
       ib=ia+iskip                                                              
       ia=ix                                                                    
       iaa=ia+ia-1                                                              
       ibb=ib+ib-1                                                              
       icc=ic+ic-1                                                              
       idd=id+id-1                                                              
       iee=ie+ie-1                                                              
       iff=if+if-1                                                              
       igg=ig+ig-1                                                              
       ihh=ih+ih-1                                                              
       iii=ii+ii-1                                                              
910     continue                                                                
        go to 1000                                                              
C                                                                               
C       CODING FOR FACTOR 11                                                    
C       _____________________                                                   
1100    continue                                                                
       go to (1,2,3,4,5,6,7,8,9,10)mu                                           
1       c1=c111                                                                 
       c2=c112                                                                  
       c3=c113                                                                  
       c4=c114                                                                  
       c5=c115                                                                  
       c6=c115                                                                  
       c8=c113                                                                  
       c9=c112                                                                  
       c10=c111                                                                 
       c12=c111                                                                 
       c15=c114                                                                 
       c16=c115                                                                 
       c20=c112                                                                 
       c25=c113                                                                 
       s1=s111                                                                  
       s2=s112                                                                  
       s3=s113                                                                  
       s4=s114                                                                  
       s5=s115                                                                  
       s6=-s115                                                                 
       s8=-s113                                                                 
       s9=-s112                                                                 
       s10=-s111                                                                
       s12=s111                                                                 
       s15=s114                                                                 
       s16=s115                                                                 
       s20=-s112                                                                
       s25=s113                                                                 
       go to 11                                                                 
2       c1=c112                                                                 
       c2=c114                                                                  
       c3=c115                                                                  
       c4=c113                                                                  
       c5=c111                                                                  
       c6=c111                                                                  
       c8=c115                                                                  
       c9=c114                                                                  
       c10=c112                                                                 
       c12=c112                                                                 
       c15=c113                                                                 
       c16=c111                                                                 
       c20=c114                                                                 
       c25=c115                                                                 
       s1=s112                                                                  
       s2=s114                                                                  
       s3=-s115                                                                 
       s4=-s113                                                                 
       s5=-s111                                                                 
       s6=s111                                                                  
       s8=s115                                                                  
       s9=-s114                                                                 
       s10=-s112                                                                
       s12=s112                                                                 
       s15=-s113                                                                
       s16=-s111                                                                
       s20=-s114                                                                
       s25=-s115                                                                
       go to 11                                                                 
3       c1=c113                                                                 
       c2=c115                                                                  
       c3=c112                                                                  
       c4=c111                                                                  
       c5=c114                                                                  
       c6=c114                                                                  
       c8=c112                                                                  
       c9=c115                                                                  
       c10=c113                                                                 
       c12=c113                                                                 
       c15=c111                                                                 
       c16=c114                                                                 
       c20=c115                                                                 
       c25=c112                                                                 
       s1=s113                                                                  
       s2=-s115                                                                 
       s3=-s112                                                                 
       s4=s111                                                                  
       s5=s114                                                                  
       s6=-s114                                                                 
       s8=s112                                                                  
       s9=s115                                                                  
       s10=-s113                                                                
       s12=s113                                                                 
       s15=s111                                                                 
       s16=s114                                                                 
       s20=s115                                                                 
       s25=-s112                                                                
       go to 11                                                                 
4       c1=c114                                                                 
       c2=c113                                                                  
       c3=c111                                                                  
       c4=c115                                                                  
       c5=c112                                                                  
       c6=c112                                                                  
       c8=c111                                                                  
       c9=c113                                                                  
       c10=c114                                                                 
       c12=c114                                                                 
       c15=c115                                                                 
       c16=c112                                                                 
       c20=c113                                                                 
       c25=c111                                                                 
       s1=s114                                                                  
       s2=-s113                                                                 
       s3=s111                                                                  
       s4=s115                                                                  
       s5=-s112                                                                 
       s6=s112                                                                  
       s8=-s111                                                                 
       s9=s113                                                                  
       s10=-s114                                                                
       s12=s114                                                                 
       s15=s115                                                                 
       s16=-s112                                                                
       s20=s113                                                                 
       s25=s111                                                                 
       go to 11                                                                 
5       c1=c115                                                                 
       c2=c111                                                                  
       c3=c114                                                                  
       c4=c112                                                                  
       c5=c113                                                                  
       c6=c113                                                                  
       c8=c114                                                                  
       c9=c111                                                                  
       c10=c115                                                                 
       c12=c115                                                                 
       c15=c112                                                                 
       c16=c113                                                                 
       c20=c111                                                                 
       c25=c114                                                                 
       s1=s115                                                                  
       s2=-s111                                                                 
       s3=s114                                                                  
       s4=-s112                                                                 
       s5=s113                                                                  
       s6=-s113                                                                 
       s8=-s114                                                                 
       s9=s111                                                                  
       s10=-s115                                                                
       s12=s115                                                                 
       s15=-s112                                                                
       s16=s113                                                                 
       s20=s111                                                                 
       s25=s114                                                                 
       go to 11                                                                 
6       c1=c115                                                                 
       c2=c111                                                                  
       c3=c114                                                                  
       c4=c112                                                                  
       c5=c113                                                                  
       c6=c113                                                                  
       c8=c114                                                                  
       c9=c111                                                                  
       c10=c115                                                                 
       c12=c115                                                                 
       c15=c112                                                                 
       c16=c113                                                                 
       c20=c111                                                                 
       c25=c114                                                                 
       s1=-s115                                                                 
       s2=s111                                                                  
       s3=-s114                                                                 
       s4=s112                                                                  
       s5=-s113                                                                 
       s6=s113                                                                  
       s8=s114                                                                  
       s9=-s111                                                                 
       s10=s115                                                                 
       s12=-s115                                                                
       s15=s112                                                                 
       s16=-s113                                                                
       s20=-s111                                                                
       s25=-s114                                                                
       go to 11                                                                 
7       c1=c114                                                                 
       c2=c113                                                                  
       c3=c111                                                                  
       c4=c115                                                                  
       c5=c112                                                                  
       c6=c112                                                                  
       c8=c111                                                                  
       c9=c113                                                                  
       c10=c114                                                                 
       c12=c114                                                                 
       c15=c115                                                                 
       c16=c112                                                                 
       c20=c113                                                                 
       c25=c111                                                                 
       s1=-s114                                                                 
       s2=s113                                                                  
       s3=-s111                                                                 
       s4=-s115                                                                 
       s5=s112                                                                  
       s6=-s112                                                                 
       s8=s111                                                                  
       s9=-s113                                                                 
       s10=s114                                                                 
       s12=-s114                                                                
       s15=-s115                                                                
       s16=s112                                                                 
       s20=-s113                                                                
       s25=-s111                                                                
       go to 11                                                                 
8       c1=c113                                                                 
       c2=c115                                                                  
       c3=c112                                                                  
       c4=c111                                                                  
       c5=c114                                                                  
       c6=c114                                                                  
       c8=c112                                                                  
       c9=c115                                                                  
       c10=c113                                                                 
       c12=c113                                                                 
       c15=c111                                                                 
       c16=c114                                                                 
       c20=c115                                                                 
       c25=c112                                                                 
       s1=-s113                                                                 
       s2=s115                                                                  
       s3=s112                                                                  
       s4=-s111                                                                 
       s5=-s114                                                                 
       s6=s114                                                                  
       s8=-s112                                                                 
       s9=-s115                                                                 
       s10=s113                                                                 
       s12=-s113                                                                
       s15=-s111                                                                
       s16=-s114                                                                
       s20=-s115                                                                
       s25=s112                                                                 
       go to 11                                                                 
9       c1=c112                                                                 
       c2=c114                                                                  
       c3=c115                                                                  
       c4=c113                                                                  
       c5=c111                                                                  
       c6=c111                                                                  
       c8=c115                                                                  
       c9=c114                                                                  
       c10=c112                                                                 
       c12=c112                                                                 
       c15=c113                                                                 
       c16=c111                                                                 
       c20=c114                                                                 
       c25=c115                                                                 
       s1=-s112                                                                 
       s2=-s114                                                                 
       s3=s115                                                                  
       s4=s113                                                                  
       s5=s111                                                                  
       s6=-s111                                                                 
       s8=-s115                                                                 
       s9=s114                                                                  
       s10=s112                                                                 
       s12=-s112                                                                
       s15=s113                                                                 
       s16=s111                                                                 
       s20=s114                                                                 
       s25=s115                                                                 
       go to 11                                                                 
10       c1=c111                                                                
       c2=c112                                                                  
       c3=c113                                                                  
       c4=c114                                                                  
       c5=c115                                                                  
       c6=c115                                                                  
       c8=c113                                                                  
       c9=c112                                                                  
       c10=c111                                                                 
       c12=c111                                                                 
       c15=c114                                                                 
       c16=c115                                                                 
       c20=c112                                                                 
       c25=c113                                                                 
       s1=-s111                                                                 
       s2=-s112                                                                 
       s3=-s113                                                                 
       s4=-s114                                                                 
       s5=-s115                                                                 
       s6=s115                                                                  
       s8=s113                                                                  
       s9=s112                                                                  
       s10=s111                                                                 
       s12=-s111                                                                
       s15=-s114                                                                
       s16=-s115                                                                
       s20=s112                                                                 
       s25=-s113                                                                
11       do 1110 l=1,m                                                          
       aa=a(iaa)                                                                
       ba=b(iaa)                                                                
       ab=a(ibb)                                                                
       bb=b(ibb)                                                                
       ac=a(icc)                                                                
       bc=b(icc)                                                                
       ad=a(idd)                                                                
       bd=b(idd)                                                                
       ae=a(iee)                                                                
       be=b(iee)                                                                
       af=a(iff)                                                                
       bf=b(iff)                                                                
       ag=a(igg)                                                                
       bg=b(igg)                                                                
        ah=a(ihh)                                                               
       bh=b(ihh)                                                                
       ai=a(iii)                                                                
       bi=b(iii)                                                                
       aj=a(ijj)                                                                
       bj=b(ijj)                                                                
       ak=a(ikk)                                                                
       bk=b(ikk)                                                                
       t1r=ab+ak                                                                
       t1i=bb+bk                                                                
       t2r=ac+aj                                                                
       t2i=bc+bj                                                                
       t3r=ad+ai                                                                
       t3i=bd+bi                                                                
       t4r=ae+ah                                                                
       t4i=be+bh                                                                
       t5r=af+ag                                                                
       t5i=bf+bg                                                                
       t6r=ab-ak                                                                
       t6i=bb-bk                                                                
       t7r=ac-aj                                                                
       t7i=bc-bj                                                                
       t8r=ad-ai                                                                
       t8i=bd-bi                                                                
       t9r=ae-ah                                                                
       t9i=be-bh                                                                
       t10r=af-ag                                                               
       t10i=bf-bg                                                               
       r1r=aa+t1r*c1+t2r*c2+t3r*c3+t4r*c4+t5r*c5                                
       r1i=ba+t1i*c1+t2i*c2+t3i*c3+t4i*c4+t5i*c5                                
       r2r=t6r*s1+t7r*s2+t8r*s3+t9r*s4+t10r*s5                                  
       r2i=t6i*s1+t7i*s2+t8i*s3+t9i*s4+t10i*s5                                  
       r3r=aa+t1r*c2+t2r*c4+t3r*c6+t4r*c8+t5r*c10                               
       r3i=ba+t1i*c2+t2i*c4+t3i*c6+t4i*c8+t5i*c10                               
       r4r=t6r*s2+t7r*s4+t8r*s6+t9r*s8+t10r*s10                                 
       r4i=t6i*s2+t7i*s4+t8i*s6+t9i*s8+t10i*s10                                 
       r5r=aa+t1r*c3+t2r*c6+t3r*c9+t4r*c12+t5r*c15                              
       r5i=ba+t1i*c3+t2i*c6+t3i*c9+t4i*c12+t5i*c15                              
       r6r=t6r*s3+t7r*s6+t8r*s9+t9r*s12+t10r*s15                                
       r6i=t6i*s3+t7i*s6+t8i*s9+t9i*s12+t10i*s15                                
       r7r=aa+t1r*c4+t2r*c8+t3r*c12+t4r*c16+t5r*c20                             
       r7i=ba+t1i*c4+t2i*c8+t3i*c12+t4i*c16+t5i*c20                             
       r8r=t6r*s4+t7r*s8+t8r*s12+t9r*s16+t10r*s20                               
       r8i=t6i*s4+t7i*s8+t8i*s12+t9i*s16+t10i*s20                               
       r9r=aa+t1r*c5+t2r*c10+t3r*c15+t4r*c20+t5r*c25                            
       r9i=ba+t1i*c5+t2i*c10+t3i*c15+t4i*c20+t5i*c25                            
       r10r=t6r*s5+t7r*s10+t8r*s15+t9r*s20+t10r*s25                             
       r10i=t6i*s5+t7i*s10+t8i*s15+t9i*s20+t10i*s25                             
       a(iaa)=aa+t1r+t2r+t3r+t4r+t5r                                            
       b(iaa)=ba+t1i+t2i+t3i+t4i+t5i                                            
       a(ibb)=r1r-r2i                                                           
       b(ibb)=r1i+r2r                                                           
       a(icc)=r3r-r4i                                                           
       b(icc)=r3i+r4r                                                           
       a(idd)=r5r-r6i                                                           
       b(idd)=r5i+r6r                                                           
       a(iee)=r7r-r8i                                                           
       b(iee)=r7i+r8r                                                           
       a(iff)=r9r-r10i                                                          
       b(iff)=r9i+r10r                                                          
       a(igg)=r9r+r10i                                                          
       b(igg)=r9i-r10r                                                          
       a(ihh)=r7r+r8i                                                           
       b(ihh)=r7i-r8r                                                           
       a(iii)=r5r+r6i                                                           
       b(iii)=r5i-r6r                                                           
       a(ijj)=r3r+r4i                                                           
       b(ijj)=r3i-r4r                                                           
       a(ikk)=r1r+r2i                                                           
       b(ikk)=r1i-r2r                                                           
       ix=ik+iskip                                                              
       ik=ij+iskip                                                              
       ij=ii+iskip                                                              
       ii=ih+iskip                                                              
       ih=ig+iskip                                                              
       ig=if+iskip                                                              
       if=ie+iskip                                                              
       ie=id+iskip                                                              
       id=ic+iskip                                                              
       ic=ib+iskip                                                              
       ib=ia+iskip                                                              
       ia=ix                                                                    
       iaa=ia+ia-1                                                              
       ibb=ib+ib-1                                                              
       icc=ic+ic-1                                                              
       idd=id+id-1                                                              
       iee=ie+ie-1                                                              
       iff=if+if-1                                                              
       igg=ig+ig-1                                                              
       ihh=ih+ih-1                                                              
       iii=ii+ii-1                                                              
       ijj=ij+ij-1                                                              
       ikk=ik+ik-1                                                              
1110    continue                                                                
C                                                                               
1000    continue                                                                
       return                                                                   
       end                                                                      
C                                                                               
