function [u]=quad_A(nex)

nnx=2*nex+1;
sk=[]; r1=[]; axpt=[]; k=[];
nop=[]; w=[]; gp=[]; ph=[]; phd=[];

discr();
axb();
u=gelim(sk,r1,nnx);
u=u';

    function discr()
        %***  define values of parameters
        k=1e-3;
        
        %***  x-coordinates
        xfirst=0.;
        xlast=20.;
        deltax=(xlast-xfirst)/nex;
        axpt(1)=xfirst;
        
        for i=2:nnx
            axpt(i)=axpt(i-1)+deltax/2.;
        end
        
        %***  nodal numbering ***
        for i=1:nex
            nop(i,1)=1+2*(i-1);
            nop(i,2)=nop(i,1)+1;
            nop(i,3)=nop(i,2)+1;
        end
        %*************************
    end

    function axb()
        for i=1:nnx
            r1(i)=0.;
            for j=1:nnx
                sk(i,j)=0.;
            end
        end
        
        for nell=1:nex
            abfind(nell);
        end
        %***  BOUNDARY CONDITIONS
        
        %left bc
        for i=1:nnx
            sk(1,i)=0.;
        end

        sk(1,1)=1.;
        r1(1)=1.5;

        %right bc
        for i=1:nnx
            sk(nnx,i)=0.;
        end
        
        sk(nnx,nnx)=1.;
        r1(nnx)=0;
    end

    function tsfun(x)
        c=x;
        ph(1)=2.*c^2-3.*c+1.;
        ph(2)=-4.*c^2+4.*c;
        ph(3)=2.*c^2-c;
        phd(1)=4.*c-3.;
        phd(2)=-8.*c+4.;
        phd(3)=4.*c-1.;
    end

    function abfind(nell)
        w  = [0.27777777777778, 0.444444444444, 0.27777777777778];
        gp = [0.1127016654    , 0.5           , 0.8872983346    ];
        
        for i=1:3
            ngl(i)=nop(nell,i);
        end
        
        for j=1:3 %Loop over qauss points
            
            tsfun(gp(j));
            
            x=0.;
            x1=0.;
            for n=1:3
                x=x+axpt(ngl(n))*ph(n);
                x1=x1+axpt(ngl(n))*phd(n);
            end
            
         
            for i=1:3
                phx(i)=phd(i)/x1;
            end
            
            for m=1:3
                m1=ngl(m);
                for n=1:3
                    n1=ngl(n);
                    sk(m1,n1)=sk(m1,n1)                 ...
                        -w(j)*x1*0.1*phx(m)*phx(n)          ...
                        -w(j)*x1*k*ph(m)*ph(n);
                end
            end
            
        end %End of loop over gauss points
    end

    function [x]=gelim(a,b,n)
        %step 1: forward elimination
        for k=1:n-1
            for i=k+1:n
                c=a(i,k)/a(k,k);
                a(i,k) = 0.0;
                b(i)=b(i)- c*b(k);
                for j=k+1:n
                    a(i,j) = a(i,j)-c*a(k,j);
                end
            end
        end
        
        %step 2: back substitution
        x(n) = b(n)/a(n,n);
        for i=n-1:-1:1
            c=0.0;
            for j=i+1:n
                c= c + a(i,j)*x(j);
            end
            x(i) = (b(i)- c)/a(i,i);
        end
    end
end
