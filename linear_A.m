function [u]=linear_A(nex)

nnx=nex+1;
sk=[]; r1=[]; axpt=[]; k=[];
nop=[]; ph=[]; phd=[]; 

discr();
axb();
u=gelim(sk,r1,nnx);
u=u';

    function discr()
        %***  define values of parameters
        k=1e-3;
        
        %***  x-coordinates
        xfirst=0.;
        xlast=20;
        deltax=(xlast-xfirst)/nex;
        axpt(1)=xfirst;
        
        for i=2:nnx
            axpt(i)=axpt(i-1)+deltax;
        end
        
        %***  nodal numbering ***
        for i=1:nex
            for j=1:2
                nop(i,j)=i+(j-1);
            end
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
        ph(1)=1.-c;
        ph(2)=c;
        phd(1)=-1;
        phd(2)=1;
    end

    function abfind(nell)
                
        for i=1:2
            ngl(i)=nop(nell,i);
        end
        
        wgp = 1.;
        gp  = 0.5;
        tsfun(gp);
        
        x=0.;
        x1=0.;
        for n=1:2
            x=x+axpt(ngl(n))*ph(n);
            x1=x1+axpt(ngl(n))*phd(n);
        end 
        
        for i=1:2
            phx(i)=phd(i)/x1;
        end
        
        for m=1:2
            m1=ngl(m);
            for n=1:2
                n1=ngl(n);
                sk(m1,n1)=sk(m1,n1)                 ...
                    -wgp*x1*0.1*phx(m)*phx(n)          ...
                    -wgp*x1*k*ph(m)*ph(n);
            end
        end
        
        
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
