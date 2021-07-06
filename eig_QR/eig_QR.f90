module my_module
    
    implicit none
    integer                         ::i,j,dim,k
    real*8,allocatable              ::H(:,:)
    real*8,ALLOCATABLE              ::E(:,:)
    real*8,ALLOCATABLE              ::Hi(:,:)
    real*8,ALLOCATABLE              ::zero(:)
    real*8,ALLOCATABLE              ::Q(:,:),R(:,:)
    real*8                          ::norm
contains
    !创建矩阵
    subroutine createH (n)
        INTEGER                     ::n
        real*8,ALLOCATABLE          ::H1(:,:)
        allocate(H1(n,n))
        do i=1,n
            do j=1,n
                if(j==i-1)then
                    H1(i,j)=1
                elseif(j==i+1)then
                    H1(i,j)=1
                else
                    H1(i,j)=0
                end if
            end do
        end do
        H1(1,n)=1
        H1(n,1)=1
        H=H1
        DEALLOCATE(H1)
        RETURN
    end subroutine
    !创建n维单位矩阵
    subroutine eye (n)
        INTEGER                     ::n
        real*8,ALLOCATABLE          ::Ei(:,:)
        ALLOCATE(Ei(n,n))
        do i = 1, n
            do j = 1, n
                if ( i==j ) then
                    Ei(i,j)=1
                else
                    Ei(i,j)=0
                end if
            end do
        end do
        E=Ei
        DEALLOCATE(Ei)
        return
    end subroutine
    !算向量的模
    subroutine qnorm (x)
        real*8,ALLOCATABLE          ::x(:)
        norm=sqrt(dot_product(x,x))
        return
    end subroutine
    !创建0向量
    subroutine zeros (n)
        INTEGER                     ::n
        real*8,ALLOCATABLE          ::z(:)
        ALLOCATE(z(n))
        do i = 1, n
            z(i)=0.
        end do
        zero=z
        DEALLOCATE(z)
        RETURN
    end subroutine
    !HOuseHolder变换 :)
    subroutine HouseHolder (x)
        real*8,ALLOCATABLE          ::x(:),ei(:),w0(:)
        real*8                      ::w(size(x),1)
        call zeros(size(x))
        ei=zero
        ei(1)=1
        if ( x(1)>0 ) then
            call qnorm(x)
            w0=x+norm*ei
            CALL qnorm(w0)
            w(:,1)=w0/norm
        else
            call qnorm(x)
            w0=x-norm*ei
            CALL qnorm(w0)
            w(:,1)=w0/norm    
        end if
        call eye(size(x))
        Hi=E-2*matmul(w,transpose(w))
        return
    end subroutine

    
    subroutine QRfact (A)
        real*8,ALLOCATABLE          ::A(:,:)
        real*8,ALLOCATABLE          ::D(:,:)
        real*8,ALLOCATABLE          ::x(:)
        real*8                      ::Qi(size(A,1),size(A,1)),Ri(size(A,1),size(A,1))
        INTEGER                     ::c
        
        Ri=A
        do c = 1, size(Ri,1)
            x=Ri(c:,c)
            call HouseHolder(x)
            Ri(c:,c:)=matmul(Hi,Ri(c:,c:))
            call eye(size(Ri,1))
            Qi=E
            
            Qi(c:,c:)=Hi

            if ( c==1 ) then
                Q=Qi
            else
                Q=matmul(Qi,Q)
            end if
        end do
        call eye(size(Ri,1))
        D=E
        do i = 1, size(D,1)
            if ( Ri(i,i)<0 ) then
                D(i,i)=-1
            end if
            
        end do
        
        R=matmul(D,Ri)
        Q=matmul(transpose(Q),D)
        RETURN
    end subroutine

    subroutine eig_QR (A)
        real*8,ALLOCATABLE          ::A(:,:),Ak(:,:),Ak0(:,:)
        LOGICAL                     ::flag=.true.
        real*8                      ::eps=1e-7,sum
        Ak=A
        do while(flag)
            Ak0=Ak
            call QRfact(Ak)
            
            Ak=matmul(R,Q)
            sum=0
            do i = 1, size(Ak,1)
                sum=sum+abs(Ak(i,i)-Ak0(i,i))
            end do
            if ( sum<eps ) then
                flag=.false.
            end if
        end do
        H=Ak

        RETURN
    end subroutine
        
end module

program QM
    use my_module
    implicit none
    ! Variables
    INTEGER                         ::n
    ! Body of QM
    n=9
    call createH(n)
    
    call eig_QR(H)
    print *,H
    pause

end program QM
    
    
    
    
    
    
    
    
