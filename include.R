if(file.exists("include.so")) file.remove("include.so")
if(file.exists("include.o")) file.remove("include.o")
system("R CMD SHLIB include.c")

if(.Platform$OS.type == "unix"){
  dyn.load("include.so")
} else if(.Platform$OS.type == "windows"){
  dyn.load("include.dll")
} else{
  warning("please check the name of your compiled shared library in this folder and enter it int include.R")
  # HERE
}


# call a c code version of the model
'ThornTimeARU.C' <- 
function(args) 
{
    abu<-rep(0,length=args$steps) 
	res <- .C("ThornTimeARU", #what the function is called in the c code
	
	as.integer(args$steps),  
	as.double(args$initials),  
	as.double(args$TAIR),  
	as.double(args$TSOIL),  
	as.double(args$M),  
	as.double(args$N),  # dummy not used
	as.double(args$FIRE), # dummy not used
    as.double(args$A),   
	as.double(args$Kl),  
	as.double(args$gs),  
	as.double(args$gr),  
	as.double(args$KM),  
	as.double(args$A0),  
	as.double(args$N0),  
	as.double(args$KA),  
	as.double(args$Jc),  
	as.double(args$Jn),  
	as.double(args$q),   
	as.double(args$RHOc), 
	as.double(args$RHOn), 
	as.double(args$Fc),   
	as.double(args$Fn),   
	as.double(args$ma1),   
	as.double(args$ma2),  
	as.double(args$tn1),  
	as.double(args$tn2),  
	as.double(args$mn1),  
	as.double(args$mn2),  
	as.double(args$mn3),  
	as.double(args$mn4),  
	as.double(args$tg1),  
	as.double(args$tg2),  
	as.double(args$tg3),   
	as.double(args$tg4),   
	as.double(args$mg1),   
	as.double(args$mg2),   
	as.double(args$f1),  
	as.double(args$f2),  
	as.double(args$tr1),  
	as.double(args$tr2),  
	as.double(args$uMs), 
	as.double(args$uMr),
	as.double(args$uCs),
	as.double(args$uCr),
	as.double(args$uNs),
	as.double(args$uNr),
	ABUTIME = as.double(abu) )
  
	return(list(abu=res$ABUTIME))	
  
}

# wrapper for calling model
callTTR <- function( steps,
        initials,
        TAIR,
        TSOIL,
        M,
        N,
        FIRE,
        PHOTO, 
        Kl,
        gs,
        gr,
        KM,
        A0, 
        N0, 
        KA,
        Jc,
        Jn,
        q,
        RHOc,
        RHOn,
        Fc,
        Fn,
        ma1,   
        ma2,  
        tn1,  
        tn2,  
        mn1,  
        mn2,  
        mn3,  
        mn4,  
        tg1,  
        tg2,  
        tg3,   
        tg4,   
        mg1,   
        mg2,   
        tr1,  
        tr2,  
        f1,   
        f2,
        uMs,
        uMr,
        uCs,
        uCr,
        uNs,
        uNr)  {
        
        args=list()
        args$steps = steps
        args$initials = initials*0.5
        args$TAIR = TAIR
        args$TSOIL = TSOIL
        args$M = M
        args$N = N
        args$FIRE = FIRE
        args$A = PHOTO 
        args$Kl = Kl
        args$gs = gs
        args$gr = gr
        args$KM = KM
        args$A0 = A0
        args$N0 = N0
        args$KA = KA
        args$Jc = Jc
        args$Jn = Jn
        args$q = q
        args$RHOc = RHOc
        args$RHOn = RHOn
        args$Fc = Fc
        args$Fn = Fn
        args$ma1 = ma1    
        args$ma2 = ma2   
        args$tn1 = tn1
        args$tn2 = tn2  
        args$mn1 = mn1  
        args$mn2 = mn2  
        args$mn3 = mn3   
        args$mn4 = mn4  
        args$tg1 = tg1  
        args$tg2 = tg2  
        args$tg3 = tg3   
        args$tg4 = tg4   
        args$mg1 = mg1   
        args$mg2 = mg2   
        args$tr1 = tr1  
        args$tr2 = tr2  
        args$f1 =  f1
        args$f2 =  f2  
        args$uMs =  uMs  
        args$uMr =  uMr
        args$uCs =  uCs
        args$uCr =  uCr
        args$uNs =  uNs
        args$uNr =  uNr
        abu<-ThornTimeARU.C(args)$abu
        abu[which( is.nan(abu))] <- 1e-08
        abu[which( is.na(abu))] <- 1e-08
        abu[which( !is.finite(abu))] <- 1e-08
        abu[abu==0]  <- 1e-08
        return(abu)
        }        

# call this to run the model from R        
ModelABU <- function(parm, Data)
{
    #--parameters
    parset <- parm[Data$pos.parset]
    sigma <- parm[Data$pos.sigma]
    ou <- parm[Data$pos.ou] # not used in this function
    x0 <- parm[Data$pos.x0]

    
    #--parameters for TTR
    ma1<-parset[1]    
	ma2<-ma1+parset[2]  
	tn1<-parset[3]     
	tn2<-tn1+parset[4]    
	mn1<-parset[5]   
	mn2<-mn1+parset[6]     
	mn3<-mn2+parset[7]     
	mn4<-mn3+parset[8]     
	tg1<-parset[9]     
	tg2<-tg1+parset[10]    
	tg3<-tg2+parset[11]    
	tg4<-tg3+parset[12]    
	mg1<-parset[13]    
	mg2<-mg1+parset[14]     
	tr1<- parset[15]     
	tr2<- tr1 + parset[16]     
	f1<- parset[17]     
	f2<- f1 + parset[18] 
	A0<- parset[19]
	N0<- parset[20]
	uMs<- sigma[1]/1e3
	uMr<- sigma[2]/1e3
	uCs<- sigma[3]/1e3
	uCr<- sigma[4]/1e3
	uNs<- sigma[5]/1e3
	uNr<- sigma[6]/1e3
	
    x <- callTTR(Data$steps,x0,Data$X[,"TAIR"],Data$X[,"TSOIL"],
                    Data$X[,"M"],Data$X[,"N"],Data$X[,"FIRE"],Data$X[,"PHOTO"],
                    Data$Kl,Data$gs,Data$gr, 
                    Data$KM,A0,N0,Data$KA,Data$Jc,Data$Jn,
                    Data$q,Data$RHOc,Data$RHOn,Data$Fc,Data$Fn,
                    ma1,ma2,tn1,tn2,mn1,mn2,mn3,mn4,
                    tg1,tg2,tg3,tg4,mg1,mg2,tr1,tr2,f1,f2,
                    uMs,uMr,uCs,uCr,uNs,uNr) 
    ABU<- x*parset[21]*1000
	return(ABU)
}	
          
