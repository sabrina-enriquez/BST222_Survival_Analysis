/*Sabrina Enriquez project 1*/
LIBNAME p1 "/folders/myfolders/p1";
PROC import DATAFILE= "/folders/myfolders/p1/vaLung.csv"
	DBMS=csv  
	 out=p1.vets
		replace;	       
	
	GETNAMES=YES; /*Option GETNAMES determines whether to generate SAS variable names from the data values in the first 
					record of the imported file.*/
RUN;

*descriptive stats;
proc SGPLOT data=p1.vets;
   VBAR trt;
run;

proc SGPLOT data=p1.vets;
   VBAR karno;
run;

proc SGPLOT data=p1.vets;
   VBAR age;
run;

proc template; 
define statgraph simplepie;  
begingraph;    
entrytitle "Treatment Group";  
layout region;      
piechart category=trt / datalabellocation=inside;    
endlayout;  
endgraph; 
end; 
run; 
proc sgrender data=p1.vets           
template=simplepie;
run;

proc template; 
define statgraph simplepie2;  
begingraph;    
entrytitle "Prior Group";  
layout region;      
piechart category=prior / datalabellocation=inside;    
endlayout;  
endgraph; 
end; 
run; 
proc sgrender data=p1.vets           
template=simplepie2;
run;

proc template; 
define statgraph simplepie3;  
begingraph;    
entrytitle "Cell Type Group";  
layout region;      
piechart category=celltype / datalabellocation=inside;    
endlayout;  
endgraph; 
end; 
run; 
proc sgrender data=p1.vets           
template=simplepie3;
run;


proc SGPLOT data=p1.vets;
   VBAR celltype;
run;

proc sgrender data=p1.vets           
template=simplepie;
run;

proc univariate data=p1.vets;
   histogram;
run;


*KM plot;
ods listing gpath= '/folders/myfolders/p1/';
ods graphics / imagename= "lifetest" imagefmt=png;
proc lifetest data=p1.vets outs=p1.KM 
plots=survival(atrisk= 0 to 1000 by 100) ;
time time*status(0); strata trt; run;

*Nelson-Aalen;
ods listing gpath= '/folders/myfolders/p1/';
ods graphics / imagename= "cumhaz" imagefmt=png;
proc phreg data=p1.vets 
plots(overlay)=cumhaz ;
class trt;
model time*status(0)=trt; run;

*Nelson-Aalen;

proc lifetest data=p1.vets nelson method=breslow outs=p1.NA 
plots=logsurv ;
time time*status(0); strata trt; run;




*K- sample test;
proc lifetest data=p1.vets 
plots=(survival) ;
time time*status(0); strata trt / 
test= (logrank tarone peto modpeto fleming(0,1) ) ; 
run;

proc lifetest data=p1.vets 
plots=(hazard) ;
time time*status(0); strata trt / 
test= (logrank tarone peto modpeto fleming(0,1) ) ; 
run;




*now trend test;
proc lifetest data=p1.vets 
plots=survival(atrisk= 0 to 1000 by 100 ) ;
time time*status(0); strata trt / 
trend test= (logrank tarone peto modpeto fleming(0,1) ) ; 
run;



*performing a proportional hazards regression with trt as the single covariate in the model;
proc phreg data =p1.vets 
plots(overlay)=(survival);
class trt;
model time*status(0) = trt;
run;
*not significant

*adjusting for age;
proc phreg data= p1.vets plots(overlay) = (survival);
class trt;
model time*status(0) = trt age ;
run;

*
*
*
*
*
* objective 2 see how other covariates relate to survival;


*cell type ;

*KM plot;
ods listing gpath= '/folders/myfolders/p1/';
ods graphics / imagename= "lifetest" imagefmt=png;
proc lifetest data=p1.vets 
plots=survival(atrisk= 0 to 1000 by 100) ;
time time*status(0); strata celltype; run;


proc lifetest data=p1.vets 
plots=hazard ;
time time*status(0); strata celltype; run;

proc lifetest data=p1.vets 
plots=logsurv ;
time time*status(0); strata celltype; run;

*doesn't pass ph assumption;
*performing a proportional hazards regression with celltype as the single covariate in the model;
proc phreg data=p1.vets plots(overlay) = (survival);
class celltype;
model time*status(0) = celltype;
run;





*adjusting cell type model for age;
proc phreg data= p1.vets plots(overlay) = (survival);
model time*status(0) = age;
strata celltype;
run;
*age not significant but still included;



*adjusting model for karno;
proc phreg data= p1.vets plots(overlay) = (survival);

model time*status(0) = age karno;
strata celltype;
run;
*karno is significant





*looking for karno hazard trend over values;
proc phreg data= p1.vets plots(overlay) = (survival);
class karno;
model time*status(0) =  karno;
run;
*karno is significant




*adjusting cell type model for diagtime;
proc phreg data= p1.vets plots(overlay) = (survival);

model time*status(0) = age karno diagtime;
strata celltype;
run;
*diagtime is not significant;




proc lifetest data=p1.vets
plots=survival(atrisk= 0 to 1000 by 100) ;
time time*status(0); strata prior; run;
*doesnt satisfy ph and is not significant.


*adjusting cell type model for strata prior just in case;
proc phreg data= p1.vets plots(overlay) = (survival);

model time*status(0) = age karno;
strata celltype prior;
run;


*prior is not significant;

proc phreg data= p1.vets plots(overlay) = (survival);
model time*status(0) =prior   ;

run;
*prio




*now we check for interaction terms- let's focus on trt since that was the purpose here; 
proc phreg data=p1.vets;
   class prior celltype(ref='large') trt(desc);
   model time*status(0) = age|trt karno|trt diagtime|trt celltype|trt prior|trt
         / risklimits alpha=0.05;
run;

*final model is:;
proc phreg data= p1.vets plots(overlay) = (survival)  ;

model time*status(0) = age karno/ risklimits alpha=0.05;
strata celltype;
run;

*question 3 let's subset our diagnosis times to 2 groups from median = 5 months;


data p1.diag;
set p1.vets;
  if diagtime>5 then group = "late "; 
  else group = "early";
run;


proc print data=p1.diag; run;


*now trend test;
proc lifetest data=p1.diag 
plots=survival(atrisk= 0 to 1000 by 100 ) ;
time time*status(0); strata group / 
trend test= (logrank tarone peto modpeto fleming(0,1) ) ; 
run;
