/*Sabrina Enriquez project 2*/
LIBNAME p2 "/folders/myfolders/p2";


PROC import DATAFILE= "/folders/myfolders/p2/finalsimp2.csv"
	DBMS=csv  
	 out=p2.sim
		replace;	       
	
	GETNAMES=YES; /*Option GETNAMES determines whether to generate SAS variable names from the data values in the first 
					record of the imported file.*/
RUN;

data p2.years;
set p2.sim;
keep Diagtime;
run;
PROC FORMAT;
  VALUE  Diagnosistime 1="1974-1983"
  					2 ="1984-1993"
					3  ="1994-2003"
						4="2003-2013";
  
RUN;

PROC FORMAT;
  VALUE  therapy 1="(R): Radiotherapy Only"
  					2 ="(S): Surgery Only"
					3  ="(RS): Radiotherapy and Surgery"
						4="(NRS): Neither R nor S";
  
RUN;
*fig1;
/* ods listing gpath= '/folders/myfolders/p2/plotsandfigs'; */
/* ods graphics / imagename= "fig1" imagefmt=png; */
title 'Incidence Pattern';
PROC SGPLOT DATA=p2.years;
   FORMAT  Diagtime Diagnosistime.;
   Vbar Diagtime / datalabel;
      label Diagtime= 'Years of Diagnosis';     
RUN;


*fig2;
/* ods listing gpath= '/folders/myfolders/p2/plotsandfigs'; */
/* ods graphics / imagename= "fig2" imagefmt=png; */
proc lifetest data=p2.sim 
plots=survival(strata=panel) ;
FORMAT  Diagtime Diagnosistime.;
*FORMAT  Therapy therapy.;
time time*censor(1); strata Therapy/ group= Diagtime; run;

/* ODS TAGSETS.csv */
/* file= '/folders/myfolders/p2/plotsandfigs/therapycox.csv' */
/* style= minimal; */
proc phreg
data= p2.sim
plots(overlay) = survival;
FORMAT  Therapy therapy.;
class Therapy;
model time*censor(1) = Therapy;
run;
/* ods tagsets.csv close; */

*fig3;
*cancer specific to classic;
data p2.cspec;
set p2.sim;
if classic=1;
run;


/* ods listing gpath= '/folders/myfolders/p2/plotsandfigs'; */
/* ods graphics / imagename= "fig3" imagefmt=png; */
proc lifetest data=p2.cspec 
plots=survival(strata=panel) ;
FORMAT  Diagtime Diagnosistime.;
time time*censor(1); strata Therapy/ group= Diagtime; run;

*fig4a
*overall survival by tumor location;
*we need to take the subset of data with diagtime ==4;
data p2.d4;
set p2.sim;
if DiagTime=4;
run;

data p2.cspecd4;
set p2.d4;
if classic=1;
run;



/* ods listing gpath= '/folders/myfolders/p2/plotsandfigs'; */
/* ods graphics / imagename= "fig4a" imagefmt=png; */
title 'Overall Survival by tumor location';
proc lifetest data=p2.d4
plots=survival ;
time time*censor(1); 
strata TumorLocation; 
label TumorLocation= 'Tumor location: 1= mobile spine, 2= sacrum';
run;
*fig4b
*cancer specific survival for classic;
/* ods listing gpath= '/folders/myfolders/p2/plotsandfigs'; */
/* ods graphics / imagename= "fig4b" imagefmt=png; */
proc lifetest data=p2.cspecd4 
plots=survival ;
time time*censor(1); 
strata TumorLocation;
title 'Cancer Specific Survival';
label TumorLocation= 'Tumor location: 1= mobile spine, 2= sacrum';
run;


PROC FORMAT;
  VALUE  gender 1="Male"
  					0 ="Female";
	Value marriage 1="Married"
  					2 ="Never Married"
  					3= "Widowed"
  					4="Divorced"
  					5="Other";
  	Value race 1="White"
  					2 ="Asian/ Pacific Islander"
  					3= "Black"
  					4="Other";
  	Value loc 1="Mobile spine"
  			2= "Sacrum";
  	Value path 1="Classic"
  				2="Chondroid"
  				3="Dedifferentiated";
  					
  
RUN;


/* ODS TAGSETS.csv */
/* file= '/folders/myfolders/p2/plotsandfigs/table1.csv' */
/* style= minimal; */
proc means
data=p2.sim;
var Age TumorSize;
run;

PROC FREQ DATA=p2.sim;
   FORMAT  Gender gender.
           Marital    marriage.
           Race race.
           TumorLocation loc.
           Pathology path.;
   TABLES Gender Marital Race TumorLocation Pathology;
RUN;
/* ods tagsets.csv close; */



*cancer specific to chondroid;
data p2.chon;
set p2.sim;
if chon=1;
run;

proc lifetest data=p2.chon 
plots=survival(strata=panel) ;
time time*censor(1); strata Therapy/ group= Diagtime; run;

*for this I need dummy variables for married vs others, white race vs others, 
surgery vs others, rad vs others, RS vs others, classic vs dediff, chondroid vs dediff. 
We made all of these except for white race and married so we do that now and add 2 columns.;

*married=1 vs others=0;
data p2.married;
set p2.d4;
if Marital=1 then married = 1;
else married = 0;
run;


*white=1 vs others=0;
data p2.white;
set p2.married;
if Race=1 then white = 1;
else white = 0;
run;

*now i want to produce the hr for table 2;
/* ODS TAGSETS.csv */
/* file= '/folders/myfolders/p2/plotsandfigs/table2.csv' */
/* style= minimal; */
proc phreg data= p2.white plots(overlay) = (survival);
FORMAT  Gender gender.
           TumorLocation loc.
           Therapy therapy.
           Pathology path.;
class married(desc) Gender(desc) white(desc) TumorLocation(desc) surgery(desc) rad(desc) RS(desc)  Pathology;
model time*censor(1) = Age TumorSize married Gender white TumorLocation surgery rad RS Pathology/
Ties= EXACT RISKLIMITS ALPHA=.05;
run;
/* ods tagsets.csv close; */

*now the multivariate model;
/* ODS TAGSETS.csv */
/* file= '/folders/myfolders/p2/plotsandfigs/table2b.csv' */
/* style= minimal; */
proc phreg data= p2.white plots(overlay) = (survival);
FORMAT     TumorLocation loc.;
class  TumorLocation(desc) surgery(desc) ;
model time*censor(1) = Age TumorLocation surgery /
Ties= EXACT RISKLIMITS ALPHA=.05;
run;
/* ods tagsets.csv close; */


*let's make table 3;
*cancer specific to classic;
data p2.cspecd4;
set p2.white;
if classic=1 ;
run;





*cancer specific to dediff;
data p2.dediff;
set p2.white;
if Pathology=3;
run;
*cancer specific to classic;
data p2.chon2;
set p2.white;
if chon=1 ;
run;

*dediff and classic;
data p2.dc;
set p2.cspecd4 p2.dediff;

*chon and classic;
data p2.cc;
set p2.chon2 p2.dediff;
run;



*generate table 3;
/* ODS TAGSETS.csv */
/* file= '/folders/myfolders/p2/plotsandfigs/table3.csv'; */
proc phreg data= p2.cspecd4 ;
FORMAT  Gender gender.
           TumorLocation loc.
           
           Pathology path.;
class married(desc) Gender(desc) white(desc) TumorLocation(desc) surgery(desc) rad(desc) RS(desc) ;
model time*censor(1) = Age TumorSize married Gender white TumorLocation surgery rad RS /
Ties= EXACT RISKLIMITS ALPHA=.05;run;

proc phreg data= p2.dc ;

class classic(desc) ;
model time*censor(1) = classic /
Ties= EXACT RISKLIMITS ALPHA=.05;

proc phreg data= p2.cc plots(overlay) = (survival);

class chon(desc) ;
model time*censor(1) = chon /
Ties= EXACT RISKLIMITS ALPHA=.05;

run;
/* ods tagsets.csv close;  */

*multivariate;
/* ODS TAGSETS.csv */
/* file= '/folders/myfolders/p2/plotsandfigs/table3b.csv'; */
proc phreg data= p2.cspecd4 ;
FORMAT  Gender gender.
           TumorLocation loc.
           
           Pathology path.;
class TumorLocation(desc) surgery(desc);
model time*censor(1) = Age TumorLocation surgery /
Ties= EXACT RISKLIMITS ALPHA=.05;run;
/* ods tagsets.csv close;  */
*that's it for regenerating the results;







*now we do trend tests for the full set;

*K- sample test;
/* ods listing gpath= '/folders/myfolders/p2/plotsandfigs'; */
/* ods graphics / imagename= "t1" imagefmt=png; */
ODS TAGSETS.csv
file= '/folders/myfolders/p2/plotsandfigs/t1.csv';
title 'KM survival curves by treatment type for all patients';
proc lifetest data=p2.sim
plots=(survival (atrisk= 0 to 765 by 100)) ;
Format Therapy therapy.;
time time*censor(1); strata Therapy / 
test= (logrank tarone peto modpeto fleming(0,1) ) ; 
run;



*now trend test;
ods listing gpath= '/folders/myfolders/p2/plotsandfigs';
ods graphics / imagename= "t2" imagefmt=png;
proc lifetest data=p2.sim
plots=survival(atrisk= 0 to 765 by 100) ;
Format Therapy therapy.;
time time*censor(1); strata Therapy /  
trend test= (logrank tarone peto modpeto fleming(0,1) ) ; 
run;
/* ods tagsets.csv close; */

*now do opposite of their plots;
ods listing gpath= '/folders/myfolders/p2/plotsandfigs';
ods graphics / imagename= "t3" imagefmt=png;
proc lifetest data=p2.sim 
plots=survival(strata=panel) ;
FORMAT  Diagtime Diagnosistime.;
FORMAT  Therapy therapy.;
time time*censor(1); strata Diagtime/ group= Therapy; run;


ods listing gpath= '/folders/myfolders/p2/plotsandfigs';
ods graphics / imagename= "t4" imagefmt=png;
proc lifetest data=p2.cspec 
plots=survival(strata=panel) ;
FORMAT  Diagtime Diagnosistime.;
FORMAT  Therapy therapy.;
time time*censor(1); strata Diagtime/ group= Therapy; run;







