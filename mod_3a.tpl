// Southern hake model

GLOBALS_SECTION  
  #include <admodel.h> 
  #include <time.h>
  #include <string.h>
  #undef depur
  #undef depuro
  #define depur(object) cout << #object "\n" << object << endl; 
  #define depuro(object) cout << #object "\n" << object << endl; exit(1);


  #undef reporte
  #define reporte(object) report << #object "\n" << object << endl;

  #if defined(WIN32) && !defined(__linux__)
      const char* PLATFORM = "Windows";
  #else
      const char* PLATFORM = "Linux";
  #endif

  adstring BaseFileName;
  adstring ReportFileName;
  adstring ResultsPath;

  adstring stripExtension(adstring fileName)
  {
    /* from Stock-Sintesis */
    const int length = fileName.size();
    for (int i=length; i>=0; --i)
    {
      if (fileName(i)=='.')
      {
        return fileName(1,i-1);
      }
    }
    return fileName;
  }


TOP_OF_MAIN_SECTION  
   arrmblsize = 50000000;
   gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
   gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
   gradient_structure::set_MAX_NVAR_OFFSET(5000);
   gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


DATA_SECTION
  !! cout<<"Modelo de Merluza del Sur corriendo bajo "<<PLATFORM<< endl;

  init_adstring DataFile;
  init_adstring ControlFile;
  init_adstring ResultsFileName;

  !! BaseFileName = stripExtension(ControlFile);
  !! cout << "dat: " << " " << DataFile << endl;
  !! cout << "ctl: " << " " << ControlFile << endl;
  !! cout << "basefileName: " << " " << BaseFileName << endl;
  !! ReportFileName = BaseFileName + adstring(".rep");
  !! ResultsPath = stripExtension(ResultsFileName);


  // |---------------------------------------------------------------------------------|
  // | MODEL DATA FROM DATA FILE
  // |---------------------------------------------------------------------------------|
  // |
  !! ad_comm::change_datafile_name(DataFile);
  init_int nyears;  
  init_int nages;
  init_ivector vyears(1,nyears);
  int styr
  !! styr = vyears(1);
  int styr_pop
  !! styr_pop = vyears(2);
  int endyr
  !! endyr = vyears(nyears);
  init_ivector vages(1,nages);
  int stage
  !! stage = vages(1);
  int endage
  !! endage = vages(nages);
  init_vector surveyindex(styr,endyr);
  init_matrix cpueindex(styr,endyr,1,3);
  vector indxtrawl(styr,endyr)
  vector indxlongline(styr,endyr)
  vector indxartisanal(styr,endyr)
  !! indxtrawl = column(cpueindex,1);
  !! indxlongline = column(cpueindex,2);
  !! indxartisanal = column(cpueindex,3);
  init_matrix landing(styr,endyr,1,3);
  vector ytrawl(styr,endyr)
  vector ylongline(styr,endyr)
  vector yartisanal(styr,endyr)
  !! ytrawl = column(landing,1);
  !! ylongline = column(landing,2);
  !! yartisanal = column(landing,3);
  init_matrix catagetrwl(styr,endyr,stage,endage)
  init_matrix catagelongl(styr,endyr,stage,endage)
  init_matrix catageartisa(styr,endyr,stage,endage)
  init_matrix natagesurvey(styr,endyr,stage,endage)
  init_matrix watagefleet(styr,endyr,stage,endage)
  vector Wm(stage,endage)
  !! Wm = extract_row(watagefleet,endyr)/1000000;
  init_number M
  init_matrix msex(styr,endyr,stage,endage)
  init_number offset

  !! ad_comm::change_datafile_name(ControlFile);
  init_int phs_init
  init_int phs_R
  init_int phs_q
  init_int phs_Sel
  init_int phs_F
  init_number h
  init_int Ptrawl_1
  init_int Ptrawl_2
  int PeT
  !! PeT = Ptrawl_2 - Ptrawl_1 + 1;
  init_int Ppal_1
  init_int Ppal_2
  int PeP
  !! PeP = Ppal_2 - Ppal_1 + 1;
  init_int Pesp_1
  init_int Pesp_2
  int PeE
  !! PeE = Pesp_2 - Pesp_1 + 1;
  init_int chQarr
  init_int inxtrawl
  init_vector strawl(1,inxtrawl)
  init_int inxpal
  init_vector spal(1,inxpal)
  init_int inxesp
  init_vector sesp(1,inxesp)
  init_int inxsurv
  init_vector ssurv(1,inxsurv)
  init_vector nss(1,4)
  init_int Pdarr_1
  init_int Pdarr_2
  init_int Pdpal_1
  init_int Pdpal_2
  init_matrix cv_matrix(styr,endyr,1,4)
  vector cv_arr(styr,endyr)
  vector cv_pal(styr,endyr)
  vector cv_cru(styr,endyr)
  !! cv_arr = column(cv_matrix,2);
  !! cv_pal = column(cv_matrix,3);
  !! cv_cru = column(cv_matrix,4);
  init_vector cv_s(1,4)
  init_number cv_p
  init_vector rango_sa(1,4)
  init_vector rango_sl(1,4)
  init_number cv_sel_a
  init_number cv_sel_c

  init_int yr_sim
  init_int nFt
  init_vector mf(1,nFt)
  init_number offsetCt
  init_number t_selart

INITIALIZATION_SECTION
  log_selA 2.7
  log_selB 1.5
  log_selC 200
  log_qpal -5
  log_qcru 0
  log_qarr -5
  log_Ro 18.5
  mu 0
  log_Farr -1.38
  log_Fesp -2.52
  log_Fpal -1.96


PARAMETER_SECTION
  init_bounded_vector log_selA(1,5,1.5,3.09,phs_Sel)
  init_bounded_vector log_selB(1,5,0.5,2.50,phs_Sel)
  init_bounded_vector log_selC(1,5,1,5.30,phs_Sel)

  init_number log_qpal(phs_q)
  init_number log_qcru(phs_q)
  init_vector log_qarr(1,2,phs_q)

  init_bounded_number log_Ro(17,20,phs_init)
  init_bounded_vector mu(styr_pop,endyr,-1,1,phs_R)
  init_bounded_vector log_Farr(Ptrawl_1,Ptrawl_2,-9,0.6,phs_F)
  init_bounded_vector log_Fesp(Pesp_1,Pesp_2,-9,0.6,phs_F)
  init_bounded_vector log_Fpal(Ppal_1,Ppal_2,-9,0.6,phs_F)
  
  likeprof_number Ro_pl
  sdreport_number Ro
  sdreport_vector SB(styr,endyr) // Biomasa desovante
  sdreport_vector R(styr,endyr) // Reclutamientos
  sdreport_vector BT(styr,endyr)
  sdreport_vector B6(styr,endyr)
  sdreport_vector S_pal(1,nages) // longline selectivity
  sdreport_vector S_arr(1,nages) // trawl selectivity
  sdreport_vector S_esp(1,nages) // artisanal selectivity
  sdreport_vector S_esp_2(1,nages) // artisanal selectivity 2
  sdreport_vector S_cru(1,nages) // survey selectivity
  sdreport_vector Farr(Ptrawl_1,Ptrawl_2)
  sdreport_vector Fpal(Ppal_1,Ppal_2)
  sdreport_vector Fesp(Pesp_1,Pesp_2)
  sdreport_vector Ftot(styr,endyr)
  sdreport_vector Fmax_Total(styr,endyr)
  sdreport_vector muArr(styr,endyr)
  sdreport_vector muPal(styr,endyr)
  sdreport_vector muEsp(styr,endyr)
  sdreport_matrix BDp(endyr+1,endyr+yr_sim,1,nFt)
  sdreport_matrix RPRp(endyr+1,endyr+yr_sim,1,nFt)
  sdreport_vector Reducc(1,nFt)
  sdreport_matrix Yproy(endyr+1,endyr+yr_sim,1,nFt)

   // -----------------------------------------------
  number alpha     // Parameter stock-recruitment relation
  number beta      // Parameter stock-recruitment relation
  number So        // Population matrix
  matrix No(styr,endyr,stage,endage) // population  at age
  matrix NS(styr,endyr,stage,endage) // Spawning abundance at age
  vector NSo(stage,endage)
  vector uno_ages(1,nages)
  vector uno_years(styr,endyr)
  vector uno_years_arr(Ptrawl_1,Ptrawl_2)
  vector uno_years_pal(Ppal_1,Ppal_2)
  vector uno_years_esp(Pesp_1,Pesp_2)
  matrix Fcr_arr(Ptrawl_1,Ptrawl_2,stage,endage)  // Fishing mortality at age trawl 
  matrix Fcr_pal(Ppal_1,Ppal_2,stage,endage)  // Fishing mortality at age longline
  matrix Fcr_esp(Pesp_1,Pesp_2,stage,endage)  // Fishing mortality at age artisanal
  matrix Fcr_total(styr,endyr,stage,endage)   // Fishing mortality at age total
  matrix Z(styr,endyr,stage,endage)   // mortality at age total
  matrix Surv(styr,endyr,stage,endage)   // Survival at age total
  matrix ND(styr,endyr,stage,endage)  // Spawning abundance to middle year
  vector BDcru(styr,endyr)
  matrix Zarr(styr,endyr,stage,endage)  // total mortality at age trawl
  matrix Zpal(styr,endyr,stage,endage)  // total mortality at age longline
  matrix Zesp(styr,endyr,stage,endage)  // total mortality at age artisanal
  vector BMVarr(styr,endyr)
  vector BMVpal(styr,endyr)
  vector BMVesp(styr,endyr)
  vector CPUEarr(styr,endyr)
  vector CPUEpal(styr,endyr)
  vector CPUEesp(styr,endyr)
  vector estBDcru(styr,endyr)
  matrix cageArr(styr,endyr,stage,endage)
  matrix cagePal(styr,endyr,stage,endage)
  matrix cageEsp(styr,endyr,stage,endage)
  matrix NSurvey(styr,endyr,stage,endage)
  vector YestArr(styr,endyr)
  vector YestPal(styr,endyr)
  vector YestEsp(styr,endyr)
  matrix pobsarr(1,inxtrawl,stage,endage)
  matrix pestarr(1,inxtrawl,stage,endage)
  matrix pobspal(1,inxpal,stage,endage)
  matrix pestpal(1,inxpal,stage,endage)
  matrix pobsesp(1,inxesp,stage,endage)
  matrix pestesp(1,inxesp,stage,endage)
  matrix pobssurv(1,inxsurv,stage,endage)
  matrix pestsurv(1,inxsurv,stage,endage)
  vector logL(1,13)
  vector penL(1,12)
  number a
  number sl
  number sr
  vector p(1,2)
  vector d_arr(1,2)
  vector d_pal(1,2)
  vector d_esp(1,2)
  vector d_cru(1,2)
  number log_qesp

  matrix S_total(styr,endyr,stage,endage)   // Selectividad total

  objective_function_value objF

  // proyecciones
  
  number Rp
  number Nplus
  number Yp

  vector Np(stage,endage)
  vector wp(stage,endage) 
  vector msp(stage,endage)
  vector Sp(stage,endage)
  vector Fp(stage,endage)
  vector Zp(stage,endage)
  vector NDp(stage,endage)
  vector Ctp(stage,endage)


PRELIMINARY_CALCS_SECTION
  uno_ages = 1;
  uno_years = 1;
  uno_years_arr = 1;
  uno_years_pal = 1;
  uno_years_esp = 1;
  Ro_pl.set_stepnumber(50);
  Ro_pl.set_stepsize(0.2);

RUNTIME_SECTION
  maximum_function_evaluations 5000, 10000, 100000, 500000
  convergence_criteria 1e-7,1e-8,1e-8, 1e-8

PROCEDURE_SECTION
  selectivity_exploitation_rate();
  initial_age_structure();
  selectivity_penalties();
  dynamics_abundance_per_fleet();
  biomass_and_mortality();
  estimates_cpue_fleet();
  catch_at_age();
  biomass_state();
  evaluate_objective_function();
  
  if(last_phase())
    {
    sim_Fcte();
    }

  if(mceval_phase())
    {
    ofstream out("mod_3.mcmc.out",ios::app);
    out << Ro << " " << objF << " " << Ftot << " " << SB << endl;
    out.close();
    }

FUNCTION selectivity_exploitation_rate
  int t,i;
     
     // Palangre
     a  = mfexp(log_selA(1));
     sl = mfexp(log_selB(1));
     sr = mfexp(log_selC(1));

     for(i=1; i<=nages; i++)
       {
         if(i <= a)
           {S_pal(i) = pow(2,-1*square((i-a)/sl));}
         else
           {S_pal(i) = pow(2,-1*square((i-a)/sr));}
       }
      S_pal = S_pal / (max(S_pal)+1e-6);
  

     // Arrastre
     a  = mfexp(log_selA(2));
     sl = mfexp(log_selB(2));
     sr = mfexp(log_selC(2));

     for(i=1; i<=nages; i++)
       {
         if(i <= a)
           {S_arr(i) = pow(2,-1*square((i-a)/sl));}
         else
           {S_arr(i) = pow(2,-1*square((i-a)/sr));}
       }
      S_arr = S_arr / (max(S_arr)+1e-6);
 

     // Espinel
     a  = mfexp(log_selA(3));
     sl = mfexp(log_selB(3));
     sr = mfexp(log_selC(3));

     for(i=1; i<=nages; i++)
       {
         if(i <= a)
           {S_esp(i) = pow(2,-1*square((i-a)/sl));}
         else
           {S_esp(i) = pow(2,-1*square((i-a)/sr));}
       }
      S_esp = S_esp / (max(S_esp)+1e-6);
 
     // Espinel2
     a  = mfexp(log_selA(4));
     sl = mfexp(log_selB(4));
     sr = mfexp(log_selC(4));

     for(i=1; i<=nages; i++)
       {
         if(i <= a)
           {S_esp_2(i) = pow(2,-1*square((i-a)/sl));}
         else
           {S_esp_2(i) = pow(2,-1*square((i-a)/sr));}
       }
      S_esp_2 = S_esp_2/ (max(S_esp_2)+1e-6);

     // Cruceros
     a  = mfexp(log_selA(5));
     sl = mfexp(log_selB(5));
     sr = mfexp(log_selC(5));

     for(i=1; i<=nages; i++)
       {
         if(i <= a)
           {S_cru(i) = pow(2,-1*square((i-a)/sl));}
         else
           {S_cru(i) = pow(2,-1*square((i-a)/sr));}
       }
      S_cru = S_cru / (max(S_cru)+1e-6);


  Fpal = mfexp(log_Fpal);
  Fesp = mfexp(log_Fesp);
  Farr = mfexp(log_Farr);

  Fcr_pal = elem_prod(outer_prod(uno_years_pal,S_pal),outer_prod(Fpal,uno_ages));
  Fcr_arr = elem_prod(outer_prod(uno_years_arr,S_arr),outer_prod(Farr,uno_ages));
  
  // Artesanal
  for(t=Pesp_1; t<=endyr; t++)
    {
    if(t<t_selart)
    {Fcr_esp(t) = S_esp*Fesp(t);}
    else 
    {Fcr_esp(t) = S_esp_2*Fesp(t);}
    }

  for(t=styr; t<=endyr; t++)
    {
    if(t < Pesp_1)
      { Fcr_total(t) = Fcr_arr(t);}
    else if (t >= Pesp_1 & t < Ppal_1)
      { Fcr_total(t) = Fcr_arr(t) + Fcr_esp(t);}
    else
      {Fcr_total(t) = Fcr_arr(t) + Fcr_esp(t) + Fcr_pal(t);}
    }
  
  Z = Fcr_total + M;
  Surv = mfexp(-1.0 * Z);
  
  for(t=styr; t<=endyr; t++)
    {
    if(t < Pesp_1)
      { Ftot(t) = Farr(t);}
    else if (t >= Pesp_1 & t < Ppal_1)
      { Ftot(t) = Farr(t) + Fesp(t);}
    else
      { Ftot(t) = Farr(t) + Fesp(t) + Fpal(t);}
  
    Fmax_Total(t)= max(Fcr_total(t));
   }



FUNCTION initial_age_structure
  int a;

  Ro = mfexp(log_Ro);
  Ro_pl = Ro;
  No(styr,stage) = Ro;
  for(a=stage+1; a<=endage; a++)
    {
    No(styr,a) = No(styr,a-1)*mfexp(-M);
       if (a==endage)
       {
       No(styr,a) += No(styr,a)/(1.0-mfexp(-M));
       }
    }
  NSo = elem_prod(elem_prod(extract_row(No,styr),msex(styr)),Wm)*mfexp(-M*9/12);
  So = sum(NSo);
  //NS.rowfill(styr,(elem_prod(elem_prod(extract_row(No,styr),msex(styr)),Wm)*mfexp(-M*9/12)));
  //So = sum(extract_row(NS,styr));
  //SB(styr) = So;
  NS.rowfill(styr,(elem_prod(elem_prod(elem_prod(extract_row(No,styr),msex(styr)),Wm),mfexp(-1.0*Z(styr)*9/12))));
  SB(styr) = sum(extract_row(NS,styr));
  alpha = (So/Ro)*(1.0-h)/(4.0*h);
  beta = (5.0*h-1.0)/(4.0*h*Ro);
  R(styr) = SB(styr)/(alpha+(beta*SB(styr))); // *mfexp(mu(styr));



FUNCTION selectivity_penalties
  int t;

     // Palangre
     a  = mfexp(log_selA(1));
     sl = mfexp(log_selB(1));
     sr = mfexp(log_selC(1));
     p(1) = a - sl;
     p(2) = a + sr;

     for(t=1; t<=2; t++)
     {
     if(p(t) <= a)
           {d_pal(t) = pow(2,-1*square((p(t)-a)/sl));}
         else
           {d_pal(t) = pow(2,-1*square((p(t)-a)/sr));}
     }   

     // Arrastre
     a  = mfexp(log_selA(2));
     sl = mfexp(log_selB(2));
     sr = mfexp(log_selC(2));
     p(1) = a - sl;
     p(2) = a + sr;

     for(t=1; t<=2; t++)
     {
     if(p(t) <= a)
           {d_arr(t) = pow(2,-1*square((p(t)-a)/sl));}
         else
           {d_arr(t) = pow(2,-1*square((p(t)-a)/sr));}
     }        

     // Espinel
     a  = mfexp(log_selA(3));
     sl = mfexp(log_selB(3));
     sr = mfexp(log_selC(3));
     p(1) = a - sl;
     p(2) = a + sr;

     for(t=1; t<=2; t++)
     {
     if(p(t) <= a)
           {d_esp(t) = pow(2,-1*square((p(t)-a)/sl));}
         else
           {d_esp(t) = pow(2,-1*square((p(t)-a)/sr));}
     }  
 
     // Cruceros  ---- ojo este es para espinel 2 (CAMBIAR **********************)
     a  = mfexp(log_selA(4));
     sl = mfexp(log_selB(4));
     sr = mfexp(log_selC(4));
     p(1) = a - sl;
     p(2) = a + sr;

     for(t=1; t<=2; t++)
     {
     if(p(t) <= a)
           {d_cru(t) = pow(2,-1*square((p(t)-a)/sl));}
         else
           {d_cru(t) = pow(2,-1*square((p(t)-a)/sr));}
     }   

    

FUNCTION dynamics_abundance_per_fleet
  int t;

  for(t=styr_pop; t<=endyr; t++) // t<=endyr
    {
    R(t) = SB(t-1)/(alpha+(beta*SB(t-1)))*mfexp(mu(t));
    No(t,1) = R(t);
    No(t)(stage+1,endage) =  ++elem_prod(No(t-1)(stage, endage - 1),Surv(t-1)(stage, endage - 1));
    No(t,endage) += No(t,endage)/(1-Surv(t-1,endage)); // plus group
    NS.rowfill(t,(elem_prod(elem_prod(elem_prod(extract_row(No,t),msex(t)),Wm),mfexp(-1.0*Z(t)*9/12))));
    SB(t) =  sum(extract_row(NS,t));    
    }


   
FUNCTION biomass_and_mortality
  int t;

  ND = elem_prod(No,msex)*mfexp(-M*9/12); 
  BDcru = rowsum(elem_prod(elem_prod(ND,outer_prod(uno_years,S_cru)),outer_prod(uno_years,Wm)));

  Zarr = Fcr_arr + M; 
  BMVarr = rowsum(elem_prod(elem_div(1-mfexp(-1.0*Zarr),Zarr),elem_prod(elem_prod(No,outer_prod(uno_years,S_arr)),outer_prod(uno_years,Wm))));
  muArr = elem_div(ytrawl+1e-6,BMVarr);
  
  Zpal = M; for(t=Ppal_1; t<=Ppal_2; t++) {Zpal(t) += Fcr_pal(t);}
  BMVpal = rowsum(elem_prod(elem_div(1-mfexp(-1.0*Zpal),Zpal),elem_prod(elem_prod(No,outer_prod(uno_years,S_pal)),outer_prod(uno_years,Wm))));
  muPal = elem_div(ylongline+1e-6,BMVpal);
  
  Zesp = M; for(t=Pesp_1; t<=Pesp_2; t++) {Zesp(t) += Fcr_esp(t);}

  for(t=styr; t<=endyr; t++) // t<=endyr
  {
  if (t<t_selart)
  {BMVesp(t) = sum( elem_prod( elem_div(1-exp(-1.0*Zesp(t)),Zesp(t)),   elem_prod(  elem_prod(No(t),S_esp) ,Wm)  ));}
  else
  {BMVesp(t) = sum( elem_prod( elem_div(1-exp(-1.0*Zesp(t)),Zesp(t)),   elem_prod( elem_prod(No(t),S_esp_2) ,Wm)  ));}
  }
  muEsp = elem_div(yartisanal+1e-6,BMVesp);



FUNCTION estimates_cpue_fleet

  CPUEarr(styr,chQarr) = mfexp(log_qarr(1))*BMVarr(styr,chQarr); 
  CPUEarr(chQarr+1,endyr) = mfexp(log_qarr(2))*BMVarr(chQarr+1,endyr);

  CPUEpal = mfexp(log_qpal)*BMVpal;

  log_qesp = log_qpal;
  CPUEesp = mfexp(log_qesp)*BMVesp;

  estBDcru = mfexp(log_qcru)*BDcru;



FUNCTION catch_at_age
  double tiny = 1e-6;
  int t;

  cageArr = elem_prod(elem_prod(No,Zarr - M),elem_div(1-mfexp(-1.0*Z),Z));
          YestArr = rowsum(elem_prod(outer_prod(uno_years,Wm),cageArr));

  cagePal = elem_prod(elem_prod(No,Zpal - M),elem_div(1-mfexp(-1.0*Z),Z));
          YestPal = rowsum(elem_prod(outer_prod(uno_years,Wm),cagePal));

  cageEsp = elem_prod(elem_prod(No,Zesp - M),elem_div(1-mfexp(-1.0*Z),Z));
          YestEsp = rowsum(elem_prod(outer_prod(uno_years,Wm),cageEsp));

  NSurvey = elem_prod(elem_prod(No,mfexp(-1.0*Z*9/12)),outer_prod(uno_years,S_cru));

  for(t=1; t<=inxtrawl; t++)
    {
    pobsarr(t) = catagetrwl(strawl(t))/sum(catagetrwl(strawl(t)+tiny)); 
    pestarr(t) = cageArr(strawl(t))/sum(cageArr(strawl(t))); 
    }

  for(t=1; t<=inxpal; t++)
    {
    pobspal(t) = catagelongl(spal(t))/sum(catagelongl(spal(t)+tiny)); 
    pestpal(t) = cagePal(spal(t))/sum(cagePal(spal(t))); 
    }

  for(t=1; t<=inxesp; t++)
    {
    pobsesp(t) = catageartisa(sesp(t))/sum(catageartisa(sesp(t)+tiny)); 
    pestesp(t) = cageEsp(sesp(t))/sum(cageEsp(sesp(t))); 
    }

  for(t=1; t<=inxsurv; t++)
    {
    pobssurv(t) = natagesurvey(ssurv(t))/sum(natagesurvey(ssurv(t)+tiny)); 
    pestsurv(t) = NSurvey(ssurv(t))/sum(NSurvey(ssurv(t))); 
    }



FUNCTION biomass_state
  int t;

  BT = rowsum(elem_prod(No,outer_prod(uno_years,Wm)));
  for (t=styr; t<=endyr; t++)
   {
   B6(t) = sum(elem_prod(No(t)(stage+5, endage),Wm(stage+5, endage)));
   }

  

FUNCTION evaluate_objective_function
  int t;

  logL.initialize();
  
  logL(1) = -1.*nss(1)*sum(elem_prod(pobsarr,log(pestarr)));
  logL(2) = -1.*nss(2)*sum(elem_prod(pobspal,log(pestpal)));
  logL(3) = -1.*nss(3)*sum(elem_prod(pobsesp,log(pestesp)));
  logL(4) = -1.*nss(4)*sum(elem_prod(pobssurv,log(pestsurv)));

   for (t=styr; t<=endyr; t++)
    {
    if (indxtrawl(t)>0)
       {logL(5) += 0.5 * square(log(indxtrawl(t))-log(CPUEarr(t)))/square(cv_arr(t));}
    if (indxlongline(t)>0)
       {logL(6) += 0.5 * square(log(indxlongline(t))-log(CPUEpal(t)))/(square((sqrt(square(cv_pal(t))+square(0.12)))));}
    if (ytrawl(t)>0)
       {logL(7) += 0.5 * square(log(ytrawl(t))-log(YestArr(t)))/square(cv_s(1));}
    if (ylongline(t)>0)
       {logL(8) += 0.5 * square(log(ylongline(t))-log(YestPal(t)))/square(cv_s(2));}
    if (yartisanal(t)>0)
       {logL(9) += 0.5 * square(log(yartisanal(t))-log(YestEsp(t)))/square(cv_s(3));}
    if (surveyindex(t)>0)
       {logL(10) += 0.5 * square(log(surveyindex(t))-log(estBDcru(t)))/square(1.38*cv_cru(t));}
    }

  logL(11) = 0.5 * norm2(mu)/square(cv_s(4));//  + size_count(mu)*log(cvs(7));
  

  // palangre arrastre espinel crucero
  penL(1) = norm2(d_pal - 0.5)/(2*square(cv_p));
  penL(2) = norm2(d_arr - 0.5)/(2*square(cv_p));
  penL(3) = norm2(d_esp - 0.5)/(2*square(cv_p));
  penL(4) = norm2(d_cru - 0.5)/(2*square(cv_p));

  penL(5) = 0.5 * (square(mfexp(log_selA(1))-rango_sa(1))/(2*square(cv_sel_a)));
  penL(6) = 0.5 * (square(mfexp(log_selA(2))-rango_sa(2))/(2*square(cv_sel_a)));
  penL(7) = 0.5 * (square(mfexp(log_selA(3))-rango_sa(3))/(2*square(cv_sel_a)));
  penL(8) = 0.5 * (square(mfexp(log_selA(4))-rango_sa(4))/(2*square(cv_sel_a)));

  penL(9) = 0.5 * (square(mfexp(log_selC(1))-rango_sl(1))/(2*square(10.59*cv_sel_c)));
  penL(10) = 0.5 * (square(mfexp(log_selC(2))-rango_sl(2))/(2*square(10.59*cv_sel_c)));
  penL(11) = 0.5 * (square(mfexp(log_selC(3))-rango_sl(3))/(2*square(cv_sel_c)));
  penL(12) = 0.5 * (square(mfexp(log_selC(4))-rango_sl(4))/(2*square(cv_sel_c)));
                     // depuro(rango_sl);
  //penL(12) += 0.5 * (square(mfexp(log_selC(5))-rango_sl(5))/(2*square(cv_sel_c)));

  objF = sum(logL) + sum(penL); //penL(9) + penL(10) + penL(11) + penL(12); 

 
FUNCTION sim_Fcte
  
  for (int j=1; j<=nFt; j++) 
      {
      Np = No(endyr); 
      Rp = mean(R(endyr-5,endyr));  
      wp = Wm; 
      Sp = Surv(endyr); 
      msp = msex(endyr);

      for (int i=endyr+1; i<=endyr+yr_sim; i++)
      	  {
	  Nplus = 1.0 - Sp(endage); //a utilizar en grupo plus
    	  Np(stage+1,endage) = ++elem_prod(Np(stage,endage-1),Sp(stage,endage-1));
    	  Np(endage) += Np(endage)/Nplus;// Grup plus
	  Np(stage) = Rp;
	  Fp = mf(j)*Fcr_total(endyr)*S_total(endyr);    //0.5	0.6	0.7	0.8	0.9	1
    	  Zp = Fp + M;
    	  Sp = exp(-1.0 * Zp);

	  NDp = elem_prod(elem_prod(Np,exp(-1.0*(9.0/12.0)*Zp)),msp); 

  	  BDp(i,j) = sum(elem_prod(NDp,wp));
          RPRp(i,j)=BDp(i,j)/So;

	  Ctp = elem_prod(elem_div(Fp,Zp),elem_prod(1.0-Sp,Np)); //Baranov
  	  Yp = sum(elem_prod(Ctp,wp));
  	  Yproy(i,j) = Yp;
	  };
       }; 
       	                                          
  Reducc = BDp(endyr+yr_sim)/(SB(endyr)+1e-6);

  
REPORT_SECTION    

 // cout<<So<<endl;exit(1);
  reporte(vages);                                                                                         
  reporte(vyears);
  reporte(PeT);
  reporte(PeP);
  reporte(PeE);                                                           
  reporte(strawl);
  reporte(spal);
  reporte(sesp);
  reporte(ssurv);                                                    
  reporte(logL);                                                                                         
  reporte(penL);
  reporte(S_arr);
  reporte(S_pal);
  reporte(S_esp);
  reporte(S_esp_2);
  reporte(S_cru);
  reporte(surveyindex);
  reporte(indxtrawl);
  reporte(indxlongline);
  reporte(indxartisanal);
  reporte(CPUEarr);
  reporte(CPUEpal);
  reporte(CPUEesp);
  reporte(estBDcru);
  reporte(ytrawl);
  reporte(ylongline);
  reporte(yartisanal);
  reporte(YestArr);
  reporte(YestPal);
  reporte(YestEsp);
  reporte(pobsarr);
  reporte(pestarr);
  reporte(pobspal);
  reporte(pestpal);
  reporte(pobsesp);
  reporte(pestesp);
  reporte(pobssurv);
  reporte(pestsurv);
  reporte(mu);
  reporte(R);
  reporte(SB);
  reporte(BDcru);
  reporte(BMVpal);
  reporte(BMVesp);
  reporte(BMVarr);
  reporte(BT);
  reporte(B6);
  reporte(Farr);
  reporte(Fpal);
  reporte(Fesp);
  reporte(Fcr_total);
  reporte(Yproy);
  reporte(BDp); 
  reporte(Reducc);
  reporte(mf);
  reporte(Ftot);
  reporte(alpha);
  reporte(beta);
  reporte(Wm);
  reporte(No);
  reporte(So);
  reporte(RPRp);
  reporte(S_total(endyr));
  reporte(Fmax_Total);
  reporte(t_selart);

  
FINAL_SECTION
  if(last_phase() && PLATFORM == "Linux")
  {
    adstring cambia = "cp mod_3a.rep " + ResultsPath + ".rep";
    system(cambia);
    
    cambia = "cp mod_3a.par " + ResultsPath + ".par";
    system(cambia);
    
    cambia = "cp mod_3a.std " + ResultsPath + ".std";
    system(cambia);
    
    cambia = "cp mod_3a.cor " + ResultsPath + ".cor";
    system(cambia);
    
    cambia = "cp mod_3a.mcmc.out " + ResultsPath + ".mcmc.out";
    system(cambia);
  }
    
