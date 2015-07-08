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
#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <mod_2.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
 cout<<"Modelo de Merluza del Sur corriendo bajo "<<PLATFORM<< endl;
  DataFile.allocate("DataFile");
  ControlFile.allocate("ControlFile");
  ResultsFileName.allocate("ResultsFileName");
 BaseFileName = stripExtension(ControlFile);
 cout << "dat: " << " " << DataFile << endl;
 cout << "ctl: " << " " << ControlFile << endl;
 cout << "basefileName: " << " " << BaseFileName << endl;
 ReportFileName = BaseFileName + adstring(".rep");
 ResultsPath = stripExtension(ResultsFileName);
 ad_comm::change_datafile_name(DataFile);
  nyears.allocate("nyears");
  nages.allocate("nages");
  vyears.allocate(1,nyears,"vyears");
 styr = vyears(1);
 styr_pop = vyears(2);
 endyr = vyears(nyears);
  vages.allocate(1,nages,"vages");
 stage = vages(1);
 endage = vages(nages);
  surveyindex.allocate(styr,endyr,"surveyindex");
  cpueindex.allocate(styr,endyr,1,3,"cpueindex");
  indxtrawl.allocate(styr,endyr);
  indxlongline.allocate(styr,endyr);
  indxartisanal.allocate(styr,endyr);
 indxtrawl = column(cpueindex,1);
 indxlongline = column(cpueindex,2);
 indxartisanal = column(cpueindex,3);
  landing.allocate(styr,endyr,1,3,"landing");
  ytrawl.allocate(styr,endyr);
  ylongline.allocate(styr,endyr);
  yartisanal.allocate(styr,endyr);
 ytrawl = column(landing,1);
 ylongline = column(landing,2);
 yartisanal = column(landing,3);
  catagetrwl.allocate(styr,endyr,stage,endage,"catagetrwl");
  catagelongl.allocate(styr,endyr,stage,endage,"catagelongl");
  catageartisa.allocate(styr,endyr,stage,endage,"catageartisa");
  natagesurvey.allocate(styr,endyr,stage,endage,"natagesurvey");
  watagefleet.allocate(styr,endyr,stage,endage,"watagefleet");
  Wm.allocate(stage,endage);
 Wm = extract_row(watagefleet,endyr)/1000000;
  M.allocate("M");
  msex.allocate(styr,endyr,stage,endage,"msex");
  offset.allocate("offset");
 ad_comm::change_datafile_name(ControlFile);
  phs_init.allocate("phs_init");
  phs_R.allocate("phs_R");
  phs_q.allocate("phs_q");
  phs_Sel.allocate("phs_Sel");
  phs_F.allocate("phs_F");
  h.allocate("h");
  Ptrawl_1.allocate("Ptrawl_1");
  Ptrawl_2.allocate("Ptrawl_2");
 PeT = Ptrawl_2 - Ptrawl_1 + 1;
  Ppal_1.allocate("Ppal_1");
  Ppal_2.allocate("Ppal_2");
 PeP = Ppal_2 - Ppal_1 + 1;
  Pesp_1.allocate("Pesp_1");
  Pesp_2.allocate("Pesp_2");
 PeE = Pesp_2 - Pesp_1 + 1;
  chQarr.allocate("chQarr");
  inxtrawl.allocate("inxtrawl");
  strawl.allocate(1,inxtrawl,"strawl");
  inxpal.allocate("inxpal");
  spal.allocate(1,inxpal,"spal");
  inxesp.allocate("inxesp");
  sesp.allocate(1,inxesp,"sesp");
  inxsurv.allocate("inxsurv");
  ssurv.allocate(1,inxsurv,"ssurv");
  nss.allocate(1,4,"nss");
  Pdarr_1.allocate("Pdarr_1");
  Pdarr_2.allocate("Pdarr_2");
  Pdpal_1.allocate("Pdpal_1");
  Pdpal_2.allocate("Pdpal_2");
  cv_matrix.allocate(styr,endyr,1,4,"cv_matrix");
  cv_arr.allocate(styr,endyr);
  cv_pal.allocate(styr,endyr);
  cv_cru.allocate(styr,endyr);
 cv_arr = column(cv_matrix,2);
 cv_pal = column(cv_matrix,3);
 cv_cru = column(cv_matrix,4);
  cv_s.allocate(1,4,"cv_s");
  cv_p.allocate("cv_p");
  rango_sa.allocate(1,4,"rango_sa");
  rango_sl.allocate(1,4,"rango_sl");
  cv_sel_a.allocate("cv_sel_a");
  cv_sel_c.allocate("cv_sel_c");
  yr_sim.allocate("yr_sim");
  nFt.allocate("nFt");
  mf.allocate(1,nFt,"mf");
  offsetCt.allocate("offsetCt");
}

void model_parameters::initializationfunction(void)
{
  log_selA.set_initial_value(2.7);
  log_selB.set_initial_value(1.5);
  log_selC.set_initial_value(200);
  log_qpal.set_initial_value(-5);
  log_qcru.set_initial_value(0);
  log_qarr.set_initial_value(-5);
  log_Ro.set_initial_value(18.5);
  mu.set_initial_value(0);
  log_Farr.set_initial_value(-1.38);
  log_Fesp.set_initial_value(-2.52);
  log_Fpal.set_initial_value(-1.96);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_selA.allocate(1,4,1.5,3.09,phs_Sel,"log_selA");
  log_selB.allocate(1,4,0.5,2.50,phs_Sel,"log_selB");
  log_selC.allocate(1,4,1,5.30,phs_Sel,"log_selC");
  log_qpal.allocate(phs_q,"log_qpal");
  log_qcru.allocate(phs_q,"log_qcru");
  log_qarr.allocate(1,2,phs_q,"log_qarr");
  log_Ro.allocate(17,20,phs_init,"log_Ro");
  mu.allocate(styr_pop,endyr,-1,1,phs_R,"mu");
  log_Farr.allocate(Ptrawl_1,Ptrawl_2,-9,0.6,phs_F,"log_Farr");
  log_Fesp.allocate(Pesp_1,Pesp_2,-9,0.6,phs_F,"log_Fesp");
  log_Fpal.allocate(Ppal_1,Ppal_2,-9,0.6,phs_F,"log_Fpal");
  Ro_pl.allocate("Ro_pl");
  Ro.allocate("Ro");
  SB.allocate(styr,endyr,"SB");
  R.allocate(styr,endyr,"R");
  BT.allocate(styr,endyr,"BT");
  B6.allocate(styr,endyr,"B6");
  S_pal.allocate(1,nages,"S_pal");
  S_arr.allocate(1,nages,"S_arr");
  S_esp.allocate(1,nages,"S_esp");
  S_cru.allocate(1,nages,"S_cru");
  Farr.allocate(Ptrawl_1,Ptrawl_2,"Farr");
  Fpal.allocate(Ppal_1,Ppal_2,"Fpal");
  Fesp.allocate(Pesp_1,Pesp_2,"Fesp");
  Ftot.allocate(styr,endyr,"Ftot");
  muArr.allocate(styr,endyr,"muArr");
  muPal.allocate(styr,endyr,"muPal");
  muEsp.allocate(styr,endyr,"muEsp");
  BDp.allocate(endyr+1,endyr+yr_sim,1,nFt,"BDp");
  Yproy.allocate(endyr+1,endyr+yr_sim,1,nFt,"Yproy");
  Reducc.allocate(1,nFt,"Reducc");
  Bdepl.allocate(styr,endyr,"Bdepl");
  RPRp.allocate(endyr+1,endyr+yr_sim,1,nFt,"RPRp");
  alpha.allocate("alpha");
  #ifndef NO_AD_INITIALIZE
  alpha.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  So.allocate("So");
  #ifndef NO_AD_INITIALIZE
  So.initialize();
  #endif
  No.allocate(styr,endyr,stage,endage,"No");
  #ifndef NO_AD_INITIALIZE
    No.initialize();
  #endif
  NS.allocate(styr,endyr,stage,endage,"NS");
  #ifndef NO_AD_INITIALIZE
    NS.initialize();
  #endif
  NSo.allocate(stage,endage,"NSo");
  #ifndef NO_AD_INITIALIZE
    NSo.initialize();
  #endif
  uno_ages.allocate(1,nages,"uno_ages");
  #ifndef NO_AD_INITIALIZE
    uno_ages.initialize();
  #endif
  uno_years.allocate(styr,endyr,"uno_years");
  #ifndef NO_AD_INITIALIZE
    uno_years.initialize();
  #endif
  uno_years_arr.allocate(Ptrawl_1,Ptrawl_2,"uno_years_arr");
  #ifndef NO_AD_INITIALIZE
    uno_years_arr.initialize();
  #endif
  uno_years_pal.allocate(Ppal_1,Ppal_2,"uno_years_pal");
  #ifndef NO_AD_INITIALIZE
    uno_years_pal.initialize();
  #endif
  uno_years_esp.allocate(Pesp_1,Pesp_2,"uno_years_esp");
  #ifndef NO_AD_INITIALIZE
    uno_years_esp.initialize();
  #endif
  Fcr_arr.allocate(Ptrawl_1,Ptrawl_2,stage,endage,"Fcr_arr");
  #ifndef NO_AD_INITIALIZE
    Fcr_arr.initialize();
  #endif
  Fcr_pal.allocate(Ppal_1,Ppal_2,stage,endage,"Fcr_pal");
  #ifndef NO_AD_INITIALIZE
    Fcr_pal.initialize();
  #endif
  Fcr_esp.allocate(Pesp_1,Pesp_2,stage,endage,"Fcr_esp");
  #ifndef NO_AD_INITIALIZE
    Fcr_esp.initialize();
  #endif
  Fcr_total.allocate(styr,endyr,stage,endage,"Fcr_total");
  #ifndef NO_AD_INITIALIZE
    Fcr_total.initialize();
  #endif
  Z.allocate(styr,endyr,stage,endage,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  Surv.allocate(styr,endyr,stage,endage,"Surv");
  #ifndef NO_AD_INITIALIZE
    Surv.initialize();
  #endif
  ND.allocate(styr,endyr,stage,endage,"ND");
  #ifndef NO_AD_INITIALIZE
    ND.initialize();
  #endif
  BDcru.allocate(styr,endyr,"BDcru");
  #ifndef NO_AD_INITIALIZE
    BDcru.initialize();
  #endif
  Zarr.allocate(styr,endyr,stage,endage,"Zarr");
  #ifndef NO_AD_INITIALIZE
    Zarr.initialize();
  #endif
  Zpal.allocate(styr,endyr,stage,endage,"Zpal");
  #ifndef NO_AD_INITIALIZE
    Zpal.initialize();
  #endif
  Zesp.allocate(styr,endyr,stage,endage,"Zesp");
  #ifndef NO_AD_INITIALIZE
    Zesp.initialize();
  #endif
  BMVarr.allocate(styr,endyr,"BMVarr");
  #ifndef NO_AD_INITIALIZE
    BMVarr.initialize();
  #endif
  BMVpal.allocate(styr,endyr,"BMVpal");
  #ifndef NO_AD_INITIALIZE
    BMVpal.initialize();
  #endif
  BMVesp.allocate(styr,endyr,"BMVesp");
  #ifndef NO_AD_INITIALIZE
    BMVesp.initialize();
  #endif
  CPUEarr.allocate(styr,endyr,"CPUEarr");
  #ifndef NO_AD_INITIALIZE
    CPUEarr.initialize();
  #endif
  CPUEpal.allocate(styr,endyr,"CPUEpal");
  #ifndef NO_AD_INITIALIZE
    CPUEpal.initialize();
  #endif
  CPUEesp.allocate(styr,endyr,"CPUEesp");
  #ifndef NO_AD_INITIALIZE
    CPUEesp.initialize();
  #endif
  estBDcru.allocate(styr,endyr,"estBDcru");
  #ifndef NO_AD_INITIALIZE
    estBDcru.initialize();
  #endif
  cageArr.allocate(styr,endyr,stage,endage,"cageArr");
  #ifndef NO_AD_INITIALIZE
    cageArr.initialize();
  #endif
  cagePal.allocate(styr,endyr,stage,endage,"cagePal");
  #ifndef NO_AD_INITIALIZE
    cagePal.initialize();
  #endif
  cageEsp.allocate(styr,endyr,stage,endage,"cageEsp");
  #ifndef NO_AD_INITIALIZE
    cageEsp.initialize();
  #endif
  NSurvey.allocate(styr,endyr,stage,endage,"NSurvey");
  #ifndef NO_AD_INITIALIZE
    NSurvey.initialize();
  #endif
  YestArr.allocate(styr,endyr,"YestArr");
  #ifndef NO_AD_INITIALIZE
    YestArr.initialize();
  #endif
  YestPal.allocate(styr,endyr,"YestPal");
  #ifndef NO_AD_INITIALIZE
    YestPal.initialize();
  #endif
  YestEsp.allocate(styr,endyr,"YestEsp");
  #ifndef NO_AD_INITIALIZE
    YestEsp.initialize();
  #endif
  pobsarr.allocate(1,inxtrawl,stage,endage,"pobsarr");
  #ifndef NO_AD_INITIALIZE
    pobsarr.initialize();
  #endif
  pestarr.allocate(1,inxtrawl,stage,endage,"pestarr");
  #ifndef NO_AD_INITIALIZE
    pestarr.initialize();
  #endif
  pobspal.allocate(1,inxpal,stage,endage,"pobspal");
  #ifndef NO_AD_INITIALIZE
    pobspal.initialize();
  #endif
  pestpal.allocate(1,inxpal,stage,endage,"pestpal");
  #ifndef NO_AD_INITIALIZE
    pestpal.initialize();
  #endif
  pobsesp.allocate(1,inxesp,stage,endage,"pobsesp");
  #ifndef NO_AD_INITIALIZE
    pobsesp.initialize();
  #endif
  pestesp.allocate(1,inxesp,stage,endage,"pestesp");
  #ifndef NO_AD_INITIALIZE
    pestesp.initialize();
  #endif
  pobssurv.allocate(1,inxsurv,stage,endage,"pobssurv");
  #ifndef NO_AD_INITIALIZE
    pobssurv.initialize();
  #endif
  pestsurv.allocate(1,inxsurv,stage,endage,"pestsurv");
  #ifndef NO_AD_INITIALIZE
    pestsurv.initialize();
  #endif
  logL.allocate(1,11,"logL");
  #ifndef NO_AD_INITIALIZE
    logL.initialize();
  #endif
  penL.allocate(1,12,"penL");
  #ifndef NO_AD_INITIALIZE
    penL.initialize();
  #endif
  a.allocate("a");
  #ifndef NO_AD_INITIALIZE
  a.initialize();
  #endif
  sl.allocate("sl");
  #ifndef NO_AD_INITIALIZE
  sl.initialize();
  #endif
  sr.allocate("sr");
  #ifndef NO_AD_INITIALIZE
  sr.initialize();
  #endif
  p.allocate(1,2,"p");
  #ifndef NO_AD_INITIALIZE
    p.initialize();
  #endif
  d_arr.allocate(1,2,"d_arr");
  #ifndef NO_AD_INITIALIZE
    d_arr.initialize();
  #endif
  d_pal.allocate(1,2,"d_pal");
  #ifndef NO_AD_INITIALIZE
    d_pal.initialize();
  #endif
  d_esp.allocate(1,2,"d_esp");
  #ifndef NO_AD_INITIALIZE
    d_esp.initialize();
  #endif
  d_cru.allocate(1,2,"d_cru");
  #ifndef NO_AD_INITIALIZE
    d_cru.initialize();
  #endif
  log_qesp.allocate("log_qesp");
  #ifndef NO_AD_INITIALIZE
  log_qesp.initialize();
  #endif
  objF.allocate("objF");
  Rp.allocate("Rp");
  #ifndef NO_AD_INITIALIZE
  Rp.initialize();
  #endif
  Nplus.allocate("Nplus");
  #ifndef NO_AD_INITIALIZE
  Nplus.initialize();
  #endif
  Yp.allocate("Yp");
  #ifndef NO_AD_INITIALIZE
  Yp.initialize();
  #endif
  Np.allocate(stage,endage,"Np");
  #ifndef NO_AD_INITIALIZE
    Np.initialize();
  #endif
  wp.allocate(stage,endage,"wp");
  #ifndef NO_AD_INITIALIZE
    wp.initialize();
  #endif
  msp.allocate(stage,endage,"msp");
  #ifndef NO_AD_INITIALIZE
    msp.initialize();
  #endif
  Sp.allocate(stage,endage,"Sp");
  #ifndef NO_AD_INITIALIZE
    Sp.initialize();
  #endif
  Fp.allocate(stage,endage,"Fp");
  #ifndef NO_AD_INITIALIZE
    Fp.initialize();
  #endif
  Zp.allocate(stage,endage,"Zp");
  #ifndef NO_AD_INITIALIZE
    Zp.initialize();
  #endif
  NDp.allocate(stage,endage,"NDp");
  #ifndef NO_AD_INITIALIZE
    NDp.initialize();
  #endif
  Ctp.allocate(stage,endage,"Ctp");
  #ifndef NO_AD_INITIALIZE
    Ctp.initialize();
  #endif
  Fproy.allocate(endyr+1,endyr+yr_sim,1,nFt,"Fproy");
  #ifndef NO_AD_INITIALIZE
    Fproy.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
  uno_ages = 1;
  uno_years = 1;
  uno_years_arr = 1;
  uno_years_pal = 1;
  uno_years_esp = 1;
  Ro_pl.set_stepnumber(50);
  Ro_pl.set_stepsize(0.2);
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{5000, 10000, 100000, 500000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-7,1e-8,1e-8, 1e-8}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::userfunction(void)
{
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
    ofstream out("mod_2.mcmc.out",ios::app);
    out << Ro << " " << objF << " " << Ftot << " " << SB << endl;
    out.close();
    }
}

void model_parameters::selectivity_exploitation_rate(void)
{
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
     // Cruceros
     a  = mfexp(log_selA(4));
     sl = mfexp(log_selB(4));
     sr = mfexp(log_selC(4));
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
  Fcr_esp = elem_prod(outer_prod(uno_years_esp,S_esp),outer_prod(Fesp,uno_ages));
  Fcr_arr = elem_prod(outer_prod(uno_years_arr,S_arr),outer_prod(Farr,uno_ages));
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
    }
}

void model_parameters::initial_age_structure(void)
{
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
}

void model_parameters::selectivity_penalties(void)
{
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
     // Cruceros
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
}

void model_parameters::dynamics_abundance_per_fleet(void)
{
  int t;
  for(t=styr_pop; t<=endyr; t++) // t<=endyr
    {
    R(t) = SB(t-1)/(alpha+(beta*SB(t-1)))*mfexp(mu(t));
    No(t,1) = R(t);
    No(t)(stage+1,endage) =  ++elem_prod(No(t-1)(stage, endage - 1),Surv(t-1)(stage, endage - 1));
    No(t,endage) += No(t,endage)/(1-Surv(t-1,endage)); // plus group
    // No(t,endage) = No(t,endage) + No(t-1,endage)*Surv(t-1,endage); // plus group
    NS.rowfill(t,(elem_prod(elem_prod(elem_prod(extract_row(No,t),msex(t)),Wm),mfexp(-1.0*Z(t)*9/12))));
    SB(t) =  sum(extract_row(NS,t));    
    }
}

void model_parameters::biomass_and_mortality(void)
{
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
  BMVesp = rowsum(elem_prod(elem_div(1-exp(-1.0*Zesp),Zesp),elem_prod(elem_prod(No,outer_prod(uno_years,S_esp)),outer_prod(uno_years,Wm))));
  muEsp = elem_div(yartisanal+1e-6,BMVesp);
}

void model_parameters::estimates_cpue_fleet(void)
{
  CPUEarr(styr,chQarr) = mfexp(log_qarr(1))*BMVarr(styr,chQarr); 
  CPUEarr(chQarr+1,endyr) = mfexp(log_qarr(2))*BMVarr(chQarr+1,endyr);
  CPUEpal = mfexp(log_qpal)*BMVpal;
  log_qesp = log_qpal;
  CPUEesp = mfexp(log_qesp)*BMVesp;
  estBDcru = mfexp(log_qcru)*BDcru;
}

void model_parameters::catch_at_age(void)
{
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
}

void model_parameters::biomass_state(void)
{
  int t;
  BT = rowsum(elem_prod(No,outer_prod(uno_years,Wm)));
  for (t=styr; t<=endyr; t++)
   {
   B6(t) = sum(elem_prod(No(t)(stage+5, endage),Wm(stage+5, endage)));
   Bdepl(t)=SB(t)/So;
  }
}

void model_parameters::evaluate_objective_function(void)
{
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
  objF = sum(logL) + penL(9) + penL(10) + penL(11) + penL(12); 
}

void model_parameters::sim_Fcte(void)
{
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
			if (i==endyr+1)
			{
				Fp = 0.85*Fcr_total(endyr);
			}
			else
			{
				Fp = mf(j)*0.24*(Fcr_total(endyr)/max(Fcr_total(endyr))); //Fcr_total(endyr);
			}
			Zp = Fp + M;
			Sp = exp(-1.0 * Zp);
			NDp 		= elem_prod(elem_prod(Np,exp(-1.0*(9.0/12.0)*Zp)),msp);
			BDp(i,j)	= sum(elem_prod(NDp,wp));
			RPRp(i,j)	= BDp(i,j)/So;
			Ctp			= elem_prod(elem_div(Fp,Zp),elem_prod(1.0-Sp,Np)); //Baranov
			Yp 			= sum(elem_prod(Ctp,wp));
			Yproy(i,j) 	= Yp;
			Fproy(i,j)	= max(Fp); 
		};
	}; 
  Reducc = BDp(endyr+yr_sim)/(SB(endyr)+1e-6);
    if(mceval_phase())
    {
    ofstream pry("proyecciones.mcmc.out",ios::app);
    for (int i=endyr+1; i<=endyr+yr_sim; i++)
    {
      pry << "Captura " << i << Yproy(i) << endl;
    }
    for (int i=endyr+1; i<=endyr+yr_sim; i++)
    {
      pry << "BD " << i << BDp(i) << endl;
    } 
    for (int i=endyr+1; i<=endyr+yr_sim; i++)
    {
      pry << "depletion " << i << RPRp(i) << endl;
    }
    pry.close();
    }
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
  reporte(Fproy);
  reporte(Fp);
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
   arrmblsize = 50000000;
   gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
   gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
   gradient_structure::set_MAX_NVAR_OFFSET(5000);
   gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  #if defined(__GNUDOS__) || defined(DOS386) || defined(__DPMI32__)  || \
     defined(__MSVC32__)
      if (!arrmblsize) arrmblsize=150000;
  #else
      if (!arrmblsize) arrmblsize=25000;
  #endif
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
