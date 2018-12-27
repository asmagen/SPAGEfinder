#include <omp.h>
#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>		/* constants */
#define MIN(a,b) a<b?a:b 
#define MAX(a,b) a>b?a:b 
#define ABS(a) (a>0.0) ? a:-a
// #include <iostream>
// #include <fstream>

using namespace std;		// shorthand
using namespace arma;		// shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat logRankMultiClass(arma::mat times, arma::uvec cls, arma::urowvec classMap, int numClass){
    // cout << 25<< endl;
    // cls should start from 1 
    // events: 0 - death; 1 - censored
    //remove rows when classMap = 0

        arma::ucolvec cls_grp = (classMap(cls -1));
        uvec sel = arma::find(cls_grp < 1);
        if(sel.size() > 0) {
            sel = arma::find(cls_grp > 0);
            times = times.rows(sel);
            cls = cls(sel);
        }

    int tLen = times.n_rows;
    // cout << 25<< endl;

    mat out(numClass, 5);
    out.fill(-1000);
    if (tLen<1)
    {
	return out;
    }
    // cout << 25<< endl;
    double maxTime = max(times.col(0));
    uvec timesSort = sort_index(times.col(0));
    // cout << 25<< endl;
    // mat tab(tLen,3 * (numClass + 1));
    // tab.zeros();
    mat eventMat(tLen, numClass + 1);
    mat deathMat(tLen, numClass + 1);
    vec timePoint(tLen);
    eventMat.zeros(); deathMat.zeros(); timePoint.zeros();
    double timeq, timecurr, eventcurr;
    int timecnt;
    int clscurr;
    // cout << trans(eventMat(span(0,10), 0)) << endl; 
    // cout << trans(deathMat(span(0,10), 0)) << endl; 

    timecnt =0; timeq = -100000;
    for (int tt = 0; tt < times.n_rows; tt++)
    {
	timecurr = times(timesSort(tt),0);
	// cout << 45<< endl;
	eventcurr = times(timesSort(tt),1);
	clscurr = classMap(cls(timesSort(tt)) -1); /// this should start with 1 because 0 corresponds to all groups
	if ((tt ==0 )|| (timeq != timecurr))
	{
	    timeq = timecurr;
	    timecnt++;
	    timePoint(timecnt-1) = timeq; 
	}
	eventMat(timecnt-1,0) = eventMat(timecnt-1,0) +1; 
	eventMat(timecnt-1,clscurr) =  eventMat(timecnt-1,clscurr) +1; 
	if (eventcurr==0){
	    deathMat(timecnt-1,0) =  deathMat(timecnt-1,0) +1; 
	    deathMat(timecnt-1,clscurr) =  deathMat(timecnt-1,clscurr) +1; 
	}
	// cout << 57<< endl;
	// cout << "tt" << tt << endl;
    } 
    // cout << tab.rows(span(0,10)) << endl;
    // cout << 61<< endl;
    // cout << "event " << trans(eventMat(span(0,10), 0)) << endl; 
    // cout << "death " << trans(deathMat(span(0,10), 0)) << endl; 
    // cout << "timePoint " << trans(timePoint(span(0,10))) << endl; 
    eventMat.resize(timecnt,numClass +1);
    deathMat.resize(timecnt,numClass +1);
    timePoint.resize(timecnt);
    mat survivalMat(timecnt, numClass + 1);
    vec chi_stat(numClass), p(numClass);
    // cout << 63<< endl;

    int tabSize = eventMat.n_rows;
    if(tabSize < 2) return out;
    int npop;
    npop = sum(eventMat.col(0));
    survivalMat(0,0) = npop;
    survivalMat(span(1,tabSize-1), 0) = npop - cumsum(eventMat(span(0, tabSize -2),0));
    vec rate = deathMat.col(0)/survivalMat.col(0);
    // cout <<  "survival" << trans(survivalMat(span(0,10), 0)) << endl; 
    // cout << "rate" << trans(rate(span(0,10))) << endl; 
    vec Exp(numClass), Obs(numClass);
    mat temp(1,1);
    for (int cc = 1; cc < numClass +1; cc++)
    {
	npop = sum(eventMat.col(cc));
	survivalMat(0,cc) = npop;
	survivalMat(span(1,tabSize-1), cc) = npop - cumsum(eventMat(span(0, tabSize -2),cc));
	temp = trans(rate )* survivalMat.col(cc);//% expected number of deaths group 1
	Exp(cc-1) =  temp(0,0);
	// cout << temp<< endl;
	Obs(cc-1) = sum(deathMat.col(cc));
	// We now calculate the significance of this result, using Chi distribution
	chi_stat(cc-1)=((Exp(cc-1)-Obs(cc-1))*(Exp(cc-1)-Obs(cc-1)))/Exp(cc-1); 
	p(cc-1) = 1- R::pchisq(chi_stat(cc-1),1, 1, 0); 
    }

    // cout << temp.n_rows << temp.n_cols<< endl;

    vec survivalRate(numClass), deltaTime(numClass), auc(numClass);
    survivalRate.ones(); deltaTime.zeros();auc.zeros();
    // double startTime = 0; // min(lastDeath1, lastDeath2);
    // lastDeath1 = lastDeath2  = startTime;
    vec lastDeath(numClass); lastDeath.fill(timePoint(0) );

    double endTime = timePoint(tabSize -1);
    for (int tt = 0; tt < tabSize; tt++){
	timecurr = timePoint(tt);
	// cout << "tt" << tt << endl;
	for (int cc = 1; cc < numClass +1; cc++){

	    // cout << cc << "cc" << tt << endl;
	    if((deathMat(tt,cc) >0) || (timecurr == endTime) ){
		deltaTime(cc-1) = timecurr - lastDeath(cc-1);
		auc(cc-1) += deltaTime(cc-1) * survivalRate(cc-1);
		survivalRate(cc -1) = survivalRate(cc-1) *( 1 -  (deathMat(tt,cc)/survivalMat(tt,cc))); 
		lastDeath(cc-1) = timecurr; 
	    }
	}
	if(timecurr >=endTime)
	break;
    }
    auc  = auc/( endTime) ; 
    // out.col(0) = p;

    out.col(0) = p;
    out.col(1) = Obs;
    out.col(2) = Exp;
    out.col(3) = chi_stat;
    out.col(4)=auc;
    uvec nan_inx=  find(p != p );
    if(nan_inx.size() >0){
	// mat init(nan_inx.size(),1);
	// init.fill(-1000);
	// p(nan_inx).fill( -1000);
	out.rows(nan_inx).fill(-1000);
    }

    return out;
}

using namespace std;            // shorthand
using namespace arma;           // shorthand
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List aggregateLogRankPairsMultiClass(arma::umat pairs, arma::umat scna1Mat, arma::umat scna2Mat, arma::mat survival, arma::vec typeInx, arma::urowvec classMap, int typeNum=1, int threads=1, int numClass = 9){
    #ifdef _OPENMP
    if ( threads > 0 )
    omp_set_num_threads( threads );
    #endif
    int pairNum = pairs.n_rows;
    int numSamples = scna1Mat.n_cols;
    mat ExpMat(pairNum, numClass), ObsMat(pairNum, numClass), aucMat(pairNum, numClass), pMat(pairNum,numClass);
    vec p_agg(pairNum);
    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < pairNum; idx++)
    {

	urowvec scna1;
	scna1 = scna1Mat.row(pairs(idx,0)-1);
	urowvec scna2;
	scna2 = scna2Mat.row(pairs(idx,1)-1);
	vec Exp(numClass), Obs(numClass), auc(numClass), chi_stat(numClass), p(numClass);
	Exp.zeros(); Obs.zeros(); auc.zeros(); chi_stat.zeros();p.zeros();
	vec nzCnt(numClass);
	nzCnt.zeros();
	for (int ii = 0; ii < typeNum; ii++)
	{
	    uvec sel = arma::find(typeInx  ==  ii);
	    // cout <<"sel" << sel.t() << endl;
	    if(sel.size()> 0 ){

		urowvec scna1CurrType = scna1.cols(sel);
		urowvec scna2CurrType =  scna2.cols(sel);
		mat survivalCurrType = survival.rows(sel);
		urowvec inx = 3 * scna2CurrType + scna1CurrType + 1;
		// cout << survivalCurrType.t() << endl;
		mat temp = logRankMultiClass(survivalCurrType, inx.t(), classMap, numClass);
		// cout<< idx << "ii"<< ii << endl;
		for (int cc = 0; cc < numClass; cc++){
		    if(temp(cc,0) >=0) {
			nzCnt(cc) = nzCnt(cc) +1;
		    }else{
			temp.row(cc).fill(0);
		    }
		}
		Exp = Exp + temp.col(2);
		Obs = Obs + temp.col(1);
		auc = auc + temp.col(4);
	    }

	}
	// cout << 772 << endl;
	double chi_stat_agg = 0;
	for (int cc = 0; cc < numClass; cc++){
	    if(nzCnt(cc) > 0 ) {
		auc(cc)  = auc(cc)/nzCnt(cc); //AUC normalization
		chi_stat(cc)=((Exp(cc)-Obs(cc))*(Exp(cc)-Obs(cc)))/Exp(cc); 
		chi_stat_agg += chi_stat(cc);
		p(cc) = 1- R::pchisq(chi_stat(cc), 1, 1, 0);  /// degree of freedom 0
		if(Exp(cc) < Obs(cc)) p(cc) = -p(cc);
	    }else{
		p(cc) = -1000;
		auc(cc) = -1000;
	    }	
	}
	p_agg(idx) = 1- R::pchisq(chi_stat_agg, numClass -1, 1, 0);  /// degree of freedom 0
	// cout << 780 << endl;

	ExpMat.row(idx) = Exp.t();
	ObsMat.row(idx) = Obs.t();
	pMat.row(idx) = p.t();
	aucMat.row(idx)=auc.t();

    }
    Rcpp::List out; 
    out["Exp"] = ExpMat; out["Obs"] = ObsMat; out["p"] = pMat; out["auc"] = aucMat;
    out["p_agg"] = p_agg;
    return(out);
}
