// [[Rcpp::depends("RcppArmadillo")]]

#include <RcppArmadillo.h>
#include <Rmath.h>


using namespace Rcpp;
using namespace arma;

RcppExport SEXP prep(SEXP EntrySEXP,
		     SEXP ExitSEXP,
		     SEXP StatusSEXP,
		     SEXP WeightSEXP,
		     SEXP XSEXP,
		     SEXP TruncationSEXP) {
  BEGIN_RCPP
    arma::vec entry = Rcpp::as<arma::vec>(EntrySEXP);
  arma::vec exit = Rcpp::as<arma::vec>(ExitSEXP);
  arma::Col<int> status = Rcpp::as<arma::Col<int> >(StatusSEXP);
  arma::Col<int>  weight = Rcpp::as<arma::Col<int> >(WeightSEXP);
  arma::mat x=Rcpp::as<arma::mat>(XSEXP);
  
  //Rcout << "x=" << x <<std::endl;
  
  bool Truncation=Rcpp::as<bool>(TruncationSEXP);
  
  unsigned n =exit.n_elem;
  if (Truncation) n *=2;
    
  mat XX(n,x.n_cols*x.n_cols);
  for (unsigned i=0; i<x.n_rows; i++) {
    rowvec Xi =x.row(i);
    XX.row(i) = vectorise(Xi.t()*Xi,1);
    if (Truncation) XX.row(i+n/2)=XX.row(i);
  }
  
  arma::Col<int> Sign;
  if (Truncation) {
    exit.insert_rows(0,entry);
    x.insert_rows(0,x);
    status.insert_rows(0,status);
    weight.insert_rows(0,weight);
    Sign.reshape(n,1); Sign.fill(1);
    for (unsigned i=0; i<(n/2);i++) Sign(i)= -1;
    status =status%(1+Sign);
   // weight =weight%(1+Sign);
  }
  
  //Rcout << "Status=" << status <<std::endl;
  arma::uvec idx0 = sort_index(status,1);
  //vec checktime = exit.elem(idx0);
  arma::uvec idx = stable_sort_index(exit.elem(idx0),0);
  idx = idx0.elem(idx);
  
  // Rcout << "idx" << idx <<std::endl;
  
  if (Truncation) {
    Sign = Sign.elem(idx);
  }
  if (x.n_rows>0) {
    XX = XX.rows(idx);
    x = x.rows(idx);
  }
  
  status = status.elem(idx);
  weight = weight.elem(idx);
  
  // Rcout << "Status=" << status <<std::endl;
  // Rcout << "Weight=" << weight <<std::endl;
  arma::uvec jumps = find(status>0);
  
  // Rcout << "jumps" << jumps <<std::endl;
  
  return(Rcpp::wrap(Rcpp::List::create(Named("XX")=XX,
				       Named("X")=x,
				       Named("jumps")=jumps,
				       Named("Sign")=Sign,
				       Named("ord")=idx,
				       Named("time")=exit,
				       Named("weight")=weight
				       )));
  END_RCPP
    }

colvec revcumsum(const colvec &a){
  unsigned n = a.n_rows;
  colvec res = a; double prev=0;
  for (unsigned i=0; i<n; i++){
    prev += a(n-i-1);
    res(n-i-1) = prev;
  }
  return(res);
}
colvec revcumsum(const colvec &a, const colvec &v1, const colvec &v2) {
  return (revcumsum(a%v1)/v2);
}


RcppExport SEXP PL(SEXP BetaSEXP, 
		   SEXP XSEXP, 
		   SEXP XXSEXP,
		   SEXP SignSEXP,
		   SEXP JumpsSEXP, 
		   SEXP WeightSEXP){
  
  BEGIN_RCPP
    
    colvec beta = Rcpp::as<colvec>(BetaSEXP);
  mat x = Rcpp::as<mat>(XSEXP);
  mat xx = Rcpp::as<mat>(XXSEXP);
  
  arma::uvec jumps = Rcpp::as<arma::uvec>(JumpsSEXP);
  arma::Col<int> Sign= Rcpp::as<arma::Col<int> >(SignSEXP);
  
  arma::Col<int> weight = Rcpp::as<arma::Col<int> >(WeightSEXP);
// Rcout <<  "weight=" <<weight <<std::endl;
  
  unsigned p=x.n_cols;
  
  colvec Xb = x*beta;
  //Rcout <<  Xb <<std::endl;
  colvec eXb = exp(Xb);
  //Rcout << eXb <<std::endl;
  if (Sign.n_rows==eXb.n_rows){
    eXb=Sign%eXb;
  }
  
  colvec S0=revcumsum(eXb);
  //Rcout << S0 <<std::endl;  
  
  //mat Wx(x.n_rows,p);
  //for (unsigned j=0;j<p;j++) {
  //  Wx.col(j)= weight%x.col(j); //weighted X
  //}
  
  // mat S1(x.n_rows,p);
  // for (unsigned j=0;j<p;j++) {
  //   S1.col(j)= revcumsum(x.col(j)%eXb);
  // }
  
  
  mat E(x.n_rows,p);
  for (unsigned j=0;j<p;j++) {
    E.col(j)= revcumsum(x.col(j),eXb,S0); // S1/S0(s)
  }
    
  // mat S2(x.n_rows,p);
  //for (unsigned j=0;j<p;j++) {
  // S2.col(j)= revcumsum(square(x.col(j))%eXb);
  //}
  //Rcout << "S2" <<S2<<std::endl; 
  
  
  mat xx2=xx;
  for (unsigned j=0; j<xx2.n_cols; j++){ // int(S2/S0 (s))
    xx2.col(j)=revcumsum(xx2.col(j),eXb,S0);
  }
  //Rcout << "xx2 =" << xx2 <<std::endl;
   
  xx2 = xx2.rows(jumps);  
 
  
  E=E.rows(jumps);
  S0=S0.elem(jumps);
  // Rcout <<  S0 <<std::endl;
  weight=weight.elem(jumps);
  // Rcout<< weight << std::endl;
  
  //arma::Mat<int> matW=arma::diagmat(weight);
  
  //Rcout << "x =" << x <<std::endl;
  //Rcout << "xjumps =" << x.rows(jumps) <<std::endl;

  //Rcout << "xx =" << xx.rows(jumps) <<std::endl;
  //Rcout << "xx2(j)=" << xx2 <<std::endl;
  //Rcout << "sum(xx2) =" << sum(xx2) <<std::endl;
  //Rcout << "reshape =" << reshape(sum(xx2),p,p) <<std::endl;
  mat xjumps=x.rows(jumps);

  mat grad = xjumps;
  for (unsigned j=0; j<grad.n_cols; j++){
    grad.col(j)=weight%(grad.col(j)-E.col(j));
  }
  
   // Rcout << "grad.old=" << matW*(x.rows(jumps)-E) <<std::endl;
   // Rcout << "grad.new=" << grad <<std::endl;
  
   mat phess=xx2;
   for (unsigned j=0; j<xx2.n_cols; j++){
     phess.col(j)=weight%xx2.col(j);
   }

 
//    Rcout << "phess =" << phess <<std::endl;
   
  mat WE=E;
  for (unsigned j=0; j<E.n_cols; j++) {
    WE.col(j) = weight%E.col(j);
  }
   
   mat hess = -(reshape(sum(phess),p,p)-WE.t()*E); 
   
   // Rcout << "WE =" << WE <<std::endl;   
   //Rcout<< "xjumps="<< xjumps <<std::endl;
   
   return(Rcpp::List::create(Rcpp::Named("jumps")=jumps,
			     Rcpp::Named("U")=grad,
			     Rcpp::Named("gradient")=sum(grad),
			     Rcpp::Named("hessian")=hess,
			     Rcpp::Named("S2S0")=xx,
			     Rcpp::Named("E")=E,
			     Rcpp::Named("S0")=S0,
			     Rcpp::Named("weight")=weight,
			     Rcpp::Named("xjumps")=xjumps
			     //   ,
			     // Named("S1")=S1,
			     // Named("S2")=S2
			     ));
   END_RCPP
     }

/* cumsumstrata - maybeuseful in the plot */

RcppExport SEXP cumsumstrataR(SEXP ia,
                              SEXP istrata,
                              SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;
  
  colvec tmpsum(nstrata); 
  //  tmpsum=tmpsum*0; 
  tmpsum.zeros(); 
  colvec res = a; 
  for (unsigned i=0; i<n; i++) {
    int ss=intstrata(i); 
    tmpsum(ss) += a(i); 
    res(i) = tmpsum(ss);
  }  
  
  List rres; 
  rres["res"]=res; 
  return(rres);
} /*}}}*/

colvec  cumsumstrata(colvec a,IntegerVector strata,int nstrata) {/*{{{*/
  unsigned n = a.n_rows;
  colvec tmpsum(nstrata); 
  tmpsum.zeros(); tmpsum.zeros(); 
  colvec res = a; 
  
  for (unsigned i=0; i<n; i++) {
    int ss=strata(i); 
    tmpsum(ss) += a(i); 
    res(i) = tmpsum(ss);
  }  
  
  return(res);
} /*}}}*/
