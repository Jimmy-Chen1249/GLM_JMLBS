#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @keywords internal
 // [[Rcpp::export]]
 List gammaUpdate(const Rcpp::List& b_, const Rcpp::List& z_, const Rcpp::List& w_,
                  const Rcpp::List& pb_, const arma::vec& haz, const Rcpp::List& v_,
                  const Rcpp::List& u_,const Rcpp::List& h_, const int& K, const int& q,
                  const int& qq, const int& nev, const arma::vec& jcount) {
   
   // Newton-Raphson updates of gamma (E-step and M-step) using an exact observed
   // information calculation
   
   // declare score, E[delta x v*], and information matrix
   arma::mat Si = arma::zeros<arma::mat>(q+qq+K, w_.size()); //Score vector
   arma::vec S = arma::zeros<arma::vec>(q+qq+K);
   arma::mat Evstari = arma::zeros<arma::mat>(q+qq+K, w_.size());
   arma::mat I = arma::zeros<arma::mat>(q+qq+K, q+qq+K); //Information matrix
   arma::mat Gammaj = arma::zeros<arma::mat>(q+qq+K, nev);
   arma::mat Gamma = arma::zeros<arma::mat>(q+qq+K, q+qq+K);
   
   // loop over subjects
   for (int i=0; i<w_.size(); i++) {
     
     Rcpp::checkUserInterrupt();
     
     // extract matrices from lists for subject i
     arma::mat b = Rcpp::as<arma::mat>(b_[i]);
     arma::mat z = Rcpp::as<arma::mat>(z_[i]);
     arma::mat w = Rcpp::as<arma::mat>(w_[i]);
     arma::vec pb = Rcpp::as<arma::vec>(pb_[i]);
     arma::vec v = Rcpp::as<arma::vec>(v_[i]);
     arma::vec u = Rcpp::as<arma::vec>(u_[i]); //this is a column vector
     Rcpp::DataFrame h = Rcpp::as<Rcpp::DataFrame>(h_[i]);
     
     // subjects who are censored before the first failure time
     // do not contribute towards \gamma estimation
     int tj_ind = h["tj.ind"];
     if (tj_ind == 0) continue;
     
     int nj = w.n_cols;  // number of failure times upto time T_i
     int delta = h["delta"]; // delta_i
     
     arma::mat Ii_int = arma::zeros<arma::mat>(q+qq+K, q+qq+K); // information matrix (uninitialized) for subject i
     arma::mat bzt = b * z; // b x t(z)
     arma::mat bztev = bzt % repmat(w, 1, K); 
     // b x t(Z) . exp{W_1 %*% gamma_{11} + Z^*_11(t) %*% gamma_{12} + Z^*_11(t) %*% b * gamma_{12}} 
     arma::mat uev = repmat(trans(u), b.n_rows, 1) % repmat(w, 1, K); //1*nj
     // Z^*_{11}(t) .* exp{W_1 %*% gamma_{11} + Z^*_11(t) %*% gamma_{12} + Z^*_11(t) %*% b * gamma_{12}}
     arma::mat Eexpvj = (mean(w.each_col() % pb, 0)) % trans(haz.subvec(0, nj-1)); // this is a row vector
     arma::mat uEexpvj = Eexpvj % trans(u); //E lam_l(t) exp{} Z^*_11(t), this is row vector
     arma::mat u2Eexpvj = Eexpvj % trans(u) % trans(u); //E lam_l(t) exp{} Z^*_11(t) ^ 2
     arma::mat Eexpv = sum(Eexpvj, 1); // lambda0 x E[exp(v*gamma)]
     arma::mat uEexpv = sum(uEexpvj, 1); //sum(E lam_l(t) exp{} Z^*_11(t))
     arma::mat u2Eexpv = sum(u2Eexpvj, 1); //sum(E lam_l(t) exp{} Z^*_11(t) ^ 2)
     arma::mat hexpand = trans(repmat(haz.subvec(0, nj-1), K, 1)); // K reps of lambda0(tj)
     arma::mat outj = (mean(bztev.each_col() % pb, 0)) % hexpand;
     arma::mat uEbexpbj = (mean(bztev.each_col() % pb, 0)) % trans(u) % hexpand; //cross-prod
     arma::mat bzt2ev = bzt % bztev; // [b x t(z)]^2 . exp(v*gamma)
     arma::mat Ii_int_Kdiag = (mean(bzt2ev.each_col() % pb, 0)) % hexpand;
     
     arma::mat Eb = mean(b.each_col() % pb, 0);
     
     // loop of K longitudinal outcomes
     for (int k=0; k<K; k++) {
       // E[delta x v*(T_i)]
       Evstari(q+qq+k, i) = delta * arma::dot(z.col(nj*(k+1)-1), Eb);
       
       // score elements for K Zb's
       Si(q+qq+k, i) = arma::as_scalar(sum(outj.cols(nj*k, nj*(k+1)-1), 1));
       
       // cross-prod (diagonal) elements for K Zb's only
       Ii_int(q+qq+k, q+qq+k) = arma::as_scalar(sum(Ii_int_Kdiag.cols(nj*k, nj*(k+1)-1), 1));
       
       // cross-prod (off-diagonal) elements for K Zb's only
       for (int k2=k+1; k2<K; k2++) {
         arma::mat bztcross = bztev.cols(nj*k, nj*(k+1)-1) % bzt.cols(nj*k2, nj*(k2+1)-1);
         Ii_int(q+qq+k, k2+qq+q) = arma::as_scalar(sum((mean(bztcross.each_col() % pb, 0)) % trans(haz.subvec(0, nj-1)), 1));
         Ii_int(k2+qq+q, q+qq+k) = Ii_int(q+qq+k, k2+qq+q);
       }
       
       // cross-prod elements for q V_i's and K Zb's
       if (q > 0) {
         for (int j=0; j<q; j++) {
           Ii_int(j, q+qq+k) = Si(q+qq+k, i) * v(j);
           Ii_int(q+qq+k, j) = Ii_int(j, q+qq+k);
         }
         for(int jj = q; jj < (q+qq); jj++)
         {
           Ii_int(jj, q+qq+k) = arma::as_scalar(sum(uEbexpbj.cols(nj*k, nj*(k+1)-1), 1));//qq=1
           Ii_int(q+qq+k, jj) = Ii_int(jj, q+qq+k);
         }
         
       }
     } // end loop over outcomes k
     
     // Gamma_j vectors (for cross-product later)
     Gammaj.submat(q+qq, 0, q+qq+K-1, nj-1) += delta * trans(reshape(trans(outj), nj, K));
     
     // elements for q V_i's (only when q > 0)
     if (q > 0) {
       // E[delta x v]
       Evstari.submat(0, i, q-1, i) = delta * v;
       
       Evstari.submat(q, i, q+qq-1, i) = delta * u.row(nj-1); //delta_i * Z^*_{11}(u_i)
       // score elements
       Si.submat(0, i, q-1, i) = v * Eexpv;  // Esum(exp{}W_1)
       Si.submat(q, i, q+qq-1, i) = uEexpv; //Esum(exp{}Z^*_{11}(t))
       // cross-prod elements
       Ii_int.submat(0, 0, q-1, q-1) = (v * trans(v)) * arma::as_scalar(Eexpv);
       for(int jj = 0; jj < q; jj++)
       {
         Ii_int.submat(jj, q, jj, q+qq-1) = uEexpv * v(jj);
       }
       Ii_int.submat(q, 0, q+qq-1, q-1) = Ii_int.submat(0, q, q-1, q+qq-1);
       Ii_int.submat(q, q, q+qq-1, q+qq-1) = u2Eexpv;
       
       // Gamma_j elements
       Gammaj.submat(0, 0, q-1, nj-1) += delta * v * Eexpvj;
       //Gammaj.submat(q, 0, q+qq-1, nj-1) += uEbexpbj;//有问题
       Gammaj.submat(q, 0, q+qq-1, nj-1) += delta * uEexpvj;
     }
     
     S += (Evstari.col(i) - delta * Si.col(i)); // NB: actual score is sum(Evstari - Si)
     I += delta * Ii_int;
     
   } // end loop over subjects i
   
   // lambda0 x Gamma_j sum term (minus from information matrix)
   for (int t=0; t<nev; t++) {
     //Rcpp::Rcout << jcount(t) << std::endl;
     Gamma += Gammaj.col(t) * trans(Gammaj.col(t)) / jcount(t);
   }
   
   return List::create(
     Named("gDelta")  = solve(I, S),
     Named("I")  = I,
     Named("S")  = S
   );
   
 }


//' @keywords internal
 // [[Rcpp::export]]
 arma::mat hazHat(const Rcpp::List& w_, const Rcpp::List& pb_, const arma::vec& nev) {
   
   // lambda0(t) for profile score function of beta
   
   arma::vec haz = arma::zeros<arma::vec>(nev.n_elem);
   
   // loop over subjects
   for (int i=0; i<w_.size(); i++) {
     
     // extract matrices from lists for subject i
     arma::mat w = Rcpp::as<arma::mat>(w_[i]);
     arma::vec pb = Rcpp::as<arma::vec>(pb_[i]);
     
     haz.subvec(0, w.n_cols-1) += arma::trans(mean(w.each_col() % pb, 0));
     
   }
   
   return(nev / haz);
   
 }

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @keywords internal
 // [[Rcpp::export]]
 arma::mat lambdaUpdate(const Rcpp::List& b_, const Rcpp::List& imat_,
                        const Rcpp::List& zt_, const Rcpp::List& pb_,
                        const Rcpp::List& v_, const Rcpp::List& u_,
                        const arma::mat& gam, const arma::vec& gam_vec, 
                        const int& q, const arma::vec& nev, const Rcpp::List& h_) {
   
   // Updates of lambda0 (E-step and M-step)
   
   arma::vec haz = arma::zeros<arma::vec>(nev.n_elem);
   
   // loop over subjects
   for (int i=0; i<b_.size(); i++) {
     
     // extract matrices from lists for subject i
     arma::mat b = Rcpp::as<arma::mat>(b_[i]);
     arma::mat I = Rcpp::as<arma::mat>(imat_[i]);
     arma::mat zt = Rcpp::as<arma::mat>(zt_[i]);
     arma::vec pb = Rcpp::as<arma::vec>(pb_[i]);
     arma::vec v = Rcpp::as<arma::vec>(v_[i]);
     arma::vec u = Rcpp::as<arma::vec>(u_[i]);
     Rcpp::DataFrame h = Rcpp::as<Rcpp::DataFrame>(h_[i]);
     
     // subjects who are censored before the first failure time
     // do not contribute towards \lambda estimation
     int tj_ind = h["tj.ind"];
     if (tj_ind == 0) continue;
     
     arma::mat expW_new = exp((b * gam) * trans(I * zt)) % repmat(trans(u), b.n_rows, 1); //nMC*tj_ind
     arma::mat EexpVstar = mean(expW_new.each_col() % pb, 0);
     if (q > 0) {
       EexpVstar *= arma::as_scalar(exp(v.t() * gam_vec.subvec(0, q-1)));
     }
     haz.subvec(0, EexpVstar.n_cols-1) += EexpVstar.t();
     
   }
   
   return(nev / haz);
   
 }


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @keywords internal
// [[Rcpp::export]]
List gammaUpdate_fixDelta(const Rcpp::List& b_,   // list: nMC x r
                          const Rcpp::List& z_,   // list: r x (nj*K)  (t(Z_sur) or similar)
                          const Rcpp::List& w_,   // list: nMC x (nj*K)  expv*(t_j) already includes exp(W^T gamma11)*exp(Z^* gamma12)*exp(Zb beta)
                          const Rcpp::List& pb_,  // list: nMC weights, mean=1 typically
                          const arma::vec& haz,   // length >= max(nj)
                          const Rcpp::List& v_,   // list: q baseline covs (vector length q)
                          const Rcpp::List& u_,   // list: time-dependent cov values at grid, length nj (or nj x 1)
                          const Rcpp::List& h_,   // list of DF with delta, tj.ind
                          const int& K,           // # longitudinal outcomes
                          const int& q,           // # baseline covs in Cox (e.g. W)
                          const int& qq,          // # time-dependent fixed covs (here you用qq=1通常)
                          const int& nev,         // # unique event times
                          const arma::vec& jcount // length nev, typically #events at each tj (or risk/event counts used by package)
) {

  const int P = q + qq + K; // total gamma dim

  arma::mat Si      = arma::zeros<arma::mat>(P, w_.size()); // integral part score contributions (no delta)
  arma::mat Evstari = arma::zeros<arma::mat>(P, w_.size()); // event part contributions (has delta)
  arma::vec S       = arma::zeros<arma::vec>(P);

  arma::mat I       = arma::zeros<arma::mat>(P, P);         // observed info approx (no delta)
  arma::mat Gammaj  = arma::zeros<arma::mat>(P, nev);
  arma::mat Gamma   = arma::zeros<arma::mat>(P, P);

  // loop over subjects
  for (int i = 0; i < w_.size(); i++) {

    Rcpp::checkUserInterrupt();

    arma::mat b  = Rcpp::as<arma::mat>(b_[i]);      // nMC x r
    arma::mat z  = Rcpp::as<arma::mat>(z_[i]);      // r x (nj*K)
    arma::mat w  = Rcpp::as<arma::mat>(w_[i]);      // nMC x (nj*K)
    arma::vec pb = Rcpp::as<arma::vec>(pb_[i]);     // nMC x 1
    arma::vec v  = Rcpp::as<arma::vec>(v_[i]);      // q x 1

    // u: allow either vec or 1-col mat
    arma::vec u;
    {
      SEXP uu = u_[i];
      if (Rf_isMatrix(uu)) {
        arma::mat um = Rcpp::as<arma::mat>(uu);
        u = um.col(0);
      } else {
        u = Rcpp::as<arma::vec>(uu);
      }
    }

    Rcpp::DataFrame h = Rcpp::as<Rcpp::DataFrame>(h_[i]);

    int tj_ind = h["tj.ind"];
    if (tj_ind == 0) continue; // censored before first failure time

    int delta = h["delta"];
    int nj = w.n_cols / K;     // number of failure-time grid points up to T_i  (per outcome block)

    // safety
    if (nj <= 0) continue;

    // ---------- core quantities ----------
    // bzt = b * z  => nMC x (nj*K)
    arma::mat bzt = b * z;

    // bztev = (bzt % w)  elementwise => nMC x (nj*K)
    arma::mat bztev = bzt % w;

    // Eexpvj for each grid time (length nj): haz_j * E[ exp(v*) | O ] at t_j
    // first compute mean(w each_col % pb,0) => rowvec length (nj*K)
    arma::rowvec Ew_all = arma::mean(w.each_col() % pb, 0); // rowvec length nj*K

    // but hazard only length nj (same for each outcome block), we need block-wise treatment
    // for baseline/int parts (W, u) we only need the "common" nj grid:
    // take the first block as representative (because w already includes all parts, and haz is same)
    arma::rowvec Ew = Ew_all.cols(0, nj - 1); // rowvec length nj

    arma::rowvec hazj = arma::trans(haz.subvec(0, nj - 1)); // rowvec length nj
    arma::rowvec Eexpvj = Ew % hazj;                        // rowvec length nj
    double Eexpv = arma::as_scalar(arma::sum(Eexpvj));      // scalar
    double uEexpv = arma::as_scalar(arma::sum(Eexpvj % arma::trans(u.subvec(0, nj - 1))));
    double u2Eexpv = arma::as_scalar(arma::sum(Eexpvj % arma::square(arma::trans(u.subvec(0, nj - 1)))));

    // ---------- build Ii_int (information) ----------
    arma::mat Ii_int = arma::zeros<arma::mat>(P, P);

    // baseline v part: q x q
    if (q > 0) {
      Ii_int.submat(0, 0, q - 1, q - 1) = (v * v.t()) * Eexpv;
    }

    // baseline v vs u (q x qq) and u vs v
    if (q > 0 && qq > 0) {
      for (int j = 0; j < q; j++) {
        Ii_int(j, q) = uEexpv * v(j); // qq assumed 1, if qq>1 you should generalize here
        Ii_int(q, j) = Ii_int(j, q);
      }
    }

    // u-u
    if (qq > 0) {
      Ii_int(q, q) = u2Eexpv; // qq assumed 1
    }

    // ---------- now K blocks for random-effect association parameters ----------
    // outj (for each k) = hazj * E[ bzt % exp(v*) | O ] at each t_j
    // We need, for each k, columns (nj*k ... nj*(k+1)-1)
    arma::mat hexpand = arma::repmat(hazj, K, 1).t(); // nj x K then transposed => nj x K? (we want nj*K as rowvec mask)
    // easier: make a rowvec of hazj repeated K times
    arma::rowvec haz_rep = arma::repmat(hazj, 1, K); // 1 x (nj*K)

    // mean(bztev each_col % pb,0) => rowvec length nj*K
    arma::rowvec Ebztev = arma::mean(bztev.each_col() % pb, 0);

    // outj rowvec length nj*K
    arma::rowvec outj = Ebztev % haz_rep;

    // For K block info: need mean( (bzt % bztev) each_col % pb,0 ) % haz_rep
    arma::mat bzt2ev = (bzt % bztev); // IMPORTANT: materialize to avoid eGlue.each_col() error
    arma::rowvec Ebzt2ev = arma::mean(bzt2ev.each_col() % pb, 0);
    arma::rowvec Ii_int_Kdiag_all = Ebzt2ev % haz_rep;

    // Eb (mean of b with weights pb)
    arma::rowvec Eb = arma::mean(b.each_col() % pb, 0); // 1 x r

    // loop over outcomes
    for (int k = 0; k < K; k++) {

      int c1 = nj * k;
      int c2 = nj * (k + 1) - 1;

      // ---- event part: delta * (Zb)(T_i) ----
      // z is r x (nj*K), take last time point in this block: col = c2
      Evstari(q + qq + k, i) = delta * arma::as_scalar( z.col(c2).t() * Eb.t() );

      // ---- integral part score: sum_j outj_j (no delta) ----
      Si(q + qq + k, i) = arma::sum(outj.subvec(c1, c2));

      // ---- info diagonal for this k: sum_j Ii_int_Kdiag_all ----
      Ii_int(q + qq + k, q + qq + k) = arma::sum(Ii_int_Kdiag_all.subvec(c1, c2));

      // ---- info off-diagonal between k and k2 ----
      for (int k2 = k + 1; k2 < K; k2++) {
        int d1 = nj * k2;
        int d2 = nj * (k2 + 1) - 1;

        arma::mat cross = (bztev.cols(c1, c2) % bzt.cols(d1, d2)); // nMC x nj
        arma::rowvec Ecross = arma::mean(cross.each_col() % pb, 0); // 1 x nj
        double val = arma::as_scalar( arma::sum(Ecross % hazj) );

        Ii_int(q + qq + k, q + qq + k2) = val;
        Ii_int(q + qq + k2, q + qq + k) = val;
      }

      // ---- cross terms with baseline v and u ----
      if (q > 0) {
        for (int j = 0; j < q; j++) {
          Ii_int(j, q + qq + k) = Si(q + qq + k, i) * v(j);
          Ii_int(q + qq + k, j) = Ii_int(j, q + qq + k);
        }
      }
      if (qq > 0) {
        // cross(u , k): sum_j outj_j * u_j  (no delta)
        arma::rowvec out_block = outj.subvec(c1, c2);
        arma::rowvec urow = arma::trans(u.subvec(0, nj - 1));
        double uk = arma::as_scalar( arma::sum(out_block % urow) );

        Ii_int(q, q + qq + k) = uk;
        Ii_int(q + qq + k, q) = uk;
      }

      // ---- Gammaj: used for Gamma correction (no delta, consistent with package) ----
      // add vector for each event time t (j=1..nj):
      // Gammaj(q+qq+k, t) += outj_{k}(t)
      // Here outj_{k}(t) = outj[c1 + t]
      for (int t = 0; t < nj; t++) {
        // t index corresponds to event grid index
        Gammaj(q + qq + k, t) += outj(c1 + t);
      }
    }

    // baseline parts (v and u) Gammaj and scores
    if (q > 0) {
      // event part
      Evstari.submat(0, i, q - 1, i) = delta * v;

      // integral score part (no delta)
      Si.submat(0, i, q - 1, i) = v * Eexpv;

      // Gammaj baseline: for each grid time t, add v * Eexpvj(t)  (no delta)
      for (int t = 0; t < nj; t++) {
        Gammaj.submat(0, t, q - 1, t) += v * Eexpvj(t);
      }
    }

    if (qq > 0) {
      // event part: delta * u(T_i)  (qq assumed 1)
      Evstari(q, i) = delta * u(nj - 1);

      // integral score: sum u(t)*Eexpvj(t) (no delta)
      Si(q, i) = uEexpv;

      // Gammaj for u: add u(t)*Eexpvj(t) (no delta)
      for (int t = 0; t < nj; t++) {
        Gammaj(q, t) += u(t) * Eexpvj(t);
      }
    }

    // final subject contribution to S and I
    S += (Evstari.col(i) - Si.col(i));
    I += Ii_int;
  }

  // build Gamma (correction term), then Newton step is solve(I - Gamma, S)
  for (int t = 0; t < nev; t++) {
    if (jcount(t) > 0) {
      Gamma += (Gammaj.col(t) * Gammaj.col(t).t()) / jcount(t);
    }
  }

  arma::vec gDelta = arma::solve(I - Gamma, S);

  return List::create(
    Named("gDelta") = gDelta,
    Named("I")      = I,
    Named("Gamma")  = Gamma,
    Named("S")      = S
  );
}
