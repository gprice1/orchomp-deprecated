/*
* Copyright (c) 2008-2014, Matt Zucker
*
* This file is provided under the following "BSD-style" License:
*
* Redistribution and use in source and binary forms, with or
* without modification, are permitted provided that the following
* conditions are met:
*
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above
* copyright notice, this list of conditions and the following
* disclaimer in the documentation and/or other materials provided
* with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
* USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/

#include "Chomp.h"
#include "ConstraintFactory.h"
#include "Constraint.h"
#include "HMC.h"
#include <float.h>
#include <cmath>
#include <iomanip>

#define debug if (0) std::cout
#define debug_assert if (0) assert

namespace chomp {

  bool isnan( MatX & mat ){
      for ( int i = 0; i < mat.rows() ; i ++ ){
          for ( int j = 0; j < mat.cols(); j ++ ){
              if ( mat(i,j) != mat(i,j) ){
                  return true;
              }
          }
      }
      return false;
  }

  const char* eventTypeString(int eventtype) {
    switch (eventtype) {
    case CHOMP_INIT: return "CHOMP_INIT";
    case CHOMP_GLOBAL_ITER: return "CHOMP_GLOBAL_ITER";
    case CHOMP_LOCAL_ITER: return "CHOMP_LOCAL_ITER";
    case CHOMP_FINISH: return "CHOMP_FINISH";
    case CHOMP_TIMEOUT: return "CHOMP_TIMEOUT";
    case CHOMP_GOALSET_ITER: return "CHOMP_GOALSET_ITER";
    default: return "[INVALID]";
    }
  }

  ChompObserver::~ChompObserver() {}

  int ChompObserver::notify(const Chomp&, 
                            ChompEventType,
                            size_t,
                            double, double, double) { 
    return 0; 
  }

  DebugChompObserver::~DebugChompObserver() {}

  int DebugChompObserver::notify(const Chomp& c, 
                                 ChompEventType e,
                                 size_t iter,
                                 double curObjective, 
                                 double lastObjective,
                                 double constraintViolation) { 
    std::cout << "chomp debug: "
              << "event=" << eventTypeString(e) << ", "
              << "iter=" << iter << ", "
              << "cur=" << std::setprecision(10) << curObjective << ", "
              << "last=" << std::setprecision(10) << lastObjective << ", "
              << "rel=" << std::setprecision(10) << ((lastObjective-curObjective)/curObjective) << ", "
              << "constraint=" << std::setprecision(10) << constraintViolation << "\n";

    if (std::isnan(curObjective) || std::isinf(curObjective) ||
        std::isnan(lastObjective) || std::isinf(lastObjective)) {
      return 1;
    }
    return 0;
  }

  ChompGradientHelper::~ChompGradientHelper() {}

  ChompCollisionHelper::ChompCollisionHelper(size_t nc, size_t nw, size_t nb):
    ncspace(nc), nwkspace(nw), nbodies(nb) {}

  ChompCollisionHelper::~ChompCollisionHelper() {}

  ChompCollGradHelper::ChompCollGradHelper(ChompCollisionHelper* h,
                                           double g):
    chelper(h), gamma(g)
  {
    dx_dq = MatX(h->nwkspace, h->ncspace);
    cgrad = MatX(h->nwkspace, 1);
  }

  ChompCollGradHelper::~ChompCollGradHelper() {}

  double ChompCollGradHelper::addToGradient(const Chomp& c,
                                            MatX& g) {

    q1 = c.getTickBorderRepeat(-1).transpose();
    q2 = c.getTickBorderRepeat(0).transpose();

    double total = 0.0;

    for (int t=0; t<c.N; ++t) {

      q0 = q1;
      q1 = q2;
      q2 = c.getTickBorderRepeat(t+1).transpose();

      cspace_vel = 0.5 * (q2 - q0) * c.inv_dt;        
      cspace_accel = (q0 - 2.0*q1 + q2) * (c.inv_dt * c.inv_dt);

      for (size_t u=0; u < chelper->nbodies; ++u) {

        float cost = chelper->getCost(q1, u, dx_dq, cgrad);
        if (cost > 0.0) {

          wkspace_vel = dx_dq * cspace_vel;

          //this prevents nans from propagating. Several lines below, 
          //    wkspace_vel /= wv_norm if wv_norm is zero, nans propogate.
          if (wkspace_vel.isZero()){ continue; }

          wkspace_accel = dx_dq * cspace_accel;
          
          float wv_norm = wkspace_vel.norm();
          wkspace_vel /= wv_norm;

          // add to total
          double scl = wv_norm * gamma / c.inv_dt;

          total += cost * scl;
          
          P = MatX::Identity(chelper->nwkspace, chelper->nwkspace)
              - (wkspace_vel * wkspace_vel.transpose());

          K = (P * wkspace_accel) / (wv_norm * wv_norm);
         

          //          scalar * M-by-W        * (WxW * Wx1   - scalar * Wx1)
          g.row(t) += (scl * (dx_dq.transpose() *
                      (P * cgrad - cost * K)).transpose());
         

        }
      }
    }
    return total;

  }

  Chomp::Chomp(ConstraintFactory* f,
               const MatX& xi_init, // should be N-by-M
               const MatX& pinit, // q0
               const MatX& pgoal, // q1
               int nmax,
               double al,
               double obstol,
               size_t mg,
               size_t ml,
               double tt,
               double timeout_seconds,
               bool use_momentum):
    factory(f),
    observer(0),
    ghelper(0),
    objective_type(MINIMIZE_ACCELERATION),
    maxN(nmax),
    xi(xi_init),
    q0(pinit),
    q1(pgoal),
    alpha(al),
    objRelErrTol(obstol),
    min_global_iter(0),
    max_global_iter(mg),
    min_local_iter(0),
    max_local_iter(ml),
    full_global_at_final(false),
    t_total(tt),
    timeout_seconds( timeout_seconds ),
    didTimeout( false ),
    use_mutex( false ),
    use_goalset( false ),
    use_momentum( use_momentum ),
    hmc( NULL )
  {

    N = xi.rows();
    N_sub = 0;

    M = xi.cols();

    assert(q0.rows() >= 1 && q0.cols() == xi.cols());
    assert(q1.rows() >= 1 && q1.cols() == xi.cols());

    minN = N;
    assert(maxN >= minN);

  }
 
  void Chomp::initMutex(){
      use_mutex = true;
      pthread_mutex_init( &trajectory_mutex, NULL );
  }


  void Chomp::clearConstraints() {

    while (!constraints.empty()) {
      delete constraints.back();
      constraints.pop_back();
    }

  }
  
  //this sets up an iteration of 'normal' chomp, as opposed to
  //    goal set chomp.
  void Chomp::prepareStandardChomp(){
    skylineChol(N, coeffs, L);

    b.resize(N,M);
    b.setZero();

    c = createBMatrix(N, coeffs, q0, q1, b, dt);

    g.resize(N,M);
  }

  void Chomp::prepareChomp() {
    
    if (objective_type == MINIMIZE_VELOCITY) {

      coeffs.resize(1,2);
      coeffs_sub.resize(1,1);

      coeffs << -1, 2;
      coeffs_sub << 2;

    } else {
    
      coeffs.resize(1,3);
      coeffs_sub.resize(1,2);
      
      coeffs << 1, -4, 6;
      coeffs_sub << 1, 6;

    }
    
    dt = t_total / (N+1);
    inv_dt = (N+1) / t_total;

    if (objective_type == MINIMIZE_VELOCITY) {
      fscl = inv_dt*inv_dt;
    } else {
      fscl = inv_dt*inv_dt*inv_dt;
    }

    cur_global_iter = 0;
    cur_local_iter = 0;

    clearConstraints();

    if (factory) {
      factory->getAll(N, constraints);
    }
    

    if( use_goalset ){ prepareGoalSet(); }
    else { prepareStandardChomp(); }


    // decide whether base case or not
    bool subsample = ( N > minN && !use_goalset );
    if (full_global_at_final && N >= maxN) {
      subsample = false;
    }

    if (!subsample) {
      N_sub = 0;
      if ( use_momentum ){
        momentum.resize( N, M );
        momentum.setZero();
      }
      if( hmc ){ hmc->setupRun(); }

    } else {
      N_sub = (N+1)/2;
      g_sub.resize(N_sub, M);
      xi_sub.resize(N_sub, M);
      skylineChol(N_sub , coeffs_sub, L_sub); 
    }


  }
  
  void Chomp::updateGradient() { 

    Ax.resize(xi.rows(), xi.cols());

    if( use_goalset ){ diagMul(coeffs, goalset_coeffs, xi, Ax); }
    else { diagMul(coeffs, xi, Ax); }

    g = Ax + b;
    if (ghelper) {
        fextra = ghelper->addToGradient(*this, g);
    } else {
        fextra = 0;
    }
  }


  // precondition: prepareChomp was called for this resolution level
  void Chomp::prepareChompIter() {
    
    if (hmc && !N_sub ) {
      hmc->iteration( cur_global_iter, xi, momentum, L, lastObjective );
    }
    
    // compute gradient, constraint & Jacobian and subsample them if
    // needed

    // get the gradient
    updateGradient();

    if (N_sub) {
      // we want even rows of xi and g
      // q0 q1 q2 q3 q4 q5 q6  with n = 7
      // q0    q2    q4    q6  

      for (int t=0; t<N_sub; t++){

        g_sub.row(t) = g.row(2*t);
        xi_sub.row(t) = xi.row(2*t);

      }

      if (factory) {
        factory->evaluate(constraints, xi, h_sub, H_sub, 2);
      }

    } else {
      // should we instead grab relevant blocks from this for local
      // smoothing?
      if (factory) {
        factory->evaluate(constraints, xi, h, H);
      }

    }

    if (h.rows()) {
      hmag = h.lpNorm<Eigen::Infinity>();
      debug << "in prepareChompIter, ||h|| = " << hmag << "\n";
    } else if (h_sub.rows()) { 
      hmag = h_sub.lpNorm<Eigen::Infinity>();
      debug << "in prepareChompIter, ||h|| = " << hmag << "\n";
    } else {
      hmag = 0;
    }

  }

  // precondition: prepareChomp was called for this resolution level
  // updates chomp equation until convergence at the current level
  
  void Chomp::runChomp(bool global, bool local) {
    
    prepareChompIter();
    lastObjective = evaluateObjective();
    //std::cout << "initial objective is " << lastObjective << "\n";

    if (notify(CHOMP_INIT, 0, lastObjective, -1, hmag)) { 
      global = false;
      local = false;
    }
    
    cur_global_iter = 0;
    while (global) {
      chompGlobal();
      ++cur_global_iter;
      ++total_global_iter;
      prepareChompIter();
      double curObjective = evaluateObjective();
      if (cur_global_iter > min_global_iter && (
          goodEnough(lastObjective, curObjective) ||
          cur_global_iter >= max_global_iter)) {
        global = false;
      }

      //check for a timeout
      if (canTimeout && stop_time < TimeStamp::now() ){
          global = false;
          didTimeout = true;
          notify(CHOMP_TIMEOUT, cur_global_iter, 
                 curObjective, lastObjective, hmag);

      } else if (notify(CHOMP_GLOBAL_ITER, cur_global_iter, 
                 curObjective, lastObjective, hmag)) {
        global = false;
      }
      lastObjective = curObjective;
    }

    if (use_goalset){ finishGoalSet(); }

    cur_global_iter = 0;

    if (full_global_at_final && N >= maxN) {
      local = false;
    }

    cur_local_iter = 0;
    while (local) {
      localSmooth();
      ++cur_local_iter;
      ++total_local_iter;
      prepareChompIter();
      double curObjective = evaluateObjective();
      if ( cur_local_iter > min_local_iter && (
          goodEnough(lastObjective, curObjective) ||
          cur_local_iter >= max_local_iter)) {
        local = false;
      }

      //check for a timeout
      if (canTimeout && stop_time < TimeStamp::now() ){
          local = false;
          didTimeout = true;
          notify(CHOMP_TIMEOUT, cur_local_iter, 
                 curObjective, lastObjective, hmag);

      } else if (notify(CHOMP_LOCAL_ITER, cur_local_iter, 
                 curObjective, lastObjective, hmag)) {
        local = false;
      }
      lastObjective = curObjective;
    }
    
    std::cout << "\nDone with global chomp" << std::endl;

    if (factory && N_sub) {
      factory->evaluate(constraints, xi, h, H);
      if (h.rows()) {
        hmag = h.lpNorm<Eigen::Infinity>();
      }
    }

    //std::cout << "after chomp, objective is " << lastObjective << "\n";

    notify(CHOMP_FINISH, 0, lastObjective, -1, hmag);

  }

  // calls runChomp and upsamples until N big enough
  // precondition: N <= maxN
  // postcondition: N >= maxN
  void Chomp::solve(bool doGlobalSmoothing, bool doLocalSmoothing) {
   
    if ( timeout_seconds < 0 ){
        canTimeout = false;
    }else {
        canTimeout = true;
        stop_time = TimeStamp::now() +
                    Duration::fromDouble( timeout_seconds );
    }
    
    if ( hmc ){ 
        use_momentum = true;
        hmc->setupHMC( objective_type, alpha );
    }

    total_global_iter = 0;
    total_local_iter = 0;
    cur_global_iter = 0;
    cur_local_iter = 0;

    //std::cout << "initial trajectory has length " << N << "\n";

    while (1) {
      
      prepareChomp();

      runChomp(doGlobalSmoothing, doLocalSmoothing); 


      if (N >= maxN) { // have already upsampled enough
        break;
      } else {
        upsample();
        //std::cout << "upsampled to " << N << "\n";
      }
      
    }

    if (factory) {
      // evaluate full constraint at end

      if (N_sub) { 
        factory->evaluate(constraints, xi, h, H);
      }

      if (h.rows()) {
        //std::cout << "after solve, ||h||_inf = " << h.lpNorm<Eigen::Infinity>() << "\n";
      }
    }

  }

  MatX Chomp::getTickBorderRepeat(int tick) const {
    if (tick < 0) { 
      return getPos(q0, (tick+1)*dt);
    } else if (tick >= xi.rows()) {
      return getPos(q1, (tick-xi.rows())*dt);
    } else {
      return xi.row( tick );
    }
  }

  // upsamples the trajectory by 2x
  void Chomp::upsample() {

    int N_up = 2*N+1; // e.g. size 3 goes to size 7
    //MatX xi_up(M*N_up, 1), q_t(M,1), q_prev(M,1), q_next(M,1);
    
    MatX xi_up(N_up, M);

    // q0    d0    d1    d2    q1   with n = 3
    // q0 u0 u1 u2 u3 u4 u5 u6 q1   with n = 7
    //
    // u0 = 0.5*(q0 + d0)
    // u1 = d0
    // u2 = 0.5*(d0 + d1)
    // u3 = d1 
    // u4 = 0.5*(d1 + d2)
    // u5 = d2
    // u6 = 0.5*(d2 + q1)

    for (int t=0; t<N_up; ++t) { // t is timestep in new (upsampled) regime

      if (t % 2 == 0) {

        assert(t == N_up-1 || (t/2) < xi.rows());
        assert(t < xi_up.rows());

        if (objective_type == MINIMIZE_VELOCITY) {

          MatX qneg1 = getTickBorderRepeat(t/2-1);
          MatX qpos1 = getTickBorderRepeat(t/2);
          xi_up.row(t) = 0.5 * (qneg1 + qpos1);

        } else { 

          MatX qneg3 = getTickBorderRepeat(t/2-2);
          MatX qneg1 = getTickBorderRepeat(t/2-1);
          MatX qpos1 = getTickBorderRepeat(t/2);
          MatX qpos3 = getTickBorderRepeat(t/2+1);

          const double c3 = -1.0/160;
          const double c1 = 81.0/160;

          xi_up.row(t) = c3*qneg3 + c1*qneg1 + c1*qpos1 + c3*qpos3;

        }


      } else {
        xi_up.row(t) = xi.row(t/2);
      }

    }

    N = N_up;

    lockTrajectory();
    xi = xi_up;
    unlockTrajectory();

    L = L_sub = g = g_sub = h = h_sub = H = H_sub = 
      P = HP = Y = W = Ax = delta = MatX();

    N_sub = 0;

  }

  // single iteration of chomp
  void Chomp::chompGlobal() { 
    
    assert(xi.rows() == N && xi.cols() == M);
    assert(Ax.rows() == N && Ax.cols() == M);
    

    // see if we're in our base case (not subsampling)
    bool subsample = N_sub != 0;
    
    const MatX& H_which = subsample ? H_sub : H;
    const MatX& g_which = subsample ? g_sub : g;
    const MatX& L_which = subsample ? L_sub : L;
    const MatX& h_which = subsample ? h_sub : h;
    const int  N_which = subsample ? N_sub : N;
    

    if (H_which.rows() == 0) {

      if (subsample) {
        assert( g_which.rows() == N_sub );
      } else {
        assert( g_which.rows() == xi.rows() );
      }
      
      skylineCholSolve(L_which, g_which);
      
      lockTrajectory();
      if (subsample) {
        for (int t=0; t<N_sub; ++t) {
          xi.row(2*t) -= g_which.row(t) * alpha;
        }
      } else if ( use_momentum ) {
        momentum += g_which * alpha;
        xi -= momentum;
      } else {
        xi -= g_which * alpha;
      }
      unlockTrajectory();
    
    //chomp update with constraints
    } else {

      P = H_which.transpose();
      
      // TODO: see if we can make this more efficient?
      for (int i=0; i<P.cols(); i++){
        skylineCholSolveMulti(L_which, P.col(i));
      }

      debug << "H = \n" << H << "\n";
      debug << "P = \n" << P << "\n";
  
      HP = H_which*P;

      cholSolver.compute(HP);
      Y = cholSolver.solve(P.transpose());

      debug << "HP*Y = \n" << HP*Y << "\n";
      debug << "P.transpose() = \n" << P.transpose() << "\n";
      debug_assert( P.transpose().isApprox( HP*Y ) );

      int newsize = H_which.cols();
      assert(newsize == N_which * M);

      assert(g_which.rows() == N_which && g_which.cols() == M);
      

      Eigen::Map<const Eigen::MatrixXd> g_flat(g_which.data(),
                                                 newsize, 1);
      W = (MatX::Identity(newsize,newsize) - H_which.transpose() * Y)
          * g_flat * alpha;
      skylineCholSolveMulti(L_which, W);

      Y = cholSolver.solve(h_which);

      debug_assert( h_which.isApprox(HP*Y) );
      
      //handle momentum if we need to.
      if (!subsample && use_momentum){
        Eigen::Map<Eigen::MatrixXd> momentum_flat(
                            momentum.data(), newsize, 1);
        momentum_flat += W;
        delta = momentum_flat + P * Y;
      }else {
        delta = W + P * Y;
      }

      assert(delta.rows() == newsize && delta.cols() == 1);

      Eigen::Map<const Eigen::MatrixXd> delta_rect(delta.data(),
                                                   N_which, M);
      
      debug << "delta = \n" << delta << "\n";
      debug << "delta_rect = \n" << delta_rect << "\n";
      
      assert(delta_rect.rows() == N_which && delta_rect.cols() == M);

      lockTrajectory();
      if (subsample) {
        
        for (int t=0; t<N_sub; ++t) {
          xi.row(2*t) -= delta_rect.row(t);
        }
        
      } else {
        
        xi -= delta_rect;
        
      }

      unlockTrajectory();
    }
    

  }

  // single iteration of local smoothing
  //
  // precondition: prepareChompIter has been called since the last
  // time xi was modified
  void Chomp::localSmooth() {

    MatX h_t, H_t, P_t, P_t_inv, delta_t;

    hmag = 0;

    for (int t=0; t<N; ++t){

      Constraint* c = constraints.empty() ? 0 : constraints[t];
      
      bool is_constrained = (c && c->numOutputs() > 0);

      if (is_constrained) { 
        c->evaluateConstraints(xi.row(t), h_t, H_t);
        is_constrained = h_t.rows() > 0;
        if (is_constrained ){
          hmag = std::max(hmag, h_t.lpNorm<Eigen::Infinity>());
        }
      }

      if (is_constrained) {

        assert(size_t(H_t.rows()) == constraints[t]->numOutputs());
        assert(H_t.cols() == M);
        assert(size_t(h_t.rows()) == constraints[t]->numOutputs());
        assert(h_t.cols() == 1);
    
        // Now we calculate, using opencv

        P_t = H_t*H_t.transpose();
        P_t_inv = P_t.inverse();

        // transpose g to be a column vector
        delta_t = ( -alpha*(MatX::Identity(M,M)-H_t.transpose()*P_t_inv*H_t)*g.row(t).transpose()
                    -H_t.transpose()*P_t_inv*h_t );

      } else {

        // transpose g to be a column vector
        delta_t = -alpha * g.row(t).transpose();
        
      }
      
      lockTrajectory();
      // transpose delta to be a row vector
      xi.row(t) += delta_t.transpose();
      unlockTrajectory();

    }

  }

  // evaluates the objective function for cur. thing.
  //
  // only works if prepareChompIter has been called since last
  // modification of xi.
  double Chomp::evaluateObjective() {

    /*
  
      K is (n+1)-by-n
  
      K = [  1  0  0  0 ... 0  0  0
      -2  1  0  0 ... 0  0  0 
      1 -2  1  0 ... 0  0  0
      0  1 -2  1 ... 0  0  0
      ...
      0  0  0  0 ... 1 -2  1 
      0  0  0  0 ... 0  1 -2 
      0  0  0  0 ... 0  0  1 ]
  
      e = [ -x0  x0   0 ... 0  x1 -x1 ]^T
  
      ||Kx + e||^2 = 
       
    */
    
    if ( false ){
        const double xi_Ax = 0.5 * mydot( xi, Ax );
        const double xi_b = mydot( xi, b );
        const double smoothness = ( xi_Ax + xi_b + c ) * fscl;
        std::cout <<  "xi * Ax    = " << xi_Ax << std::endl;
        std::cout <<  "xi * b     = " << xi_b << std::endl;
        std::cout <<  "c          = " << c << std::endl;
        std::cout <<  "fextra     = " << fextra << std::endl;
        std::cout <<  "fscale     = " << fscl << std::endl;
        std::cout <<  "smoothness = " << smoothness << std::endl;

        return smoothness + fextra;
    }

    return (0.5 * mydot( xi, Ax ) + mydot( xi, b ) + c ) + fextra;

  }

  // returns true if performance has converged
  bool Chomp::goodEnough(double oldObjective, double newObjective) {
    
    return (fabs((oldObjective-newObjective)/newObjective)<objRelErrTol);

  }

  void Chomp::constrainedUpsampleTo(int Nmax, double htol, double hstep) {

    MatX h, H, delta;
  
    while (N < Nmax) { 

      upsample();
      prepareChomp();

      double hinit = 0, hfinal = 0;
      
      lockTrajectory();
      for (int i=0; i<N; i+=2) {

        Constraint* c = constraints.empty() ? 0 : constraints[i];

        if (!c || !c->numOutputs()) { continue; }
        
        for (int iter=0; ; ++iter) { 
          c->evaluateConstraints(xi.row(i), h, H);
          if (h.rows()) {
            double hn = h.lpNorm<Eigen::Infinity>();
            if (iter == 0) { hinit = std::max(hn, hinit); }
            if (hn < htol) { hfinal = std::max(hn, hfinal); break; }
            delta = H.colPivHouseholderQr().solve(h);
            xi.row(i) -= hstep * delta.transpose();
          }
        }
      
      }
      unlockTrajectory();

      prepareChompIter();
      double f = evaluateObjective();
      if (0) { f = f ? f : f; } // shut up
      //std::cout << "after upsample to " << N << ", objective is " << f << " and ||h|| was reduced from " << hinit << " to " << hfinal << "\n";

    }

  }

  int Chomp::notify(ChompEventType event,
                    size_t iter,
                    double curObjective,
                    double lastObjective,
                    double constraintViolation) const {

    if (observer) {

      return observer->notify(*this, event, iter, 
                              curObjective, lastObjective, constraintViolation);

    } else {

      return 0;

    }

  }

  ////////////////GOAL SET FUNCTIONS//////////////////////////////////

  void Chomp::useGoalSet( Constraint * goalset ){
      this->goalset = goalset;
      use_goalset = true;
  }


  void Chomp::prepareGoalSet(){
    
    //resize xi, and add q1 into it.
    xi.conservativeResize( xi.rows() + 1, xi.cols() );
    xi.row( xi.rows() - 1 ) = q1;
    int n = xi.rows();

    if (objective_type == MINIMIZE_VELOCITY) {
      goalset_coeffs.resize(1,1);
      goalset_coeffs << 1;

    } else {
      goalset_coeffs.resize(2,2);
      goalset_coeffs << 6, -3,
                       -3,  2 ;
    }
    
    constraints.push_back( goalset );

    skylineChol( n, coeffs, goalset_coeffs, L);

    b.resize(n,M);
    b.setZero();

    c = createBMatrix(n, coeffs, q0, b, dt);

    g.resize(n,M);

    //if we are doing goal set chomp, we should not subsample
    N_sub = 0;

    N = n;
  }

  void Chomp::finishGoalSet(){
    
    std::cout << "Finishing goal set" << std::endl;

    use_goalset = false;

    q1 = xi.row( xi.rows() - 1 );
    xi.conservativeResize( xi.rows() -1, xi.cols() );
    
    N = xi.rows();

    //remove the goal constraint, so that it is not deleted along
    //  with the other constraints.
    constraints.pop_back();
    
    std::cout << "Preparing chomp" <<std::endl;

    prepareChomp();
    std::cout << "Done Preparing chomp" <<std::endl;
  }
}// namespace
