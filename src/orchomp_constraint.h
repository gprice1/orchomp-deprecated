
#ifndef _ORCHOMP_CONSTRAINT_H_
#define _ORCHOMP_CONSTRAINT_H_

#include "chomp-multigrid/chomp/Chomp.h"
#include "chomp-multigrid/chomp/Constraint.h"
#include "chomp-multigrid/chomp/ConstraintFactory.h"

namespace orchomp
{

class mod;

class UnifiedConstraint : public chomp::Constraint{
  public:
    
    std::vector< chomp::Constraint * > constraints;
    
    void addConstraint( chomp::Constraint * c ){
        constraints.push_back ( c );
    }
    virtual void evaluateConstraints(const chomp::MatX& qt, 
                                     chomp::MatX& h, 
                                     chomp::MatX& H);

    virtual size_t numOutputs(){ return 1; }
    
    //empty constructor
    UnifiedConstraint(){}
    ~UnifiedConstraint(){}

};

class ORTSRConstraint : public chomp::TSRConstraint {
  public:
    mod * module;

    ORTSRConstraint( mod * module, chomp::Transform & pose_0_w, 
                     chomp::MatX & Bw, chomp::Transform & pose_w_e ) :
         TSRConstraint( pose_0_w, Bw, pose_w_e ),
         module( module ) {}

    //this function takes in a robot state, qt, and returns the position of
    // the relevant end-effector in the world frame. This is equivalent
    //  to the transformation from the end-effector frame to the world
    //  frame.
    virtual void forwardKinematics( const chomp::MatX& qt,
                                    chomp::Transform & pos );
    
    //this takes in a state and it gets the jacobian 
    virtual void computeJacobian( const chomp::MatX& qt,
                                  const chomp::Transform & pose_world_ee,
                                  chomp::MatX & jacobian
                                  );
};



class ORJointLimitConstraint : public chomp::Constraint {
  public:
    mod * module;
    int n_outputs;

    ORJointLimitConstraint( mod * module)
        : module( module ) , n_outputs(1) {}
    virtual void evaluateConstraints(const chomp::MatX& qt, 
                                     chomp::MatX& h, 
                                     chomp::MatX& H);
    virtual size_t numOutputs(){
        return n_outputs;
    }
};


//this does nothing right now.
class ORConstraintFactory : public chomp::ConstraintFactory {
  public: 
    typedef std::pair< double, double > pair_d;
    mod * module;
    std::vector< chomp::Constraint * > constraints;
    std::vector< pair_d > times;

    ORConstraintFactory( mod * module );
    ~ORConstraintFactory();

    virtual chomp::Constraint* getConstraint(size_t t, size_t total);
    
    void addConstraint( chomp::Constraint * c, double start, double end );

    virtual void evaluate( const std::vector<chomp::Constraint*>& constraints, 
                   const chomp::MatX& xi, chomp::MatX& h_tot,
                   chomp::MatX& H_tot, int step);

};

}//namespace

#endif
